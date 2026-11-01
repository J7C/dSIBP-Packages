// FlintNDE C++ 引擎实现：请求校验、部分分式极点递推、prefix 认证输运与批量线程池。
//
// 语义镜像 flintnde：
// - PartialFractionSystem.solution_taylor_coefficients 的极点状态递推（offset 取
//   pole-center、逐极点 (state-current)/offset、R*state 累加、除以 n+1）；
// - _transport_ordinary_segment 的逐段半径校验（|delta| < radius_fraction*最近极点
//   距离，中点 double 比较）与段后中点重启；
// - certified 双链与 direct_user_point_path 的“用户点即节点、不插点不删点”路线。
// 输运只保留 prefix 复用：两链同精度只差阶数，同一步内 reference 级数的前
// primary_order+1 阶系数就是 primary 级数（前缀性质），逐列算术与独立双链
// 完全相同。留数乘法内核按 --residues 选择：sparse 走 CSR 非零元跳乘，
// dense 走 acb_mat_mul+acb_mat_add（精确实零项不改变球，两路线结果一致），
// auto 按留数稀疏度自适应。
#include "engine.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <set>
#include <thread>
#include <vector>

#if defined(_WIN32)
#include <process.h>
#endif

#include "flint/flint.h"
#include "util.hpp"

namespace fnde {

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_seconds(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double>(end - start).count();
}

// ---------- JSON 构造小工具 ----------

fjson::Value jstr(std::string text) {
    fjson::Value v;
    v.type = fjson::Value::Type::String;
    v.text = std::move(text);
    return v;
}

fjson::Value jnum(double x) {
    fjson::Value v;
    v.type = fjson::Value::Type::Number;
    v.number = x;
    return v;
}

fjson::Value jbool(bool x) {
    fjson::Value v;
    v.type = fjson::Value::Type::Bool;
    v.boolean = x;
    return v;
}

fjson::Value jobj() {
    fjson::Value v;
    v.type = fjson::Value::Type::ObjectT;
    return v;
}

fjson::Value jarr() {
    fjson::Value v;
    v.type = fjson::Value::Type::ArrayT;
    return v;
}

// 与桥接 _acb_record 对应：中点字符串 + 完整十进制球。
fjson::Value acb_record(const acb_t value, long digits) {
    fjson::Value record = jobj();
    record.fields["real"] = jstr(arb_midpoint_string(acb_realref(value), digits));
    record.fields["imag"] = jstr(arb_midpoint_string(acb_imagref(value), digits));
    for (const char* name : {"real_ball", "imag_ball"}) {
        arb_srcptr part = (name[0] == 'r') ? acb_realref(value) : acb_imagref(value);
        std::string midpoint;
        std::string radius;
        long exponent = 0;
        arb_ball_fields(part, digits, midpoint, radius, exponent);
        fjson::Value ball = jobj();
        ball.fields["midpoint"] = jstr(midpoint);
        ball.fields["radius"] = jstr(radius);
        ball.fields["exponent"] = jnum(static_cast<double>(exponent));
        record.fields[name] = ball;
    }
    return record;
}

fjson::Value vector_records(const acb_mat_t vector, long digits) {
    fjson::Value items = jarr();
    for (long row = 0; row < acb_mat_nrows(vector); ++row) {
        items.items.push_back(acb_record(acb_mat_entry(vector, row, 0), digits));
    }
    return items;
}

std::string arb_ball_string(const arb_t value, long digits) {
    char* text = arb_get_str(value, digits, ARB_STR_MORE);
    if (text == nullptr) throw std::runtime_error("arb_get_str failed");
    std::string result(text);
    flint_free(text);
    return result;
}

fjson::Value refusal(
    const std::string& status, const std::string& message, const fjson::Value* request) {
    fjson::Value v = jobj();
    v.fields["status"] = jstr(status);
    if (request != nullptr) {
        const fjson::Value* schema = request->find("schema");
        if (schema != nullptr && schema->is_string()) v.fields["schema"] = *schema;
        const fjson::Value* digest = request->find("masterDigest");
        if (digest != nullptr && digest->is_string()) v.fields["masterDigest"] = *digest;
    }
    v.fields["message"] = jstr(message);
    return v;
}

void require_exact_keys(
    const fjson::Value& value, const std::set<std::string>& expected, const char* label) {
    if (!value.is_object()) {
        throw std::runtime_error(std::string(label) + " must be an object");
    }
    std::set<std::string> actual;
    for (const auto& entry : value.fields) actual.insert(entry.first);
    if (actual != expected) {
        std::string message = std::string(label) + " fields do not match the current schema";
        for (const auto& key : expected) {
            if (actual.find(key) == actual.end()) message += "; missing=" + key;
        }
        for (const auto& key : actual) {
            if (expected.find(key) == actual.end()) message += "; unexpected=" + key;
        }
        throw std::runtime_error(message);
    }
}

// ---------- 系统构造：letters -> 合并极点的部分分式系统 ----------

// 稠密留数的 CSR 稀疏视图：向量递推按非零元跳乘。
// 零元以精确零球参与稠密 acb_mat_mul 时本就不改变结果，跳过不损失数值。
struct SparseResidue {
    long dimension = 0;
    long nonzero_count = 0;
    std::vector<long> row_ptr;              // dimension + 1
    std::vector<long> cols;                 // nonzero_count
    std::vector<std::unique_ptr<Acb>> vals; // nonzero_count
    // state_rows：R_j 非零列的并集，即累加 R_j*state 时实际会读到的 state 行。
    // 稀疏内核只推进这些行；其余行的 state 值从不被读取，跳过推进不改变结果。
    std::vector<long> state_rows;
    bool empty() const { return nonzero_count == 0; }
};

// auto 模式的稀疏判定阈值：留数非零元占比 <= 该值时用稀疏内核。
constexpr double kSparseDensityThreshold = 0.25;

struct SegmentSystem {
    long dimension = 0;
    std::vector<std::unique_ptr<AcbMat>> residues;
    std::vector<std::unique_ptr<Acb>> poles;
    std::vector<SparseResidue> sparse_residues;
    bool use_sparse = true;        // 该段实际使用的留数内核
    double residue_density = 0.0;  // 非零元 / (极点数*d*d)，auto 判据与报告用
};

struct SegmentJob {
    std::unique_ptr<Acb> start;
    std::vector<std::unique_ptr<Acb>> points;
    std::vector<long> user_indices;
    SegmentSystem system;
};

std::unique_ptr<GaussQ> parse_gaussian_record(const fjson::Value& record, const char* label) {
    require_exact_keys(record, {"real", "imag"}, label);
    if (!record.at("real").is_string() || !record.at("imag").is_string()) {
        throw std::runtime_error(std::string(label) + " real/imag must be strings");
    }
    std::unique_ptr<GaussQ> value(new GaussQ());
    parse_gaussian_part(record.at("real").text, value->re);
    parse_gaussian_part(record.at("imag").text, value->im);
    return value;
}

bool gaussian_is_zero(const GaussQ& value) {
    return fmpq_is_zero(value.re) && fmpq_is_zero(value.im);
}

void parse_scalar_pair(const fjson::Value& record, acb_t out, const char* label) {
    require_exact_keys(record, {"real", "imag"}, label);
    if (!record.at("real").is_string() || !record.at("imag").is_string()) {
        throw std::runtime_error(std::string(label) + " real/imag must be strings");
    }
    parse_rational_or_decimal(record.at("real").text, acb_realref(out));
    parse_rational_or_decimal(record.at("imag").text, acb_imagref(out));
}

// Wolfram 标量文本（segments[].start，如 "0"）：镜像 _mma_scalar 规整后按实数解析。
void parse_scalar_text(const std::string& raw, acb_t out) {
    std::string cleaned;
    cleaned.reserve(raw.size());
    bool in_backtick_suffix = false;
    for (size_t index = 0; index < raw.size(); ++index) {
        const char character = raw[index];
        if (character == '`') {
            // 丢弃反引号精度后缀中的数字与点。
            in_backtick_suffix = true;
            continue;
        }
        if (in_backtick_suffix) {
            if (std::isdigit(static_cast<unsigned char>(character)) || character == '.') {
                continue;
            }
            in_backtick_suffix = false;
        }
        if (character == '*' && index + 1 < raw.size() && raw[index + 1] == '^') {
            cleaned.push_back('e');
            ++index;
            continue;
        }
        cleaned.push_back(character);
    }
    parse_rational_or_decimal(cleaned, acb_realref(out));
    arb_zero(acb_imagref(out));
}

// 镜像桥接 _letters_to_system：beta==0 跳过；相同极点合并留数。
std::vector<SparseResidue> build_sparse_residues(const SegmentSystem& system);

void build_segment_system(
    const fjson::Value& letters, long dimension, ResidueMode residue_mode,
    const std::string& name, SegmentSystem& out) {
    out.dimension = dimension;
    out.residues.clear();
    out.poles.clear();
    std::vector<std::unique_ptr<GaussQ>> pole_exact;
    for (size_t index = 0; index < letters.items.size(); ++index) {
        const fjson::Value& record = letters.items[index];
        require_exact_keys(record, {"alpha", "beta", "residue"}, "letters[]");
        std::unique_ptr<GaussQ> beta = parse_gaussian_record(record.at("beta"), "letters[].beta");
        if (gaussian_is_zero(*beta)) continue;
        std::unique_ptr<GaussQ> alpha =
            parse_gaussian_record(record.at("alpha"), "letters[].alpha");
        if (gaussian_is_zero(*alpha)) {
            throw std::runtime_error(
                "dlog letter vanishes at the segment anchor (" + name + ")");
        }
        const fjson::Value& residue_records = record.at("residue");
        if (!residue_records.is_array() ||
            residue_records.items.size() != static_cast<size_t>(dimension)) {
            throw std::runtime_error("letters[].residue must be a d by d matrix");
        }
        std::unique_ptr<AcbMat> residue(new AcbMat(dimension, dimension));
        for (long row = 0; row < dimension; ++row) {
            const fjson::Value& row_records = residue_records.items[row];
            if (!row_records.is_array() ||
                row_records.items.size() != static_cast<size_t>(dimension)) {
                throw std::runtime_error("letters[].residue must be a d by d matrix");
            }
            for (long col = 0; col < dimension; ++col) {
                parse_scalar_pair(row_records.items[col],
                                  acb_mat_entry(residue->matrix, row, col),
                                  "letters[].residue[][]");
            }
        }
        std::unique_ptr<GaussQ> pole(new GaussQ());
        gaussian_divide_pole(alpha->re, alpha->im, beta->re, beta->im, *pole);
        long merge_index = -1;
        for (size_t previous = 0; previous < pole_exact.size(); ++previous) {
            if (pole_exact[previous]->equals(*pole)) {
                merge_index = static_cast<long>(previous);
                break;
            }
        }
        if (merge_index >= 0) {
            acb_mat_add(out.residues[merge_index]->matrix,
                        out.residues[merge_index]->matrix, residue->matrix, prec());
        } else {
            pole_exact.push_back(std::move(pole));
            std::unique_ptr<Acb> pole_acb(new Acb());
            gauss_to_acb(*pole_exact.back(), pole_acb->value);
            out.poles.push_back(std::move(pole_acb));
            out.residues.push_back(std::move(residue));
        }
    }
    out.sparse_residues = build_sparse_residues(out);
    // 留数内核选择：Sparse/Dense 强制，auto 按非零元密度判定。
    long nonzeros = 0;
    for (const auto& item : out.sparse_residues) nonzeros += item.nonzero_count;
    const long cells =
        static_cast<long>(out.residues.size()) * dimension * dimension;
    out.residue_density =
        cells == 0 ? 0.0 : static_cast<double>(nonzeros) / static_cast<double>(cells);
    switch (residue_mode) {
        case ResidueMode::Sparse:
            out.use_sparse = true;
            break;
        case ResidueMode::Dense:
            out.use_sparse = false;
            break;
        case ResidueMode::Auto:
        default:
            out.use_sparse =
                cells == 0 || out.residue_density <= kSparseDensityThreshold;
            break;
    }
}

// 从合并后的稠密留数建 CSR：精确零球（"0"/"0" 解析产物）全部跳过。
std::vector<SparseResidue> build_sparse_residues(const SegmentSystem& system) {
    std::vector<SparseResidue> sparse(system.residues.size());
    for (size_t j = 0; j < system.residues.size(); ++j) {
        SparseResidue& item = sparse[j];
        item.dimension = system.dimension;
        item.row_ptr.assign(system.dimension + 1, 0);
        std::vector<bool> column_used(system.dimension, false);
        for (long row = 0; row < system.dimension; ++row) {
            for (long col = 0; col < system.dimension; ++col) {
                if (!acb_is_zero(
                        acb_mat_entry(system.residues[j]->matrix, row, col))) {
                    item.cols.push_back(col);
                    item.vals.emplace_back(new Acb(
                        acb_mat_entry(system.residues[j]->matrix, row, col)));
                    column_used[col] = true;
                }
            }
            item.row_ptr[row + 1] = static_cast<long>(item.cols.size());
        }
        item.nonzero_count = static_cast<long>(item.cols.size());
        for (long col = 0; col < system.dimension; ++col) {
            if (column_used[col]) item.state_rows.push_back(col);
        }
    }
    return sparse;
}

// total += R_j · state，只遍历 R_j 非零元；等价于稠密和式中剔除 +0 项。
void sparse_pole_accumulate(
    acb_mat_t total, const SparseResidue& sparse, const acb_mat_t state) {
    const long columns = acb_mat_ncols(state);
    acb_t product;
    acb_init(product);
    for (long row = 0; row < sparse.dimension; ++row) {
        for (long idx = sparse.row_ptr[row]; idx < sparse.row_ptr[row + 1]; ++idx) {
            const long col = sparse.cols[idx];
            for (long k = 0; k < columns; ++k) {
                acb_mul(product, sparse.vals[idx]->value,
                        acb_mat_entry(state, col, k), prec());
                acb_add(acb_mat_entry(total, row, k),
                        acb_mat_entry(total, row, k), product, prec());
            }
        }
    }
    acb_clear(product);
}

// ---------- 数值内核 ----------

// 极点状态递推：解 Taylor 系数 solution[0..order]，cols 列初值并排放在 v0。
// 镜像 PartialFractionSystem.solution_taylor_coefficients 的运算序列。
// 留数累加按 system.use_sparse 分派：sparse 为 CSR 非零元跳乘（state 只在
// state_rows 行上推进），dense 为 acb_mat_mul+acb_mat_add（全行推进）。
std::vector<std::unique_ptr<AcbMat>> pole_recurrence(
    const SegmentSystem& system, const acb_t center, long order, const acb_mat_t v0,
    long ragged_stop = -1) {
    const long dimension = system.dimension;
    const long columns = acb_mat_ncols(v0);
    const size_t pole_count = system.poles.size();
    // ragged_stop：0 号列只需要到该阶（含）的系数块（前缀性质：主链前缀）。
    // 尾部阶（degree+1 > ragged_stop）不再推进/相乘 0 号列，改用 1 列临时量。
    const bool ragged =
        ragged_stop >= 0 && columns >= 2 && ragged_stop < order;
    std::vector<std::unique_ptr<Acb>> offsets(pole_count);
    for (size_t j = 0; j < pole_count; ++j) {
        offsets[j].reset(new Acb());
        acb_sub(offsets[j]->value, system.poles[j]->value, center, prec());
        if (acb_contains_zero(offsets[j]->value)) {
            throw std::runtime_error(
                "expansion center coincides with a pole; use local singular data");
        }
    }
    std::vector<std::unique_ptr<AcbMat>> states(pole_count);
    for (auto& state : states) state.reset(new AcbMat(dimension, columns));
    std::vector<std::unique_ptr<AcbMat>> solution(order + 1);
    solution[0].reset(new AcbMat(v0));
    AcbMat total(dimension, columns);
    AcbMat dense_workspace(dimension, columns);  // dense 内核乘积临时量
    // 尾部单列工作区：states_tail[j] 为各极点第 1 列状态，宽 1。
    std::vector<std::unique_ptr<AcbMat>> states_tail;
    AcbMat total_tail(dimension, 1);
    AcbMat dense_workspace_tail(dimension, 1);
    acb_t difference;
    acb_init(difference);
    for (long degree = 0; degree < order; ++degree) {
        const bool tail = ragged && degree >= ragged_stop;
        if (tail && states_tail.empty()) {
            states_tail.resize(pole_count);
            for (size_t j = 0; j < pole_count; ++j) {
                const SparseResidue& sparse = system.sparse_residues[j];
                if (sparse.empty()) continue;
                states_tail[j].reset(new AcbMat(dimension, 1));
                if (system.use_sparse) {
                    // 只拷贝该极点会读到的行；其余行从不被读取。
                    for (long row : sparse.state_rows) {
                        acb_set(acb_mat_entry(states_tail[j]->matrix, row, 0),
                                acb_mat_entry(states[j]->matrix, row, 1));
                    }
                } else {
                    for (long row = 0; row < dimension; ++row) {
                        acb_set(acb_mat_entry(states_tail[j]->matrix, row, 0),
                                acb_mat_entry(states[j]->matrix, row, 1));
                    }
                }
            }
        }
        if (tail) {
            acb_mat_zero(total_tail.matrix);
            for (size_t j = 0; j < pole_count; ++j) {
                const SparseResidue& sparse = system.sparse_residues[j];
                if (sparse.empty()) continue;  // 全零留数极点无贡献
                AcbMat& state = *states_tail[j];
                if (system.use_sparse) {
                    for (long row : sparse.state_rows) {
                        acb_sub(difference,
                                acb_mat_entry(state.matrix, row, 0),
                                acb_mat_entry(solution[degree]->matrix, row, 1), prec());
                        acb_div(acb_mat_entry(state.matrix, row, 0),
                                difference, offsets[j]->value, prec());
                    }
                    sparse_pole_accumulate(total_tail.matrix, sparse, state.matrix);
                } else {
                    for (long row = 0; row < dimension; ++row) {
                        acb_sub(difference,
                                acb_mat_entry(state.matrix, row, 0),
                                acb_mat_entry(solution[degree]->matrix, row, 1), prec());
                        acb_div(acb_mat_entry(state.matrix, row, 0),
                                difference, offsets[j]->value, prec());
                    }
                    acb_mat_mul(dense_workspace_tail.matrix,
                                system.residues[j]->matrix, state.matrix, prec());
                    acb_mat_add(total_tail.matrix, total_tail.matrix,
                                dense_workspace_tail.matrix, prec());
                }
            }
            solution[degree + 1].reset(new AcbMat(dimension, columns));
            acb_mat_zero(solution[degree + 1]->matrix);
            for (long row = 0; row < dimension; ++row) {
                acb_div_ui(acb_mat_entry(solution[degree + 1]->matrix, row, 1),
                           acb_mat_entry(total_tail.matrix, row, 0),
                           static_cast<ulong>(degree + 1), prec());
            }
            continue;
        }
        acb_mat_zero(total.matrix);
        for (size_t j = 0; j < pole_count; ++j) {
            const SparseResidue& sparse = system.sparse_residues[j];
            if (sparse.empty()) continue;  // 全零留数极点无贡献
            AcbMat& state = *states[j];
            if (system.use_sparse) {
                for (long row : sparse.state_rows) {
                    for (long col = 0; col < columns; ++col) {
                        acb_sub(difference,
                                acb_mat_entry(state.matrix, row, col),
                                acb_mat_entry(solution[degree]->matrix, row, col), prec());
                        acb_div(acb_mat_entry(state.matrix, row, col),
                                difference, offsets[j]->value, prec());
                    }
                }
                sparse_pole_accumulate(total.matrix, sparse, state.matrix);
            } else {
                for (long row = 0; row < dimension; ++row) {
                    for (long col = 0; col < columns; ++col) {
                        acb_sub(difference,
                                acb_mat_entry(state.matrix, row, col),
                                acb_mat_entry(solution[degree]->matrix, row, col), prec());
                        acb_div(acb_mat_entry(state.matrix, row, col),
                                difference, offsets[j]->value, prec());
                    }
                }
                acb_mat_mul(dense_workspace.matrix,
                            system.residues[j]->matrix, state.matrix, prec());
                acb_mat_add(total.matrix, total.matrix, dense_workspace.matrix, prec());
            }
        }
        solution[degree + 1].reset(new AcbMat(dimension, columns));
        for (long row = 0; row < dimension; ++row) {
            for (long col = 0; col < columns; ++col) {
                acb_div_ui(acb_mat_entry(solution[degree + 1]->matrix, row, col),
                           acb_mat_entry(total.matrix, row, col),
                           static_cast<ulong>(degree + 1), prec());
            }
        }
    }
    acb_clear(difference);
    return solution;
}

// Horner 求值（单列提取）：按列取级数系数，value = value*delta + c。
void horner_evaluate_column(
    const std::vector<std::unique_ptr<AcbMat>>& coefficients, long count, long column,
    const acb_t delta, acb_mat_t out) {
    const long rows = coefficients[0]->nrows();
    acb_t product;
    acb_init(product);
    for (long row = 0; row < rows; ++row) {
        acb_set(acb_mat_entry(out, row, 0),
                acb_mat_entry(coefficients[count - 1]->matrix, row, column));
    }
    for (long degree = count - 2; degree >= 0; --degree) {
        for (long row = 0; row < rows; ++row) {
            acb_mul(product, acb_mat_entry(out, row, 0), delta, prec());
            acb_add(acb_mat_entry(out, row, 0),
                    product,
                    acb_mat_entry(coefficients[degree]->matrix, row, column), prec());
        }
    }
    acb_clear(product);
}

// 段后中点重启：镜像 acb_midpoint_matrix。
void midpoint_restart(acb_mat_t matrix) {
    acb_t mid;
    acb_init(mid);
    for (long row = 0; row < acb_mat_nrows(matrix); ++row) {
        for (long col = 0; col < acb_mat_ncols(matrix); ++col) {
            acb_get_mid(mid, acb_mat_entry(matrix, row, col));
            acb_set(acb_mat_entry(matrix, row, col), mid);
        }
    }
    acb_clear(mid);
}

// 最近极点距离（镜像 nearest_singularity_distance：按中点 double 取最小）。
void nearest_pole_distance(const SegmentSystem& system, const acb_t point, arb_t out) {
    if (system.poles.empty()) {
        throw std::runtime_error("no poles declared for automatic radius control");
    }
    arb_t distance;
    arb_init(distance);
    bool first = true;
    for (const auto& pole : system.poles) {
        acb_t shifted;
        acb_init(shifted);
        acb_sub(shifted, point, pole->value, prec());
        acb_abs(distance, shifted, prec());
        acb_clear(shifted);
        if (first || arb_mid_double(distance) < arb_mid_double(out)) {
            arb_set(out, distance);
            first = false;
        }
    }
    arb_clear(distance);
}

// 列向量无穷范数（镜像 vector_norm_inf：按中点 double 取最大）。
void vector_norm_inf(const acb_mat_t vector, arb_t out) {
    if (acb_mat_ncols(vector) != 1) {
        throw std::runtime_error("vector_norm_inf requires a column vector");
    }
    arb_t magnitude;
    arb_init(magnitude);
    bool first = true;
    for (long row = 0; row < acb_mat_nrows(vector); ++row) {
        acb_abs(magnitude, acb_mat_entry(vector, row, 0), prec());
        if (first || arb_mid_double(magnitude) > arb_mid_double(out)) {
            arb_set(out, magnitude);
            first = false;
        }
    }
    arb_clear(magnitude);
}

void relative_difference_inf(const acb_mat_t primary, const acb_mat_t reference, arb_t out) {
    if (acb_mat_nrows(primary) != acb_mat_nrows(reference) ||
        acb_mat_ncols(primary) != acb_mat_ncols(reference)) {
        throw std::runtime_error("relative-difference dimensions do not match");
    }
    arb_t norm;
    arb_init(norm);
    vector_norm_inf(reference, norm);
    if (arb_contains_zero(norm)) {
        arb_clear(norm);
        throw std::runtime_error("reference matrix norm contains zero");
    }
    acb_mat_t difference;
    acb_mat_init(difference, acb_mat_nrows(primary), acb_mat_ncols(primary));
    acb_mat_sub(difference, reference, primary, prec());
    arb_t difference_norm;
    arb_init(difference_norm);
    vector_norm_inf(difference, difference_norm);
    arb_div(out, difference_norm, norm, prec());
    arb_clear(difference_norm);
    acb_mat_clear(difference);
    arb_clear(norm);
}

// a < b 是否严格成立（球语义：b-a 整个为正）。
bool arb_definitely_less(const arb_t a, const arb_t b) {
    arb_t gap;
    arb_init(gap);
    arb_sub(gap, b, a, prec());
    bool result = arb_is_positive(gap);
    arb_clear(gap);
    return result;
}

struct SegmentOutcome {
    std::vector<std::unique_ptr<AcbMat>> snapshots;  // 主链各节点值（d x 1）
    std::unique_ptr<AcbMat> final_primary;
    std::unique_ptr<AcbMat> final_reference;
    double primary_seconds = 0.0;
    double reference_seconds = 0.0;
    double kernel_seconds = 0.0;
    long node_count = 0;
    std::vector<double> step_ratios;
};

struct RadiusGuard {
    arb_t nearest;
    arb_t radius;
    arb_t delta_abs;
    acb_t delta;

    RadiusGuard() {
        arb_init(nearest);
        arb_init(radius);
        arb_init(delta_abs);
        acb_init(delta);
    }
    ~RadiusGuard() {
        arb_clear(nearest);
        arb_clear(radius);
        arb_clear(delta_abs);
        acb_clear(delta);
    }

    // 返回 |delta|/nearest 的中点比值；不满足半径条件时抛错。
    double check(
        const SegmentSystem& system, const acb_t start, const acb_t target,
        const arb_t radius_fraction, long segment_index) {
        nearest_pole_distance(system, start, nearest);
        arb_mul(radius, nearest, radius_fraction, prec());
        acb_sub(delta, target, start, prec());
        acb_abs(delta_abs, delta, prec());
        if (!(arb_mid_double(delta_abs) < arb_mid_double(radius))) {
            throw std::runtime_error(
                "segment leaves its Cauchy disk; refine the path (segment " +
                std::to_string(segment_index + 1) + ")");
        }
        return arb_mid_double(delta_abs) / arb_mid_double(nearest);
    }
};

// 前缀复用输运（唯一保留路线）：两链同精度只差阶数 → 同一步内 reference 级数的
// 0..primary_order 阶系数就是 primary 级数（前缀性质）。首步两链同起点，单列递推到
// 参考阶即可，主链取系数前缀；此后各步两链各自中点重启、起点不同，用 ragged 双列
// 推进（primary 列只推到 primary_order，尾部乘法只算 1 列）。
// 逐列算术与独立双链完全相同，结果与逐链完整重算逐位一致。
SegmentOutcome run_transport(
    const SegmentJob& job, const acb_mat_t primary_vector,
    const acb_mat_t reference_vector, long primary_order, long reference_order,
    const arb_t radius_fraction) {
    SegmentOutcome outcome;
    const long dimension = job.system.dimension;
    AcbMat primary(dimension, 1);
    acb_mat_set(primary.matrix, primary_vector);
    AcbMat reference(dimension, 1);
    acb_mat_set(reference.matrix, reference_vector);
    AcbMat pair(dimension, 2);
    const long node_count = static_cast<long>(job.points.size()) + 1;
    outcome.snapshots.resize(node_count);
    outcome.snapshots[0].reset(new AcbMat(primary.matrix));
    outcome.node_count = node_count;
    RadiusGuard guard;
    Clock::time_point kernel_start = Clock::now();
    bool chains_coincide = true;
    for (long index = 0; index < node_count - 1; ++index) {
        acb_srcptr start =
            index == 0 ? job.start->value : job.points[index - 1]->value;
        acb_srcptr target = job.points[index]->value;
        outcome.step_ratios.push_back(
            guard.check(job.system, start, target, radius_fraction, index));
        if (chains_coincide) {
            auto coefficients =
                pole_recurrence(job.system, start, reference_order, primary.matrix);
            horner_evaluate_column(
                coefficients, primary_order + 1, 0, guard.delta, primary.matrix);
            horner_evaluate_column(
                coefficients, reference_order + 1, 0, guard.delta, reference.matrix);
        } else {
            for (long row = 0; row < dimension; ++row) {
                acb_set(acb_mat_entry(pair.matrix, row, 0),
                        acb_mat_entry(primary.matrix, row, 0));
                acb_set(acb_mat_entry(pair.matrix, row, 1),
                        acb_mat_entry(reference.matrix, row, 0));
            }
            auto coefficients = pole_recurrence(
                job.system, start, reference_order, pair.matrix, primary_order);
            horner_evaluate_column(
                coefficients, primary_order + 1, 0, guard.delta, primary.matrix);
            horner_evaluate_column(
                coefficients, reference_order + 1, 1, guard.delta, reference.matrix);
        }
        midpoint_restart(primary.matrix);
        midpoint_restart(reference.matrix);
        outcome.snapshots[index + 1].reset(new AcbMat(primary.matrix));
        // 前缀截断与全阶求值的中点重启值一般不同，此后两链起点分叉。
        chains_coincide = false;
    }
    outcome.kernel_seconds = elapsed_seconds(kernel_start, Clock::now());
    outcome.primary_seconds = outcome.kernel_seconds;
    outcome.reference_seconds = outcome.kernel_seconds;
    outcome.final_primary.reset(new AcbMat(primary.matrix));
    outcome.final_reference.reset(new AcbMat(reference.matrix));
    return outcome;
}

const std::set<std::string>& request_keys() {
    static const std::set<std::string> keys = {
        "schema", "backendPackagePath", "masterDigest", "dimension", "segments",
        "pathPlanning", "singularityMode", "boundary", "workingPrecisionDigits",
        "primaryOrder", "referenceOrder", "targetRelativeError", "certificationMode",
        "messageLanguage", "columnVectorConvention", "dlogStatus",
    };
    return keys;
}

}  // namespace

ResidueMode parse_residue_mode(const std::string& text) {
    if (text == "auto") return ResidueMode::Auto;
    if (text == "sparse") return ResidueMode::Sparse;
    if (text == "dense") return ResidueMode::Dense;
    throw std::runtime_error(
        "unknown residue kernel '" + text + "'; expected auto, sparse or dense");
}

std::string residue_mode_text(ResidueMode mode) {
    switch (mode) {
        case ResidueMode::Sparse: return "sparse";
        case ResidueMode::Dense: return "dense";
        case ResidueMode::Auto:
        default: return "auto";
    }
}

fjson::Value run_request(const fjson::Value& request, const EngineOptions& options) {
    require_exact_keys(request, request_keys(), "request");
    if (!request.at("schema").is_string() ||
        request.at("schema").text != "madstree_flintnde_evaluate_v1") {
        throw std::runtime_error("unsupported MadStree-FlintNDE schema");
    }
    if (!request.at("pathPlanning").is_bool()) {
        throw std::runtime_error("pathPlanning must be a JSON boolean");
    }
    if (request.at("pathPlanning").boolean) {
        return refusal(
            "pathPlanningNotImplemented",
            "the FlintNDE engine currently only implements the direct user-point "
            "route; enable path planning on the Python side or supply explicit nodes",
            &request);
    }
    const std::string& singularity_mode = request.at("singularityMode").text;
    if (singularity_mode != "avoid") {
        return refusal(
            "singularityJumpNotImplemented",
            "direct user-point transport supports avoid mode only in the FlintNDE "
            "engine",
            &request);
    }
    const std::string& certification_mode = request.at("certificationMode").text;
    if (certification_mode != "certified") {
        return refusal(
            "embeddedCertificationNotImplemented",
            "the FlintNDE engine currently implements the certified double chain "
            "only",
            &request);
    }
    const long dimension = request.at("dimension").as_integer("dimension");
    if (dimension < 1) throw std::runtime_error("dimension must be positive");
    const long digits = request.at("workingPrecisionDigits").as_integer("digits");
    if (digits < 20) throw std::runtime_error("workingPrecisionDigits must be >= 20");
    const long primary_order = request.at("primaryOrder").as_integer("primaryOrder");
    const long reference_order = request.at("referenceOrder").as_integer("referenceOrder");
    if (reference_order <= primary_order) {
        throw std::runtime_error("referenceOrder must exceed primaryOrder");
    }
    set_thread_prec(working_bits(digits));

    const fjson::Value& boundary = request.at("boundary");
    if (!boundary.is_object() || !boundary.find("kind") ||
        boundary.at("kind").text != "finite") {
        return refusal(
            "regularSingularBoundaryNotImplemented",
            "the FlintNDE engine currently accepts finite boundary vectors only; "
            "generate the anchor boundary with FlintNDE or supply finite values",
            &request);
    }
    require_exact_keys(boundary, {"kind", "values"}, "boundary");
    AcbMat primary_initial(dimension, 1);
    const fjson::Value& values = boundary.at("values");
    if (!values.is_array() || values.items.size() != static_cast<size_t>(dimension)) {
        throw std::runtime_error("finite boundary vector dimension mismatch");
    }
    for (long row = 0; row < dimension; ++row) {
        parse_scalar_pair(
            values.items[row], acb_mat_entry(primary_initial.matrix, row, 0),
            "boundary.values[]");
    }

    const fjson::Value& segments = request.at("segments");
    if (!segments.is_array() || segments.items.empty()) {
        throw std::runtime_error("segments must be a nonempty list");
    }

    arb_t radius_fraction;
    arb_init(radius_fraction);
    arb_set_str(radius_fraction, "0.6", prec());

    std::vector<SegmentJob> jobs;
    jobs.reserve(segments.items.size());
    for (size_t index = 0; index < segments.items.size(); ++index) {
        const fjson::Value& record = segments.items[index];
        require_exact_keys(
            record, {"start", "points", "letters", "fromUserIndex", "userIndices"},
            "segments[]");
        SegmentJob job;
        job.start.reset(new Acb());
        if (!record.at("start").is_string()) {
            throw std::runtime_error("segments[].start must be a scalar string");
        }
        parse_scalar_text(record.at("start").text, job.start->value);
        const fjson::Value& points = record.at("points");
        if (!points.is_array() || points.items.empty()) {
            throw std::runtime_error("segments[].points must be a nonempty list");
        }
        for (const auto& point : points.items) {
            std::unique_ptr<Acb> target(new Acb());
            parse_scalar_pair(point, target->value, "segments[].points[]");
            job.points.push_back(std::move(target));
        }
        const fjson::Value& user_indices = record.at("userIndices");
        if (!user_indices.is_array() ||
            user_indices.items.size() != job.points.size()) {
            throw std::runtime_error("segments[].userIndices must match points");
        }
        for (const auto& item : user_indices.items) {
            job.user_indices.push_back(item.as_integer("userIndices[]"));
        }
        build_segment_system(
            record.at("letters"), dimension, options.residue_mode,
            request.at("masterDigest").text + "-" + std::to_string(index),
            job.system);
        // direct_user_point_path 语义：起点与用户点不得落在极点上。
        for (const auto& pole : job.system.poles) {
            acb_t shifted;
            acb_init(shifted);
            acb_sub(shifted, job.start->value, pole->value, prec());
            bool start_on_pole = acb_contains_zero(shifted);
            acb_clear(shifted);
            if (start_on_pole) {
                arb_clear(radius_fraction);
                throw std::runtime_error(
                    "the path start coincides with a pole and cannot be transported "
                    "as an ordinary start");
            }
            for (const auto& point : job.points) {
                acb_t offset;
                acb_init(offset);
                acb_sub(offset, point->value, pole->value, prec());
                bool point_on_pole = acb_contains_zero(offset);
                acb_clear(offset);
                if (point_on_pole) {
                    arb_clear(radius_fraction);
                    throw std::runtime_error(
                        "a user point coincides with a pole and cannot be used as a "
                        "transport node");
                }
            }
        }
        jobs.push_back(std::move(job));
    }

    const std::string& language = request.at("messageLanguage").text;
    std::unique_ptr<AcbMat> primary_vector(new AcbMat(dimension, 1));
    acb_mat_set(primary_vector->matrix, primary_initial.matrix);
    std::unique_ptr<AcbMat> reference_vector(new AcbMat(dimension, 1));
    acb_mat_set(reference_vector->matrix, primary_initial.matrix);

    fjson::Value segments_out = jarr();
    double primary_seconds_total = 0.0;
    double reference_seconds_total = 0.0;
    double kernel_seconds_total = 0.0;
    bool any_sparse = false;
    bool any_dense = false;
    bool all_met = true;
    arb_t difference;
    arb_t target;
    arb_init(difference);
    arb_init(target);
    if (arb_set_str(target, request.at("targetRelativeError").text.c_str(), prec()) != 0) {
        arb_clear(difference);
        arb_clear(target);
        arb_clear(radius_fraction);
        throw std::runtime_error("cannot parse targetRelativeError");
    }

    for (size_t index = 0; index < jobs.size(); ++index) {
        const SegmentJob& job = jobs[index];
        Clock::time_point segment_wall_start = Clock::now();
        SegmentOutcome outcome = run_transport(
            job, primary_vector->matrix, reference_vector->matrix,
            primary_order, reference_order, radius_fraction);
        if (job.system.use_sparse) {
            any_sparse = true;
        } else {
            any_dense = true;
        }
        double segment_wall = elapsed_seconds(segment_wall_start, Clock::now());
        relative_difference_inf(
            outcome.final_primary->matrix, outcome.final_reference->matrix, difference);
        bool met = arb_definitely_less(difference, target);
        all_met = all_met && met;

        fjson::Value segment_record = jobj();
        segment_record.fields["segmentIndex"] = jnum(static_cast<double>(index));
        segment_record.fields["pathPlanning"] = jbool(false);
        fjson::Value nodes = jarr();
        fjson::Value node_point = jobj();
        node_point.fields["real"] = jstr(
            arb_midpoint_string(acb_realref(job.start->value), 40));
        node_point.fields["imag"] = jstr(
            arb_midpoint_string(acb_imagref(job.start->value), 40));
        nodes.items.push_back(node_point);
        for (const auto& point : job.points) {
            fjson::Value item = jobj();
            item.fields["real"] = jstr(arb_midpoint_string(acb_realref(point->value), 40));
            item.fields["imag"] = jstr(arb_midpoint_string(acb_imagref(point->value), 40));
            nodes.items.push_back(item);
        }
        segment_record.fields["actualNodes"] = nodes;
        segment_record.fields["nodeCount"] = jnum(static_cast<double>(outcome.node_count));
        segment_record.fields["coveredSampleCount"] = jnum(0.0);
        fjson::Value plan_report = jobj();
        plan_report.fields["schema"] = jstr("flintnde_direct_user_nodes_v1");
        plan_report.fields["planning_action"] = jstr("disabled_use_user_points_as_nodes");
        plan_report.fields["singularity_mode"] = jstr("avoid");
        plan_report.fields["engine"] = jstr("flintnde_engine");
        plan_report.fields["residueKernel"] =
            jstr(job.system.use_sparse ? "sparse" : "dense");
        plan_report.fields["residueDensity"] = jnum(job.system.residue_density);
        plan_report.fields["segment_count"] =
            jnum(static_cast<double>(outcome.node_count - 1));
        plan_report.fields["singularity_jump_count"] = jnum(0.0);
        plan_report.fields["covered_sample_count"] = jnum(0.0);
        plan_report.fields["working_precision_bits"] =
            jnum(static_cast<double>(prec()));
        fjson::Value step_ratios = jarr();
        for (double ratio : outcome.step_ratios) step_ratios.items.push_back(jnum(ratio));
        plan_report.fields["step_over_convergence_radius"] = step_ratios;
        segment_record.fields["planReport"] = plan_report;
        fjson::Value assignments = jarr();
        fjson::Value point_values = jarr();
        for (size_t point_index = 0; point_index < job.points.size(); ++point_index) {
            fjson::Value assignment = jobj();
            assignment.fields["userIndex"] =
                jnum(static_cast<double>(job.user_indices[point_index]));
            assignment.fields["localPointIndex"] = jnum(static_cast<double>(point_index));
            assignment.fields["source"] = jstr("node_snapshot");
            assignment.fields["nodeIndex"] = jnum(static_cast<double>(point_index + 1));
            assignment.fields["segmentIndex"] = jnum(-1.0);
            assignments.items.push_back(assignment);
            fjson::Value point_record = jobj();
            point_record.fields["userIndex"] =
                jnum(static_cast<double>(job.user_indices[point_index]));
            point_record.fields["localPointIndex"] = jnum(static_cast<double>(point_index));
            point_record.fields["source"] = jstr("node_snapshot");
            point_record.fields["evaluationAlgorithm"] = fjson::Value();
            point_record.fields["values"] =
                vector_records(outcome.snapshots[point_index + 1]->matrix, digits);
            point_values.items.push_back(point_record);
        }
        segment_record.fields["pointAssignments"] = assignments;
        segment_record.fields["pointValues"] = point_values;
        segment_record.fields["relativeDifferenceInf"] =
            jstr(arb_ball_string(difference, digits));
        segment_record.fields["relativeDifferenceMidpoint"] =
            jnum(arb_mid_double(difference));
        segment_record.fields["targetRelativeErrorMet"] = jbool(met);
        segment_record.fields["certificationMode"] = jstr("certified");
        segment_record.fields["continuationReady"] = jbool(true);
        segment_record.fields["endpointValues"] =
            vector_records(outcome.final_primary->matrix, digits);
        segment_record.fields["wallSeconds"] = jnum(segment_wall);
        segments_out.items.push_back(segment_record);
        primary_seconds_total += outcome.primary_seconds;
        reference_seconds_total += outcome.reference_seconds;
        kernel_seconds_total += outcome.kernel_seconds;
        primary_vector = std::move(outcome.final_primary);
        reference_vector = std::move(outcome.final_reference);
    }

    std::string message =
        language == "CN"
            ? "FlintNDE 引擎已把用户给定点严格作为节点，在本进程内完成认证双链输运。"
            : "The FlintNDE engine evaluated every user point as a node with the "
              "certified double chain in one process.";

    fjson::Value result = jobj();
    result.fields["status"] = jstr("success");
    result.fields["schema"] = request.at("schema");
    result.fields["workerPid"] = jnum(static_cast<double>(_getpid()));
    result.fields["executionAction"] = jstr("execute_user_nodes");
    result.fields["masterDigest"] = request.at("masterDigest");
    result.fields["dimension"] = jnum(static_cast<double>(dimension));
    result.fields["backendPackagePath"] = request.at("backendPackagePath");
    result.fields["messageLanguage"] = request.at("messageLanguage");
    result.fields["message"] = jstr(message);
    result.fields["pathPlanning"] = jbool(false);
    result.fields["segmentCount"] = jnum(static_cast<double>(jobs.size()));
    result.fields["segments"] = segments_out;
    result.fields["singularityClassifications"] = jarr();
    fjson::Value boundary_report = jobj();
    boundary_report.fields["kind"] = jstr("finite");
    boundary_report.fields["path"] = jarr();
    boundary_report.fields["seconds"] = jnum(0.0);
    result.fields["boundary"] = boundary_report;
    result.fields["planningSeconds"] = jnum(0.0);
    result.fields["primarySeconds"] = jnum(primary_seconds_total);
    result.fields["referenceSeconds"] = jnum(reference_seconds_total);
    result.fields["targetRelativeErrorMet"] = jbool(all_met);
    result.fields["finalValues"] = vector_records(primary_vector->matrix, digits);
    fjson::Value metadata = jobj();
    metadata.fields["engine"] = jstr("flintnde_engine");
    metadata.fields["flintVersion"] = jstr(flint_version);
    metadata.fields["flintThreads"] = jnum(static_cast<double>(options.flint_threads));
    metadata.fields["residuesRequested"] = jstr(residue_mode_text(options.residue_mode));
    metadata.fields["residueKernel"] = jstr(
        any_sparse && any_dense ? "mixed" : (any_sparse ? "sparse" : "dense"));
    metadata.fields["kernelSeconds"] = jnum(kernel_seconds_total);
    metadata.fields["workingPrecisionBits"] =
        jnum(static_cast<double>(prec()));
    result.fields["flintndeEngine"] = metadata;

    arb_clear(difference);
    arb_clear(target);
    arb_clear(radius_fraction);
    return result;
}

fjson::Value run_batch(const std::string& raw, const EngineOptions& options) {
    // 批量请求体积可达数百 MB：顶层除 tasks 外的小字段建 DOM，tasks 元素只记
    // 文本 span，worker 各自解析、用完即弃，主进程不保留整批任务 DOM。
    fjson::TopLevelScan scan = fjson::scan_object_streaming(raw, "tasks");
    if (!scan.has_streamed_array) {
        throw std::runtime_error("ep batch request must contain the tasks array");
    }
    require_exact_keys(
        scan.fields, {"schema", "parallelTaskCount", "messageLanguage"},
        "ep batch request");
    if (scan.fields.at("schema").text != "madstree_flintnde_ep_batch_v1") {
        throw std::runtime_error("unsupported MadStree-FlintNDE ep batch schema");
    }
    const long requested =
        scan.fields.at("parallelTaskCount").as_integer("parallelTaskCount");
    const std::string& language = scan.fields.at("messageLanguage").text;
    const size_t task_count = scan.spans.size();
    if (task_count == 0) {
        throw std::runtime_error("ep batch tasks must be a nonempty list");
    }
    const long effective = std::min(
        options.worker_override > 0 ? options.worker_override : requested,
        static_cast<long>(task_count));
    const long prec_bits = prec();

    std::vector<fjson::Value> results(task_count);
    std::vector<fjson::Value> task_eps(task_count);
    std::atomic<size_t> next_index(0);
    Clock::time_point wall_start = Clock::now();
    auto worker = [&]() {
        set_thread_prec(prec_bits);
        while (true) {
            size_t index = next_index.fetch_add(1);
            if (index >= task_count) break;
            try {
                fjson::Value task = fjson::parse_span(
                    raw.c_str(), scan.spans[index].first, scan.spans[index].second);
                require_exact_keys(task, {"ep", "request"}, "tasks[]");
                task_eps[index] = task.at("ep");
                results[index] = run_request(task.at("request"), options);
            } catch (const std::exception& error) {
                // 单任务失败降级为该任务的 failed 结果，不拖垮整批线程池。
                results[index] = refusal("failed", error.what(), nullptr);
            }
        }
    };
    std::vector<std::thread> pool;
    const long worker_count = std::max<long>(1, effective);
    for (long index = 0; index < worker_count; ++index) pool.emplace_back(worker);
    for (auto& thread : pool) thread.join();
    double wall_seconds = elapsed_seconds(wall_start, Clock::now());

    fjson::Value results_array = jarr();
    size_t failed_tasks = 0;
    for (size_t index = 0; index < task_count; ++index) {
        fjson::Value item = jobj();
        item.fields["ep"] = task_eps[index];
        item.fields["result"] = results[index];
        const fjson::Value* status = results[index].find("status");
        if (status == nullptr || !status->is_string() || status->text != "success") {
            ++failed_tasks;
        }
        results_array.items.push_back(item);
    }
    std::string message =
        language == "CN"
            ? "FlintNDE 引擎 ep 任务池：请求并行数 " + std::to_string(requested) +
                  "，实际并行数 " + std::to_string(effective) + "。"
            : "FlintNDE engine ep task pool used " + std::to_string(effective) +
                  " of " + std::to_string(requested) + " requested workers.";
    fjson::Value result = jobj();
    result.fields["status"] = jstr(failed_tasks == 0 ? "success" : "failed");
    result.fields["schema"] = scan.fields.at("schema");
    result.fields["message"] = jstr(message);
    result.fields["messageLanguage"] = scan.fields.at("messageLanguage");
    result.fields["parallelTaskCountRequested"] = jnum(static_cast<double>(requested));
    result.fields["parallelTaskCountEffective"] = jnum(static_cast<double>(effective));
    result.fields["taskCount"] = jnum(static_cast<double>(task_count));
    if (failed_tasks != 0) {
        result.fields["failedTaskCount"] = jnum(static_cast<double>(failed_tasks));
    }
    result.fields["results"] = results_array;
    fjson::Value metadata = jobj();
    metadata.fields["engine"] = jstr("flintnde_engine");
    metadata.fields["flintVersion"] = jstr(flint_version);
    metadata.fields["flintThreads"] = jnum(static_cast<double>(options.flint_threads));
    metadata.fields["residuesRequested"] = jstr(residue_mode_text(options.residue_mode));
    metadata.fields["wallSeconds"] = jnum(wall_seconds);
    result.fields["flintndeEngine"] = metadata;
    return result;
}

}  // namespace fnde
