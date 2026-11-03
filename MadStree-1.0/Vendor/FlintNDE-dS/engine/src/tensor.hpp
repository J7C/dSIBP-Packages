// FlintNDE 引擎张量重建：把 MadStree 的 tensorSectors＋结构化 letter 贡献在引擎内
// 重建为精确 Q(i) 稠密留数，免 Python 侧展平与整批稠密 residue 传输。
//
// 移植自 flintnde.tensor_dlog 并保留其全部约定：局部基矩阵取 msPaperT/
// msHadamard/msIdentity2 去掉公共 1/sqrt(2) 后的整数高斯矩阵；state 编号
// IntegerDigits[m-1, 2, nSlots]，slot 1 为最高位，U = P (kron locals) S；
// U/U^{-1} 只按 sector 自身维数物化，全局 2^n 矩阵从不出现。
// 算术全部走 fmpq 精确车道，与 Python 车道（Fraction 解析＋精确求和后单次
// 取球）逐位一致；校验镜像适配器 flintnde_transport 的字段检查。
#ifndef FLINTNDE_ENGINE_TENSOR_HPP
#define FLINTNDE_ENGINE_TENSOR_HPP

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "json.hpp"
#include "util.hpp"

namespace fnde {

enum class TensorKind { Diagonal, Conjugated };
enum class TensorScope { Sector, Global };

struct TensorSectorData {
    std::string key;
    long dimension = 0;
    long offset = 0;      // 1-based 全局起始列
    long raw_count = 0;   // 2^slots
    std::vector<std::string> slot_bases;
    std::vector<long> raw_to_master_column;        // rawCount 个，0 表示商掉
    std::vector<long> master_representative_raw;   // 每个本地列的代表 raw（1-based）

    long first_global_index() const { return offset - 1; }
};

struct TensorContributionData {
    TensorKind kind = TensorKind::Diagonal;
    TensorScope scope = TensorScope::Sector;
    std::string sector_key;
    std::vector<GaussVal> vector;             // diagonal：sector 维或全局维
    std::vector<GaussVal> pre_vector;         // conjugated sector scope
    std::vector<long> letter_indicator;       // conjugated：本地维 0/1
    std::vector<long> middle_rows;            // conjugated global：0-based 本地行
    std::vector<long> middle_columns;         // conjugated global：0-based 全局列
    std::vector<GaussVal> middle_values;
    GaussVal post_factor;
};

// 局部基标签 -> 2x2 精确矩阵（正向）与逐元素逆。公共 1/sqrt(2) 因子在
// U^{-1} D U 中整体抵消，全程只用整数高斯矩阵（与 tensor_dlog 相同）。
inline const std::vector<std::vector<GaussVal>>& local_base(
    const std::string& tag, bool inverse) {
    static const std::map<std::string, std::vector<std::vector<GaussVal>>> bases = {
        {"massiveT", {{GaussVal::rational(1), GaussVal::gaussian(0, 1, -1, 1)},
                      {GaussVal::gaussian(0, 1, -1, 1), GaussVal::rational(1)}}},
        {"masslessH", {{GaussVal::rational(1), GaussVal::rational(1)},
                       {GaussVal::rational(1), GaussVal::rational(-1)}}},
        {"identity", {{GaussVal::rational(1), GaussVal::rational(0)},
                      {GaussVal::rational(0), GaussVal::rational(1)}}},
    };
    static const std::map<std::string, std::vector<std::vector<GaussVal>>> inverses = {
        {"massiveT", {{GaussVal::gaussian(1, 2, 0, 2), GaussVal::gaussian(0, 2, 1, 2)},
                      {GaussVal::gaussian(0, 2, 1, 2), GaussVal::gaussian(1, 2, 0, 2)}}},
        {"masslessH", {{GaussVal::gaussian(1, 2, 0, 2), GaussVal::gaussian(1, 2, 0, 2)},
                       {GaussVal::gaussian(1, 2, 0, 2), GaussVal::gaussian(-1, 2, 0, 2)}}},
        {"identity", {{GaussVal::rational(1), GaussVal::rational(0)},
                      {GaussVal::rational(0), GaussVal::rational(1)}}},
    };
    const auto& table = inverse ? inverses : bases;
    const auto item = table.find(tag);
    if (item == table.end()) {
        throw std::runtime_error("unknown local base tag: " + tag);
    }
    return item->second;
}

// 槽轴链式作用：axis 0 对最高位，复杂度 O(n 2^n)（镜像 apply_local_chain）。
inline void apply_local_chain(
    const std::vector<std::string>& slot_bases, bool inverse,
    std::vector<GaussVal>& vector) {
    const long rank = static_cast<long>(slot_bases.size());
    if (rank == 0) return;
    if (static_cast<long>(vector.size()) != (1 << rank)) {
        throw std::runtime_error("tensor chain dimension mismatch");
    }
    for (long axis = 0; axis < rank; ++axis) {
        const long step = 1 << (rank - 1 - axis);
        const long block = step << 1;
        const auto& matrix = local_base(slot_bases[static_cast<size_t>(axis)], inverse);
        const GaussVal& m00 = matrix[0][0];
        const GaussVal& m01 = matrix[0][1];
        const GaussVal& m10 = matrix[1][0];
        const GaussVal& m11 = matrix[1][1];
        for (long base = 0; base < static_cast<long>(vector.size()); base += block) {
            for (long offset = 0; offset < step; ++offset) {
                const long first = base + offset;
                const long second = first + step;
                GaussVal lower(vector[static_cast<size_t>(first)]);
                GaussVal upper(vector[static_cast<size_t>(second)]);
                GaussVal term_a;
                GaussVal term_b;
                gauss_mul(term_a, m00, lower);
                gauss_mul(term_b, m01, upper);
                gauss_add(vector[static_cast<size_t>(first)], term_a, term_b);
                gauss_mul(term_a, m10, lower);
                gauss_mul(term_b, m11, upper);
                gauss_add(vector[static_cast<size_t>(second)], term_a, term_b);
            }
        }
    }
}

class TensorReconstructor {
  public:
    TensorReconstructor(std::vector<TensorSectorData> sectors, long dimension)
        : sectors_(std::move(sectors)), dimension_(dimension) {
        for (size_t index = 0; index < sectors_.size(); ++index) {
            if (by_key_.count(sectors_[index].key)) {
                throw std::runtime_error(
                    "duplicate tensor sector " + sectors_[index].key);
            }
            by_key_[sectors_[index].key] = index;
        }
        // 平铺校验：按 offset 排序后必须无重叠覆盖 1..dimension。
        std::vector<size_t> ordered(sectors_.size());
        for (size_t index = 0; index < ordered.size(); ++index) ordered[index] = index;
        std::sort(ordered.begin(), ordered.end(), [this](size_t a, size_t b) {
            return sectors_[a].offset < sectors_[b].offset;
        });
        long expected = 1;
        for (size_t index : ordered) {
            if (sectors_[index].offset != expected) {
                throw std::runtime_error(
                    "tensor sectors must tile the global master order");
            }
            expected += sectors_[index].dimension;
        }
        if (expected - 1 != dimension_) {
            throw std::runtime_error(
                "tensor sector dimensions do not match the master count");
        }
    }

    long dimension() const { return dimension_; }
    const std::vector<TensorSectorData>& sectors() const { return sectors_; }

    size_t sector_index(const std::string& key, const std::string& label) const {
        const auto item = by_key_.find(key);
        if (item == by_key_.end()) {
            throw std::runtime_error(label + " refers to unknown tensor sector " + key);
        }
        return item->second;
    }

    // U 或 U^{-1} 作用于 sector 内 master 基向量（镜像 _apply）。
    std::vector<GaussVal> apply_sector(
        size_t sector_index_value, const std::vector<GaussVal>& vector,
        bool inverse) const {
        const TensorSectorData& sector = sectors_[sector_index_value];
        if (static_cast<long>(vector.size()) != sector.dimension) {
            throw std::runtime_error(sector.key + ": apply dimension mismatch");
        }
        std::vector<GaussVal> raw;
        raw.reserve(static_cast<size_t>(sector.raw_count));
        for (long column : sector.raw_to_master_column) {
            if (column == 0) {
                raw.emplace_back();
            } else {
                raw.push_back(vector[static_cast<size_t>(column - 1)]);
            }
        }
        apply_local_chain(sector.slot_bases, inverse, raw);
        std::vector<GaussVal> result;
        result.reserve(sector.master_representative_raw.size());
        for (long raw_index : sector.master_representative_raw) {
            result.push_back(raw[static_cast<size_t>(raw_index - 1)]);
        }
        return result;
    }

    // U 的稠密表示（仅 sector 自身维数），按列物化后转置缓存。
    const std::vector<std::vector<GaussVal>>& sector_u(const std::string& key) const {
        const auto cached = u_.find(key);
        if (cached != u_.end()) return cached->second;
        const size_t index = sector_index(key, "sector_u");
        const TensorSectorData& sector = sectors_[index];
        std::vector<std::vector<GaussVal>> columns;
        for (long column = 0; column < sector.dimension; ++column) {
            std::vector<GaussVal> unit(static_cast<size_t>(sector.dimension));
            for (long row = 0; row < sector.dimension; ++row) {
                unit[static_cast<size_t>(row)] =
                    GaussVal::rational(row == column ? 1 : 0);
            }
            columns.push_back(apply_sector(index, unit, false));
        }
        std::vector<std::vector<GaussVal>> matrix(
            static_cast<size_t>(sector.dimension),
            std::vector<GaussVal>(static_cast<size_t>(sector.dimension)));
        for (long column = 0; column < sector.dimension; ++column) {
            for (long row = 0; row < sector.dimension; ++row) {
                matrix[static_cast<size_t>(row)][static_cast<size_t>(column)] =
                    columns[static_cast<size_t>(column)][static_cast<size_t>(row)];
            }
        }
        return u_.emplace(key, std::move(matrix)).first->second;
    }

    const std::vector<std::vector<GaussVal>>& sector_u_inverse(const std::string& key) const {
        const auto cached = u_inverse_.find(key);
        if (cached != u_inverse_.end()) return cached->second;
        const size_t index = sector_index(key, "sector_u_inverse");
        const TensorSectorData& sector = sectors_[index];
        std::vector<std::vector<GaussVal>> columns;
        for (long column = 0; column < sector.dimension; ++column) {
            std::vector<GaussVal> unit(static_cast<size_t>(sector.dimension));
            for (long row = 0; row < sector.dimension; ++row) {
                unit[static_cast<size_t>(row)] =
                    GaussVal::rational(row == column ? 1 : 0);
            }
            columns.push_back(apply_sector(index, unit, true));
        }
        std::vector<std::vector<GaussVal>> matrix(
            static_cast<size_t>(sector.dimension),
            std::vector<GaussVal>(static_cast<size_t>(sector.dimension)));
        for (long column = 0; column < sector.dimension; ++column) {
            for (long row = 0; row < sector.dimension; ++row) {
                matrix[static_cast<size_t>(row)][static_cast<size_t>(column)] =
                    columns[static_cast<size_t>(column)][static_cast<size_t>(row)];
            }
        }
        return u_inverse_.emplace(key, std::move(matrix)).first->second;
    }

    // K = U^{-1} diag(indicator) U；只有 indicator 为 1 的行参与（镜像实现）。
    const std::vector<std::vector<GaussVal>>& conjugated_block(
        const std::string& key, const std::vector<long>& indicator) const {
        const auto cached = blocks_.find(std::make_pair(key, indicator));
        if (cached != blocks_.end()) return cached->second;
        const auto& u = sector_u(key);
        const auto& inverse = sector_u_inverse(key);
        const long size = static_cast<long>(u.size());
        std::vector<long> active;
        for (long index = 0; index < size; ++index) {
            if (indicator[static_cast<size_t>(index)]) active.push_back(index);
        }
        std::vector<std::vector<GaussVal>> result(
            static_cast<size_t>(size), std::vector<GaussVal>(static_cast<size_t>(size)));
        for (long row = 0; row < size; ++row) {
            const auto& inverse_row = inverse[static_cast<size_t>(row)];
            std::vector<GaussVal> accumulator(
                static_cast<size_t>(size));
            for (long source : active) {
                const GaussVal& coefficient = inverse_row[static_cast<size_t>(source)];
                if (coefficient.is_zero()) continue;
                const auto& basis_row = u[static_cast<size_t>(source)];
                for (long column = 0; column < size; ++column) {
                    const GaussVal& value = basis_row[static_cast<size_t>(column)];
                    if (value.is_zero()) continue;
                    GaussVal product;
                    gauss_mul(product, coefficient, value);
                    GaussVal sum;
                    gauss_add(sum, accumulator[static_cast<size_t>(column)], product);
                    accumulator[static_cast<size_t>(column)] = sum;
                }
            }
            result[static_cast<size_t>(row)] = std::move(accumulator);
        }
        return blocks_.emplace(std::make_pair(key, indicator), std::move(result))
            .first->second;
    }

    // 把一组结构化贡献累加成 master 基上的稀疏留数（镜像 residue_entries）。
    // 缓存成员声明为 mutable：重建结果只依赖构造期输入，const 调用可安全写入。
    void residue_entries(
        const std::vector<TensorContributionData>& contributions,
        std::map<std::pair<long, long>, GaussVal>& entries) const {
        GaussVal zero;
        for (const auto& contribution : contributions) {
            if (contribution.kind == TensorKind::Diagonal) {
                long offset = 0;
                if (contribution.scope == TensorScope::Sector) {
                    const size_t index =
                        sector_index(contribution.sector_key, "tensor diagonal");
                    offset = sectors_[index].first_global_index();
                }
                for (size_t position = 0; position < contribution.vector.size();
                     ++position) {
                    const GaussVal& value = contribution.vector[position];
                    if (value.is_zero()) continue;
                    const long index = offset + static_cast<long>(position);
                    const std::pair<long, long> position_key(index, index);
                    GaussVal sum;
                    const auto existing = entries.find(position_key);
                    gauss_add(sum, existing == entries.end() ? zero : existing->second,
                              value);
                    entries[position_key] = sum;
                }
                continue;
            }
            const size_t index =
                sector_index(contribution.sector_key, "tensor conjugated");
            const TensorSectorData& sector = sectors_[index];
            const auto& block = conjugated_block(sector.key, contribution.letter_indicator);
            const long offset = sector.first_global_index();
            if (contribution.scope == TensorScope::Sector) {
                for (size_t column = 0; column < contribution.pre_vector.size(); ++column) {
                    GaussVal factor;
                    gauss_mul(factor, contribution.post_factor,
                              contribution.pre_vector[column]);
                    if (factor.is_zero()) continue;
                    for (long row = 0; row < sector.dimension; ++row) {
                        const GaussVal& value =
                            block[static_cast<size_t>(row)][column];
                        if (value.is_zero()) continue;
                        const std::pair<long, long> position_key(
                            offset + row, offset + static_cast<long>(column));
                        GaussVal product;
                        gauss_mul(product, factor, value);
                        GaussVal sum;
                        const auto existing = entries.find(position_key);
                        gauss_add(sum,
                                  existing == entries.end() ? zero : existing->second,
                                  product);
                        entries[position_key] = sum;
                    }
                }
            } else {
                for (size_t item = 0; item < contribution.middle_values.size(); ++item) {
                    GaussVal factor;
                    gauss_mul(factor, contribution.post_factor,
                              contribution.middle_values[item]);
                    if (factor.is_zero()) continue;
                    const long local_row = contribution.middle_rows[item];
                    const long column = contribution.middle_columns[item];
                    for (long row = 0; row < sector.dimension; ++row) {
                        const GaussVal& value =
                            block[static_cast<size_t>(row)][static_cast<size_t>(local_row)];
                        if (value.is_zero()) continue;
                        const std::pair<long, long> position_key(offset + row, column);
                        GaussVal product;
                        gauss_mul(product, factor, value);
                        GaussVal sum;
                        const auto existing = entries.find(position_key);
                        gauss_add(sum,
                                  existing == entries.end() ? zero : existing->second,
                                  product);
                        entries[position_key] = sum;
                    }
                }
            }
        }
        for (auto item = entries.begin(); item != entries.end();) {
            if (item->second.is_zero()) {
                item = entries.erase(item);
            } else {
                ++item;
            }
        }
    }

  private:
    std::vector<TensorSectorData> sectors_;
    long dimension_ = 0;
    std::map<std::string, size_t> by_key_;
    std::map<std::string, std::vector<std::vector<GaussVal>>> mutable u_;
    std::map<std::string, std::vector<std::vector<GaussVal>>> mutable u_inverse_;
    std::map<std::pair<std::string, std::vector<long>>,
             std::vector<std::vector<GaussVal>>>
        mutable blocks_;
};

// ---------- 请求侧解析与校验（镜像适配器 flintnde_transport 的检查） ----------

inline void require_exact_keys_tensor(
    const fjson::Value& value, const std::set<std::string>& expected,
    const std::string& label) {
    if (!value.is_object()) {
        throw std::runtime_error(label + " must be a JSON object");
    }
    std::set<std::string> actual;
    for (const auto& entry : value.fields) actual.insert(entry.first);
    if (actual != expected) {
        std::string message = label + " fields do not match the current schema";
        for (const auto& key : expected) {
            if (!actual.count(key)) message += "; missing=" + key;
        }
        for (const auto& key : actual) {
            if (!expected.count(key)) message += "; unexpected=" + key;
        }
        throw std::runtime_error(message);
    }
}

inline GaussVal parse_gauss_val(const fjson::Value& record, const std::string& label) {
    require_exact_keys_tensor(record, {"real", "imag"}, label);
    if (!record.at("real").is_string() || !record.at("imag").is_string()) {
        throw std::runtime_error(label + " real/imag must be strings");
    }
    GaussQ parsed;
    parse_gauss_record_exact(record.at("real").text, record.at("imag").text, parsed,
                             label.c_str());
    GaussVal value;
    fmpq_set(value.re, parsed.re);
    fmpq_set(value.im, parsed.im);
    return value;
}

inline long tensor_integer(const fjson::Value& value, const std::string& label) {
    return static_cast<long>(value.as_integer(label.c_str()));
}

inline TensorReconstructor parse_tensor_sectors(
    const fjson::Value& tensor_sectors, long dimension) {
    if (!tensor_sectors.is_array() || tensor_sectors.items.empty()) {
        throw std::runtime_error("tensorSectors must be a nonempty list");
    }
    static const std::set<std::string> known_bases = {"massiveT", "masslessH", "identity"};
    std::vector<TensorSectorData> sectors;
    std::set<std::string> keys;
    std::vector<std::pair<long, long>> offsets;
    for (size_t index = 0; index < tensor_sectors.items.size(); ++index) {
        const std::string label = "tensorSectors[" + std::to_string(index) + "]";
        const fjson::Value& record = tensor_sectors.items[index];
        require_exact_keys_tensor(
            record,
            {"sectorKey", "dimension", "offset", "rawCount", "slotBases",
             "rawToMasterColumn", "masterRepresentativeRaw"},
            label);
        TensorSectorData sector;
        const fjson::Value& key_record = record.at("sectorKey");
        // 无根传播子的 sector 键是空串，属合法取值（镜像适配器 _sector_key）。
        if (!key_record.is_string()) {
            throw std::runtime_error(label + ".sectorKey must be a string");
        }
        sector.key = key_record.text;
        if (!keys.insert(sector.key).second) {
            throw std::runtime_error(label + ": duplicate sectorKey " + sector.key);
        }
        sector.dimension = tensor_integer(record.at("dimension"), label + ".dimension");
        sector.offset = tensor_integer(record.at("offset"), label + ".offset");
        sector.raw_count = tensor_integer(record.at("rawCount"), label + ".rawCount");
        if (sector.dimension < 1 || sector.offset < 1 || sector.raw_count < 1) {
            throw std::runtime_error(
                label + ": dimension/offset/rawCount must be positive");
        }
        const fjson::Value& bases = record.at("slotBases");
        if (!bases.is_array()) {
            throw std::runtime_error(label + ".slotBases must be a list");
        }
        for (const auto& base : bases.items) {
            if (!base.is_string() || !known_bases.count(base.text)) {
                throw std::runtime_error(
                    label + ".slotBases must be known local base tags");
            }
            sector.slot_bases.push_back(base.text);
        }
        if (sector.raw_count != (1L << sector.slot_bases.size())) {
            throw std::runtime_error(
                label + ".rawCount must equal 2**len(slotBases)");
        }
        const fjson::Value& columns = record.at("rawToMasterColumn");
        if (!columns.is_array() ||
            columns.items.size() != static_cast<size_t>(sector.raw_count)) {
            throw std::runtime_error(label + ".rawToMasterColumn length mismatch");
        }
        for (const auto& column : columns.items) {
            const long value =
                tensor_integer(column, label + ".rawToMasterColumn[]");
            if (value < 0 || value > sector.dimension) {
                throw std::runtime_error(
                    label + ".rawToMasterColumn[] out of range");
            }
            sector.raw_to_master_column.push_back(value);
        }
        const fjson::Value& representatives = record.at("masterRepresentativeRaw");
        if (!representatives.is_array() ||
            representatives.items.size() != static_cast<size_t>(sector.dimension)) {
            throw std::runtime_error(
                label + ".masterRepresentativeRaw length mismatch");
        }
        for (size_t position = 0; position < representatives.items.size(); ++position) {
            const long raw = tensor_integer(
                representatives.items[position],
                label + ".masterRepresentativeRaw[]");
            if (raw < 1 || raw > sector.raw_count) {
                throw std::runtime_error(
                    label + ".masterRepresentativeRaw[] out of range");
            }
            if (sector.raw_to_master_column[static_cast<size_t>(raw - 1)] !=
                static_cast<long>(position) + 1) {
                throw std::runtime_error(
                    label + ": raw representative does not map to its master column");
            }
            sector.master_representative_raw.push_back(raw);
        }
        offsets.emplace_back(sector.offset, sector.dimension);
        sectors.push_back(std::move(sector));
    }
    return TensorReconstructor(std::move(sectors), dimension);
}

inline std::vector<TensorContributionData> parse_tensor_contributions(
    const fjson::Value& value, const TensorReconstructor& tensor,
    const std::string& label) {
    if (!value.is_array() || value.items.empty()) {
        throw std::runtime_error(label + " must be a nonempty list");
    }
    std::vector<TensorContributionData> contributions;
    for (size_t index = 0; index < value.items.size(); ++index) {
        const std::string item_label = label + "[" + std::to_string(index) + "]";
        const fjson::Value& record = value.items[index];
        if (!record.is_object()) {
            throw std::runtime_error(item_label + " must be a JSON object");
        }
        const fjson::Value* kind = record.find("kind");
        TensorContributionData contribution;
        if (kind != nullptr && kind->is_string() && kind->text == "diagonal") {
            require_exact_keys_tensor(
                record, {"kind", "scope", "sectorKey", "vector"}, item_label);
            contribution.kind = TensorKind::Diagonal;
            const fjson::Value& scope = record.at("scope");
            long length = 0;
            if (scope.is_string() && scope.text == "sector") {
                contribution.scope = TensorScope::Sector;
                const fjson::Value& key = record.at("sectorKey");
                if (!key.is_string()) {
                    throw std::runtime_error(
                        item_label + ".sectorKey must be a string");
                }
                contribution.sector_key = key.text;
                length = static_cast<long>(tensor.sectors()[
                    tensor.sector_index(contribution.sector_key, item_label)]
                                               .dimension);
            } else if (scope.is_string() && scope.text == "global") {
                contribution.scope = TensorScope::Global;
                if (!record.at("sectorKey").is_null()) {
                    throw std::runtime_error(
                        item_label + ".sectorKey must be null for global scope");
                }
                length = tensor.dimension();
            } else {
                throw std::runtime_error(
                    item_label + ".scope must be \"sector\" or \"global\"");
            }
            const fjson::Value& vector = record.at("vector");
            if (!vector.is_array() ||
                vector.items.size() != static_cast<size_t>(length)) {
                throw std::runtime_error(item_label + ".vector length mismatch");
            }
            for (size_t position = 0; position < vector.items.size(); ++position) {
                contribution.vector.push_back(parse_gauss_val(
                    vector.items[position],
                    item_label + ".vector[" + std::to_string(position) + "]"));
            }
            contributions.push_back(std::move(contribution));
            continue;
        }
        if (kind == nullptr || !kind->is_string() || kind->text != "conjugated") {
            throw std::runtime_error(
                item_label + ".kind must be \"diagonal\" or \"conjugated\"");
        }
        require_exact_keys_tensor(
            record,
            {"kind", "scope", "sectorKey", "preVector", "letterIndicator", "middle",
             "postFactor"},
            item_label);
        contribution.kind = TensorKind::Conjugated;
        const fjson::Value& key = record.at("sectorKey");
        if (!key.is_string()) {
            throw std::runtime_error(item_label + ".sectorKey must be a string");
        }
        contribution.sector_key = key.text;
        const long local = static_cast<long>(tensor.sectors()[
            tensor.sector_index(contribution.sector_key, item_label)].dimension);
        const fjson::Value& indicator = record.at("letterIndicator");
        if (!indicator.is_array() ||
            indicator.items.size() != static_cast<size_t>(local)) {
            throw std::runtime_error(item_label + ".letterIndicator length mismatch");
        }
        for (const auto& bit : indicator.items) {
            const long value = tensor_integer(bit, item_label + ".letterIndicator[]");
            if (value != 0 && value != 1) {
                throw std::runtime_error(
                    item_label + ".letterIndicator[] must be 0 or 1");
            }
            contribution.letter_indicator.push_back(value);
        }
        contribution.post_factor =
            parse_gauss_val(record.at("postFactor"), item_label + ".postFactor");
        const fjson::Value& scope = record.at("scope");
        if (scope.is_string() && scope.text == "sector") {
            contribution.scope = TensorScope::Sector;
            if (record.at("preVector").is_null()) {
                throw std::runtime_error(
                    item_label + ".preVector is required for sector scope");
            }
            const fjson::Value& pre_vector = record.at("preVector");
            if (!pre_vector.is_array() ||
                pre_vector.items.size() != static_cast<size_t>(local)) {
                throw std::runtime_error(item_label + ".preVector length mismatch");
            }
            for (size_t position = 0; position < pre_vector.items.size(); ++position) {
                contribution.pre_vector.push_back(parse_gauss_val(
                    pre_vector.items[position],
                    item_label + ".preVector[" + std::to_string(position) + "]"));
            }
            if (!record.at("middle").is_null()) {
                throw std::runtime_error(
                    item_label + ".middle must be null for sector scope");
            }
        } else if (scope.is_string() && scope.text == "global") {
            contribution.scope = TensorScope::Global;
            if (!record.at("preVector").is_null()) {
                throw std::runtime_error(
                    item_label + ".preVector must be null for global scope");
            }
            if (record.at("middle").is_null()) {
                throw std::runtime_error(
                    item_label + ".middle is required for global scope");
            }
            const fjson::Value& middle = record.at("middle");
            require_exact_keys_tensor(middle, {"rows", "columns", "entries"},
                                      item_label + ".middle");
            const fjson::Value& rows = middle.at("rows");
            const fjson::Value& columns = middle.at("columns");
            const fjson::Value& entries = middle.at("entries");
            if (!rows.is_array() || !columns.is_array() || !entries.is_array() ||
                rows.items.size() != entries.items.size() ||
                columns.items.size() != entries.items.size()) {
                throw std::runtime_error(
                    item_label + ".middle sparse fields must have equal length");
            }
            for (size_t position = 0; position < entries.items.size(); ++position) {
                const long row = tensor_integer(
                    rows.items[position], item_label + ".middle.rows[]");
                if (row < 1 || row > local) {
                    throw std::runtime_error(
                        item_label + ".middle.rows[] out of range");
                }
                const long column = tensor_integer(
                    columns.items[position], item_label + ".middle.columns[]");
                if (column < 1 || column > tensor.dimension()) {
                    throw std::runtime_error(
                        item_label + ".middle.columns[] out of range");
                }
                contribution.middle_rows.push_back(row - 1);
                contribution.middle_columns.push_back(column - 1);
                contribution.middle_values.push_back(parse_gauss_val(
                    entries.items[position],
                    item_label + ".middle.entries[" + std::to_string(position) + "]"));
            }
        } else {
            throw std::runtime_error(
                item_label + ".scope must be \"sector\" or \"global\"");
        }
        contributions.push_back(std::move(contribution));
    }
    return contributions;
}

}  // namespace fnde

#endif  // FLINTNDE_ENGINE_TENSOR_HPP
