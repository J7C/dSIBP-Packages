// FlintNDE 引擎数值工具：FLINT RAII 包装、Wolfram 数值文本规整、letter 部件的
// 有理/十进制双态解析与精度换算。
#ifndef FLINTNDE_ENGINE_UTIL_HPP
#define FLINTNDE_ENGINE_UTIL_HPP

#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "flint/acb.h"
#include "flint/acb_mat.h"
#include "flint/arb.h"
#include "flint/fmpq.h"
#include "flint/fmpz.h"

namespace fnde {

// 与 flintnde.core.configure_working_precision 相同的换算：位数 + 32 guard bits。
inline long working_bits(long digits) {
    return static_cast<long>(std::ceil(digits * std::log2(10.0))) + 32;
}

// FLINT C API 无全局精度：FlintNDE 引擎用线程局部精度镜像 python ctx.prec 语义。
inline long& prec_slot() {
    thread_local long precision = 128;
    return precision;
}
inline void set_thread_prec(long bits) { prec_slot() = bits; }
inline long prec() { return prec_slot(); }

struct Acb {
    acb_t value;
    Acb() { acb_init(value); }
    explicit Acb(const acb_t source) {
        acb_init(value);
        acb_set(value, source);
    }
    ~Acb() { acb_clear(value); }
    Acb(const Acb&) = delete;
    Acb& operator=(const Acb&) = delete;
    operator acb_ptr() { return value; }
    operator acb_srcptr() const { return value; }
};

struct Arb {
    arb_t value;
    Arb() { arb_init(value); }
    ~Arb() { arb_clear(value); }
    Arb(const Arb&) = delete;
    Arb& operator=(const Arb&) = delete;
    operator arb_ptr() { return value; }
    operator arb_srcptr() const { return value; }
};

struct AcbMat {
    acb_mat_t matrix;
    explicit AcbMat(long rows, long cols) { acb_mat_init(matrix, rows, cols); }
    explicit AcbMat(const acb_mat_t source) {
        acb_mat_init(matrix, acb_mat_nrows(source), acb_mat_ncols(source));
        acb_mat_set(matrix, source);
    }
    ~AcbMat() { acb_mat_clear(matrix); }
    AcbMat(const AcbMat&) = delete;
    AcbMat& operator=(const AcbMat&) = delete;
    long nrows() const { return acb_mat_nrows(matrix); }
    long ncols() const { return acb_mat_nrows(matrix) == 0 ? 0 : acb_mat_ncols(matrix); }
    acb_ptr at(long row, long col) { return acb_mat_entry(matrix, row, col); }
    acb_srcptr at(long row, long col) const { return acb_mat_entry(matrix, row, col); }
};

// Q(i) 有理数对：极点合并与 exact 判等需要保留 fmpq 层。
struct GaussQ {
    fmpq_t re;
    fmpq_t im;
    GaussQ() {
        fmpq_init(re);
        fmpq_init(im);
    }
    ~GaussQ() {
        fmpq_clear(re);
        fmpq_clear(im);
    }
    GaussQ(const GaussQ&) = delete;
    GaussQ& operator=(const GaussQ&) = delete;
    bool equals(const GaussQ& other) const {
        return fmpq_equal(re, other.re) && fmpq_equal(im, other.im);
    }
};

// 可拷贝的 Q(i) 值：张量重建需要在容器里存放精确值（GaussQ 禁拷贝）。
struct GaussVal {
    fmpq_t re;
    fmpq_t im;
    GaussVal() {
        fmpq_init(re);
        fmpq_init(im);
    }
    ~GaussVal() {
        fmpq_clear(re);
        fmpq_clear(im);
    }
    GaussVal(const GaussVal& other) {
        fmpq_init(re);
        fmpq_init(im);
        fmpq_set(re, other.re);
        fmpq_set(im, other.im);
    }
    GaussVal& operator=(const GaussVal& other) {
        if (this != &other) {
            fmpq_set(re, other.re);
            fmpq_set(im, other.im);
        }
        return *this;
    }
    GaussVal(GaussVal&& other) {
        fmpq_init(re);
        fmpq_init(im);
        fmpq_set(re, other.re);
        fmpq_set(im, other.im);
    }
    GaussVal& operator=(GaussVal&& other) {
        if (this != &other) {
            fmpq_set(re, other.re);
            fmpq_set(im, other.im);
        }
        return *this;
    }
    static GaussVal rational(long long numerator, long long denominator = 1) {
        GaussVal value;
        fmpq_set_si(value.re, numerator, denominator);
        fmpq_zero(value.im);
        return value;
    }
    static GaussVal gaussian(long long re_num, long long re_den, long long im_num,
                             long long im_den) {
        GaussVal value;
        fmpq_set_si(value.re, re_num, re_den);
        fmpq_set_si(value.im, im_num, im_den);
        return value;
    }
    bool is_zero() const {
        return fmpq_is_zero(re) && fmpq_is_zero(im);
    }
};

inline void gauss_add(GaussVal& out, const GaussVal& a, const GaussVal& b) {
    fmpq_add(out.re, a.re, b.re);
    fmpq_add(out.im, a.im, b.im);
}

inline void gauss_sub(GaussVal& out, const GaussVal& a, const GaussVal& b) {
    fmpq_sub(out.re, a.re, b.re);
    fmpq_sub(out.im, a.im, b.im);
}

inline void gauss_mul(GaussVal& out, const GaussVal& a, const GaussVal& b) {
    fmpq_t real, imaginary;
    fmpq_init(real);
    fmpq_init(imaginary);
    fmpq_mul(real, a.re, b.re);
    fmpq_t cross;
    fmpq_init(cross);
    fmpq_mul(cross, a.im, b.im);
    fmpq_sub(real, real, cross);
    fmpq_mul(imaginary, a.re, b.im);
    fmpq_mul(cross, a.im, b.re);
    fmpq_add(imaginary, imaginary, cross);
    fmpq_set(out.re, real);
    fmpq_set(out.im, imaginary);
    fmpq_clear(real);
    fmpq_clear(imaginary);
    fmpq_clear(cross);
}

// 十进制串（可带 e 指数）精确解析为 fmpq，语义对齐 Python Fraction(text)：
// "1.5e-3" -> 3/2000。输入需已经过 clean_wolfram_number_text 规整。
inline bool parse_exact_decimal_fmpq(const std::string& cleaned, fmpq_t out) {
    if (cleaned.empty()) return false;
    size_t exponent_position = std::string::npos;
    for (size_t index = 0; index < cleaned.size(); ++index) {
        const char character = cleaned[index];
        if (character == 'e' || character == 'E') {
            if (exponent_position != std::string::npos) return false;
            exponent_position = index;
        }
    }
    std::string mantissa = cleaned;
    long exponent = 0;
    if (exponent_position != std::string::npos) {
        mantissa = cleaned.substr(0, exponent_position);
        const std::string exponent_text =
            cleaned.substr(exponent_position + 1);
        if (exponent_text.empty()) return false;
        try {
            size_t consumed = 0;
            exponent = std::stol(exponent_text, &consumed);
            if (consumed != exponent_text.size()) return false;
        } catch (const std::exception&) {
            return false;
        }
    }
    size_t start = 0;
    if (!mantissa.empty() && (mantissa[0] == '+' || mantissa[0] == '-')) ++start;
    if (start >= mantissa.size()) return false;
    bool negative = !mantissa.empty() && mantissa[0] == '-';
    std::string digits;
    long fraction_count = 0;
    bool seen_point = false;
    for (size_t index = start; index < mantissa.size(); ++index) {
        const char character = mantissa[index];
        if (character == '.') {
            if (seen_point) return false;
            seen_point = true;
            continue;
        }
        if (character < '0' || character > '9') return false;
        digits.push_back(character);
        if (seen_point) ++fraction_count;
    }
    if (digits.empty()) return false;
    fmpz_t numerator;
    fmpz_t denominator;
    fmpz_t power;
    fmpz_init(numerator);
    fmpz_init(denominator);
    fmpz_init(power);
    bool ok = true;
    if (fmpz_set_str(numerator, digits.c_str(), 10) != 0) ok = false;
    if (ok) {
        fmpz_set_ui(denominator, 1);
        const long shift = exponent - fraction_count;
        fmpz_ui_pow_ui(power, 10, static_cast<ulong>(shift >= 0 ? shift : -shift));
        if (shift >= 0) {
            fmpz_mul(numerator, numerator, power);
        } else {
            fmpz_mul(denominator, denominator, power);
        }
        if (negative) fmpz_neg(numerator, numerator);
        fmpq_set_fmpz_frac(out, numerator, denominator);
        fmpq_canonicalise(out);
    }
    fmpz_clear(numerator);
    fmpz_clear(denominator);
    fmpz_clear(power);
    return ok;
}

inline std::string clean_wolfram_number_text(const std::string& raw);

// 把 {real,imag} 记录解析为精确 Q(i)：有理串走 fmpq_set_str，十进制串走
// parse_exact_decimal_fmpq。语义对齐 Python gaussian_rational（Fraction 解析，
// 小数按精确有理处理），与 parse_letter_pair 的十进制球语义不同——张量重建
// 需要精确车道才能与 Python 车道逐位一致。
inline void parse_gauss_record_exact(const std::string& real_text,
                                     const std::string& imag_text, GaussQ& out,
                                     const char* label) {
    const std::string re_cleaned = clean_wolfram_number_text(real_text);
    const std::string im_cleaned = clean_wolfram_number_text(imag_text);
    bool ok = true;
    if (!re_cleaned.empty() && fmpq_set_str(out.re, re_cleaned.c_str(), 10) != 0) {
        ok = parse_exact_decimal_fmpq(re_cleaned, out.re);
    }
    if (ok && !im_cleaned.empty() &&
        fmpq_set_str(out.im, im_cleaned.c_str(), 10) != 0) {
        ok = parse_exact_decimal_fmpq(im_cleaned, out.im);
    }
    if (!ok) {
        throw std::runtime_error(std::string("cannot parse exact gaussian scalar (") +
                                 label + "): " + real_text + "," + imag_text);
    }
}

// Wolfram 数值文本规整：丢弃 `prec 反引号精度后缀（后随数字与点）、把 *^ 换成
// e、十进制尾部裸点补 0。有理串（"3/7"、"-5"）与普通十进制串原样通过，
// 规整结果保持与旧输入等价。
inline std::string clean_wolfram_number_text(const std::string& raw) {
    std::string cleaned;
    cleaned.reserve(raw.size());
    bool in_backtick_suffix = false;
    for (size_t index = 0; index < raw.size(); ++index) {
        const char character = raw[index];
        if (character == '`') {
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
    if (!cleaned.empty() && cleaned.back() == '.') cleaned.push_back('0');
    return cleaned;
}

// 解析 Wolfram InputForm 有理数（"3/7"、"−5"）或十进制串（可带 `prec 后缀与
// *^ 记数）；十进制按当前线程工作精度取球。失败抛异常。
inline void parse_rational_or_decimal(const std::string& text, arb_t out) {
    const std::string cleaned = clean_wolfram_number_text(text);
    if (cleaned == "0") {
        // 稀疏留数矩阵里 "0" 占绝大多数：快速路径绕过 fmpq 机制。
        arb_zero(out);
        return;
    }
    fmpq_t rational;
    fmpq_init(rational);
    if (fmpq_set_str(rational, cleaned.c_str(), 10) == 0) {
        // 与 python-flint ctx.prec 舍入一致：按当前线程工作精度取球。
        arb_set_fmpq(out, rational, prec());
    } else if (arb_set_str(out, cleaned.c_str(), prec()) != 0) {
        fmpq_clear(rational);
        throw std::runtime_error("cannot parse exact scalar string: " + text);
    }
    fmpq_clear(rational);
}

// letter alpha/beta 单个实部/虚部的双态载体：exact 为真时 rational 持有 fmpq
// 精确值；为假时该部件是十进制串。ball 在两种情形都填好（exact 时按当前精度
// 取球），供球语义判零与比较直接使用。
struct LetterPart {
    bool exact = true;
    fmpq_t rational;
    arb_t ball;
    LetterPart() {
        fmpq_init(rational);
        arb_init(ball);
    }
    ~LetterPart() {
        fmpq_clear(rational);
        arb_clear(ball);
    }
    LetterPart(const LetterPart&) = delete;
    LetterPart& operator=(const LetterPart&) = delete;
};

// 双态解析：有理串走 fmpq（exact=true，与精确 lane 语义逐位一致），十进制串
//（可带 `prec 后缀与 *^ 记数）按当前工作精度取球（exact=false）。失败抛异常。
inline void parse_letter_part(const std::string& text, LetterPart& out) {
    const std::string cleaned = clean_wolfram_number_text(text);
    if (cleaned == "0") {
        fmpq_zero(out.rational);
        arb_zero(out.ball);
        out.exact = true;
        return;
    }
    if (fmpq_set_str(out.rational, cleaned.c_str(), 10) == 0) {
        out.exact = true;
        arb_set_fmpq(out.ball, out.rational, prec());
        return;
    }
    if (arb_set_str(out.ball, cleaned.c_str(), prec()) != 0) {
        throw std::runtime_error("cannot parse letter coefficient: " + text);
    }
    out.exact = false;
}

// pole = -alpha/beta 于 Q(i)：(a+bi)/(c+di) = ((ac+bd)+(bc-ad)i)/(c^2+d^2)。
inline void gaussian_divide_pole(
    const fmpq_t a, const fmpq_t b, const fmpq_t c, const fmpq_t d, GaussQ& out) {
    fmpq_t ac, bd, bc, ad, c2, d2, denom;
    fmpq_init(ac);
    fmpq_init(bd);
    fmpq_init(bc);
    fmpq_init(ad);
    fmpq_init(c2);
    fmpq_init(d2);
    fmpq_init(denom);
    fmpq_mul(ac, a, c);
    fmpq_mul(bd, b, d);
    fmpq_mul(bc, b, c);
    fmpq_mul(ad, a, d);
    fmpq_mul(c2, c, c);
    fmpq_mul(d2, d, d);
    fmpq_add(denom, c2, d2);
    fmpq_add(out.re, ac, bd);  // -(pole) 实部分子 ac+bd；调用方再取负
    fmpq_neg(out.re, out.re);
    fmpq_sub(out.im, ad, bc);  // -(pole) 虚部分子 = ad-bc = -(bc-ad)
    fmpq_div(out.re, out.re, denom);
    fmpq_div(out.im, out.im, denom);
    fmpq_clear(ac);
    fmpq_clear(bd);
    fmpq_clear(bc);
    fmpq_clear(ad);
    fmpq_clear(c2);
    fmpq_clear(d2);
    fmpq_clear(denom);
}

// 把 Q(i) 值写入 acb（real/imag 分别取 fmpq 的球逼近）。
inline void gauss_to_acb(const GaussQ& source, acb_t out) {
    arb_t re;
    arb_t im;
    arb_init(re);
    arb_init(im);
    arb_set_fmpq(re, source.re, prec());
    arb_set_fmpq(im, source.im, prec());
    acb_set_arb_arb(out, re, im);
    arb_clear(re);
    arb_clear(im);
}

// 中点字符串：对应 python-flint str(digits, radius=False, more=True)。
inline std::string arb_midpoint_string(const arb_t value, long digits) {
    arb_t mid;
    arb_init(mid);
    arb_get_mid_arb(mid, value);
    char* text = arb_get_str(mid, digits, ARB_STR_MORE | ARB_STR_NO_RADIUS);
    arb_clear(mid);
    if (text == nullptr) throw std::runtime_error("arb_get_str failed");
    std::string result(text);
    flint_free(text);
    return result;
}

// 与 flintnde.core.arb_ball_to_json 对应的十进制球记录字段。
inline void arb_ball_fields(const arb_t value, long digits, std::string& midpoint,
                            std::string& radius, long& exponent) {
    fmpz_t m;
    fmpz_t r;
    fmpz_t e;
    fmpz_init(m);
    fmpz_init(r);
    fmpz_init(e);
    arb_get_fmpz_mid_rad_10exp(m, r, e, value, digits);
    char* mtext = fmpz_get_str(nullptr, 10, m);
    char* rtext = fmpz_get_str(nullptr, 10, r);
    midpoint = mtext;
    radius = rtext;
    exponent = static_cast<long>(fmpz_get_si(e));
    flint_free(mtext);
    flint_free(rtext);
    fmpz_clear(m);
    fmpz_clear(r);
    fmpz_clear(e);
}

// 以中点 double 做排序/比较（镜像 flintnde 的 float(x.mid()) 用法）。
inline double arb_mid_double(const arb_t value) {
    arb_t mid;
    arb_init(mid);
    arb_get_mid_arb(mid, value);
    double result = arf_get_d(arb_midref(mid), ARF_RND_NEAR);
    arb_clear(mid);
    return result;
}

}  // namespace fnde

#endif  // FLINTNDE_ENGINE_UTIL_HPP
