// FlintNDE 引擎数值工具：FLINT RAII 包装、Gaussian rational 解析与精度换算。
#ifndef FLINTNDE_ENGINE_UTIL_HPP
#define FLINTNDE_ENGINE_UTIL_HPP

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

// 解析 Wolfram InputForm 有理数（"3/7"、"−5"）或十进制串；失败抛异常。
inline void parse_rational_or_decimal(const std::string& text, arb_t out) {
    if (text == "0") {
        // 稀疏留数矩阵里 "0" 占绝大多数：快速路径绕过 fmpq 机制。
        arb_zero(out);
        return;
    }
    fmpq_t rational;
    fmpq_init(rational);
    if (fmpq_set_str(rational, text.c_str(), 10) == 0) {
        // 与 python-flint ctx.prec 舍入一致：按当前线程工作精度取球。
        arb_set_fmpq(out, rational, prec());
    } else if (arb_set_str(out, text.c_str(), prec()) != 0) {
        fmpq_clear(rational);
        throw std::runtime_error("cannot parse exact scalar string: " + text);
    }
    fmpq_clear(rational);
}

inline void parse_gaussian_part(const std::string& text, fmpq_t out) {
    if (fmpq_set_str(out, text.c_str(), 10) != 0) {
        throw std::runtime_error("cannot parse Gaussian rational part: " + text);
    }
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
