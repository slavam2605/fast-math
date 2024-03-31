//
// Created by Vyacheslav.Moklev on 24/03/2024.
//

#include "big_int/toom_cook.h"
#include "big_int/big_int_ops.h"

using namespace big_int;
using namespace big_int_impl;

std::tuple<bint_t, bint_t, bint_t> split_abs_3way(const bint_t& a, const int a_limit, const int m) {
    if (a_limit <= m) {
        bint_t copy;
        copy.data.resize(a_limit);
        std::copy_n(a.data.begin(), a_limit, copy.data.begin());
        return std::make_tuple(bint_t(0ll), bint_t(0ll), copy);
    }

    bint_t low;
    low.data.resize(m);
    std::copy_n(a.data.begin(), m, low.data.begin());
    normalize(low);

    bint_t mid;
    mid.data.resize(std::min(m, a_limit - m));
    std::copy_n(a.data.begin() + m, mid.data.size(), mid.data.begin());
    normalize(mid);

    if (a_limit <= 2 * m)
        return std::make_tuple(bint_t(0ll), mid, low);

    bint_t high;
    high.data.resize(std::min(m, a_limit - 2 * m));
    std::copy_n(a.data.begin() + 2 * m, high.data.size(), high.data.begin());
    normalize(high);

    return std::make_tuple(high, mid, low);
}

std::pair<bint_t, bint_t> split_abs_3way_inplace(bint_t& a, const int a_limit, const int m) {
    if (a_limit <= 2 * m) {
        throw std::runtime_error("split_abs_3way_inplace: a_limit [" + std::to_string(a_limit) +
            "] <= 2 * m [" + std::to_string(2 * m) + "]");
    }

    bint_t mid;
    mid.data.resize(std::min(m, a_limit - m));
    std::copy_n(a.data.begin() + m, mid.data.size(), mid.data.begin());
    normalize(mid);

    bint_t high;
    high.data.resize(std::min(m, a_limit - 2 * m));
    std::copy_n(a.data.begin() + 2 * m, high.data.size(), high.data.begin());
    normalize(high);

    a.sign = false;
    a.data.resize(m);
    normalize(a);
    return std::make_pair(high, mid);
}

void div_abs_inplace_3_assert_rem(bint_t& a) {
    __uint128_t current = 0;
    for (int i = a.data.size() - 1; i >= 0; i--) {
        current = (current << 64) + a.data[i];
        a.data[i] = current / 3;
        current %= 3;
    }
    while (a.data.size() > 1 && a.data.back() == 0) {
        a.data.pop_back();
    }
    if (current != 0) throw std::runtime_error("div_abs_inplace_3_assert_rem: rem was " +
        std::to_string(static_cast<uint64_t>(current)) + ", a = " + a.to_string());
}

bint_t toom3(const bint_t& a, const bint_t& b, int a_limit, int b_limit) {
    if (a_limit < 0) a_limit = a.data.size();
    if (b_limit < 0) b_limit = b.data.size();

    int m = (std::max(a_limit, b_limit) + 2) / 3;
    auto [a2, a1, a0] = split_abs_3way(a, a_limit, m);
    auto [b2, b1, b0] = split_abs_3way(b, b_limit, m);

    auto p0 = a0 + a2;
    auto p_1 = p0 + a1;
    auto p_m1 = p0 - a1;
    auto p_m2 = p_m1 + a2;
    p_m2 <<= 1;
    p_m2 = p_m2 - a0;

    auto q0 = b0 + b2;
    auto q_1 = q0 + b1;
    auto q_m1 = q0 - b1;
    auto q_m2 = q_m1 + b2;
    q_m2 <<= 1;
    q_m2 = q_m2 - b0;

    auto r_0 = multiply(a0, b0);
    auto r_1 = multiply(p_1, q_1);
    auto r_m1 = multiply(p_m1, q_m1);
    auto r_m2 = multiply(p_m2, q_m2);
    auto r_inf = multiply(a2, b2);

    auto r4 = r_inf;
    auto r3 = r_m2 - r_1;
    div_abs_inplace_3_assert_rem(r3);
    auto r1 = r_1 - r_m1;
    r1 >>= 1;
    auto r2 = r_m1 - r_0;
    auto diff = r2 - r3;
    r3 = r2 - r3;
    r3 >>= 1;
    r_inf <<= 1;
    r3 = r3 + r_inf;
    r2 = r2 + r1 - r4;
    r1 = r1 - r3;

    add_abs_inplace(r_0, r1, 0, -1, m);
    add_abs_inplace(r_0, r2, 0, -1, 2 * m);
    add_abs_inplace(r_0, r3, 0, -1, 3 * m);
    add_abs_inplace(r_0, r4, 0, -1, 4 * m);
    r_0.sign = a.sign ^ b.sign;
    normalize(r_0);
    return r_0;
}

void toom3_square(bint_t& a, int a_limit) {
    if (a_limit < 0) a_limit = a.data.size();

    const int m = (a_limit + 2) / 3;
    auto [a2, a1] = split_abs_3way_inplace(a, a_limit, m);

    const auto p0 = a + a2;
    auto p_1 = p0 + a1;
    auto p_m1 = p0 - a1;
    auto p_m2 = p_m1 + a2;
    p_m2 <<= 1;
    p_m2 = p_m2 - a;

    square(a);
    square(p_1);
    square(p_m1);
    square(p_m2);
    square(a2);

    const auto r4 = a2;
    auto r3 = p_m2 - p_1;
    div_abs_inplace_3_assert_rem(r3);
    auto r1 = p_1 - p_m1;
    r1 >>= 1;
    auto r2 = p_m1 - a;
    auto diff = r2 - r3;
    r3 = r2 - r3;
    r3 >>= 1;
    a2 <<= 1;
    r3 = r3 + a2;
    r2 = r2 + r1 - r4;
    r1 = r1 - r3;

    add_abs_inplace(a, r1, 0, -1, m);
    add_abs_inplace(a, r2, 0, -1, 2 * m);
    add_abs_inplace(a, r3, 0, -1, 3 * m);
    add_abs_inplace(a, r4, 0, -1, 4 * m);
    normalize(a); // a.sign is already false from split_abs_3way_inplace
}