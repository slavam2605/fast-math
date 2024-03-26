//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int/karatsuba.h"
#include "big_int/big_int_ops.h"

using namespace big_int_impl;

std::pair<bint_t, int> split_abs(const bint_t& a, int m, int a_limit) {
    if (a_limit <= m)
        return std::make_pair(bint_t(0ll), a_limit);

    bint_t high;
    high.data.resize(a_limit - m);
    std::copy_n(a.data.begin() + m, a_limit - m, high.data.begin());
    return std::make_pair(high, m);
}

bint_t karatsuba(const bint_t& a, const bint_t& b, int a_limit, int b_limit) { // NOLINT(*-no-recursion)
    if (a_limit < 0) a_limit = a.data.size();
    if (b_limit < 0) b_limit = b.data.size();

    if (a_limit <= big_int::KARATSUBA_THRESHOLD || b_limit <= big_int::KARATSUBA_THRESHOLD)
        return schoolbook_multiply(a, b, a_limit, b_limit);

    int m = (std::max(a_limit, b_limit) + 1) / 2;
    auto [a1, a0_limit] = split_abs(a, m, a_limit);
    auto [b1, b0_limit] = split_abs(b, m, b_limit);

    auto z2 = karatsuba(a1, b1);
    auto z0 = karatsuba(a, b, a0_limit, b0_limit);
    add_abs_inplace(a1, a, a0_limit);
    add_abs_inplace(b1, b, b0_limit);
    auto z1 = karatsuba(a1, b1);
    sub_abs_inplace(z1, z2);
    sub_abs_inplace(z1, z0);

    add_abs_inplace(z0, z1, -1, m);
    add_abs_inplace(z0, z2, -1, 2 * m);
    z0.sign = a.sign ^ b.sign;
    normalize(z0);

    return z0;
}
