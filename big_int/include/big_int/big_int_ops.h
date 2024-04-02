//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#ifndef BIG_INT_OPS_H
#define BIG_INT_OPS_H

#include "big_int.h"

namespace big_int {
    constexpr int KARATSUBA_THRESHOLD = 50;
    constexpr int KARATSUBA_SQUARE_THRESHOLD = 50;
    constexpr int TOOM_COOK_THRESHOLD = 220;
    constexpr int TOOM_COOK_SQUARE_THRESHOLD = 220;

    bint_t add(const bint_t& a, const bint_t& b);
    bint_t sub(const bint_t& a, const bint_t& b);
    bint_t multiply(const bint_t& a, const bint_t& b, int a_limit = -1, int b_limit = -1);
    void divide_abs(bint_t& a, const bint_t& b, bint_t& rem);

    void divide_knuth_abs(bint_t& a, const bint_t& b, bint_t& rem);
    void fast_pow_inplace(bint_t& a, uint64_t n);
    void shift_left_inplace(bint_t& a, int64_t shift);
    void shift_right_inplace(bint_t& a, int64_t shift);

    std::strong_ordering compare(const bint_t& a, const bint_t& b);
}

namespace big_int_impl {
    void square(bint_t& a, int a_limit = -1);
    bint_t schoolbook_square(const bint_t& a, int a_limit);
    bint_t schoolbook_multiply(const bint_t& a, const bint_t& b, int a_limit = -1, int b_limit = -1);
    std::strong_ordering compare_abs(const bint_t&a, const bint_t& b, int shift = 0);
    void add_abs_inplace(bint_t& a, const bint_t& b, int b_start = 0, int b_limit = -1, int b_shift = 0, bool let_overflow = false);
    void sub_abs_inplace(bint_t& a, const bint_t& b, int b_limit = -1);
    void div_abs_inplace(bint_t& a, uint64_t b, uint64_t& rem);
    uint64_t count_bits(const bint_t& a);
    void normalize(bint_t& a);
    bool is_normalized(const bint_t& a);
}

#endif //BIG_INT_OPS_H
