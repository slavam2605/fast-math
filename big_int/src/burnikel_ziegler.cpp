//
// Created by Vyacheslav.Moklev on 27/03/2024.
//

#include <big_int/burnikel_ziegler.h>
#include <big_int/big_int_ops.h>

uint64_t count_bits(const bint_t& a) {
    if (a.data.size() > 1 && a.data.back() == 0) {
        throw std::runtime_error("count_bits: value is not normalized");
    }

    return a.data.size() * 64 - std::countl_zero(a.data.back());
}

bint_t get_chunk(const bint_t& a, const int chunk_size, const int chunk_index) {
    bint_t result;
    const int leftover_size = static_cast<int>(a.data.size()) - chunk_index * chunk_size;
    if (leftover_size < 0)
        return bint_t(0);

    result.data.resize(std::min(chunk_size, leftover_size));
    std::copy(
        a.data.begin() + chunk_index * chunk_size,
        a.data.begin() + std::min(static_cast<int>(a.data.size()), (chunk_index + 1) * chunk_size),
        result.data.begin()
    );
    big_int_impl::normalize(result);
    return result;
}

void set_ones(bint_t& a, int n) {
    a.data.resize(n);
    std::ranges::fill(a.data, static_cast<uint64_t>(-1));
}

// Forward declaration for a recursion
void divide2n1n(bint_t& a, const bint_t& b, bint_t& rem);

void divide3n2n(bint_t& a, const bint_t& b, bint_t& rem) { // NOLINT(*-no-recursion)
    const int n_half = b.data.size() / 2;
    const auto a3 = get_chunk(a, n_half, 0);
    auto b1 = get_chunk(b, n_half, 1);
    a >>= n_half * 64; // now a = [a1, a2]
    if (big_int_impl::compare_abs(a, b1, n_half) < 0) {
        divide2n1n(a, b1, rem); // now rem = r1
    } else {
        rem = a;
        big_int_impl::add_abs_inplace(rem, b1);
        b1 <<= n_half * 64;
        big_int_impl::sub_abs_inplace(rem, b1);
        set_ones(a, n_half);
    }
    // now a = q_hat
    const auto d = big_int::multiply(a, b, -1, n_half);
    rem <<= n_half * 64;
    big_int_impl::add_abs_inplace(rem, a3);
    while (big_int_impl::compare_abs(rem, d) < 0) {
        rem += b;
        a -= bint_t(1);
    }
    rem -= d;
}

void divide2n1n(bint_t& a, const bint_t& b, bint_t& rem) { // NOLINT(*-no-recursion)
    // let a = [a1, a2, a3, a4] (a1 is the highest)
    // let b = [b1, b2] (b1 is the highest)
    const int n = b.data.size();
    if (n % 2 != 0 || n < BURNIKEL_ZIEGLER_THRESHOLD) {
        return big_int::divide_knuth_abs(a, b, rem);
    }

    const int n_half = n / 2;
    auto a4 = get_chunk(a, n_half, 0);
    a >>= n_half * 64; // now a = [a1, a2, a3]

    divide3n2n(a, b, rem); // now a = q1 (high part of result)
    big_int_impl::add_abs_inplace(a4, rem, 0, -1, n_half); // now a4 = [r1, r2, a4]
    divide3n2n(a4, b, rem); // now a4 = q2 (low part of result)

    a <<= n_half * 64;
    big_int_impl::add_abs_inplace(a, a4); // now a = [q1, q2]
}

void divide_burnikel_ziegler(bint_t& a, const bint_t& b, bint_t& rem) {
    const int r = a.data.size();
    const int s = b.data.size();
    if (r < s) {
        rem = a;
        a = bint_t(0);
        return;
    }

    // let m = min{2^k}, where 2^k * BURNIKEL_ZIEGLER_THRESHOLD > s
    const int m = 1 << (64 - std::countl_zero(b.data.size() / BURNIKEL_ZIEGLER_THRESHOLD));

    const int j = (s + m - 1) / m;  // j = ceil(s / m)
    const int n = j * m;            // new length for b
    const int64_t n64 = 64ll * n;

    // amount of bits to shift `a` and `b` left
    const int sigma = static_cast<int>(n64 - count_bits(b));
    const bint_t new_b = b << sigma;
    a <<= sigma;

    // amount of blocks to split `a` into, but at least 2
    const int t = std::max(2, static_cast<int>((count_bits(a) + n64) / n64));

    bint_t result(0ll);
    auto z = get_chunk(a, n, t - 2);
    big_int_impl::add_abs_inplace(z, a, n * (t - 1), std::min<int>(n * t, a.data.size()), n);
    for (int i = t - 2; i >= 0; i--) {
        divide2n1n(z, new_b, rem);
        big_int_impl::add_abs_inplace(result, z, 0, -1, n * i);
        if (i > 0) {
            z = std::move(rem);
            z <<= n64;
            big_int_impl::add_abs_inplace(z, a, n * (i - 1), n * i);
        }
    }
    rem >>= sigma;
    a = result;
    big_int_impl::normalize(a);
    big_int_impl::normalize(rem);
}