//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int/big_int.h"
#include "big_int/big_int_ops.h"

bint_t::bint_t() : sign(false) {}

bint_t::bint_t(bint_t&& other) noexcept : sign(other.sign), data(std::move(other.data)) {}

bint_t::bint_t(const uint64_t value) : sign(false) {
    data.push_back(value);
}

bint_t::bint_t(const int64_t value) : sign(value < 0) {
    data.push_back(abs(value));
}

bint_t::bint_t(const int value) : bint_t(static_cast<int64_t>(value)) {}

bint_t& bint_t::pow(const int n) {
    big_int::fast_pow_inplace(*this, n);
    return *this;
}

bint_t bint_t::pow(int n) const {
    auto result = *this;
    result.pow(n);
    return result;
}

std::strong_ordering bint_t::operator<=>(const bint_t& other) const {
    return big_int::compare(*this, other);
}

bool bint_t::operator==(const bint_t& other) const {
    return (*this <=> other) == 0;
}

bint_t bint_t::operator-() const {
    auto result = *this;
    result.sign = !result.sign;
    return result;
}

bint_t bint_t::operator+(const bint_t& other) const {
    return big_int::add(*this, other);
}

bint_t bint_t::operator-(const bint_t& other) const {
    return big_int::sub(*this, other);
}

bint_t bint_t::operator*(const bint_t& other) const {
    return big_int::multiply(*this, other);
}

bint_t bint_t::operator/(const bint_t& other) const {
    bint_t a = *this;
    a /= other;
    return a;
}

bint_t bint_t::operator%(const bint_t& other) const {
    bint_t a = *this;
    bint_t rem;
    big_int::divide_abs(a, other, rem);
    if (rem == bint_t(0ll)) return rem;
    if (this->sign != other.sign) {
        const bint_t r = rem;
        rem = other;
        big_int_impl::sub_abs_inplace(rem, r);
    }
    rem.sign = other.sign;
    return rem;
}

bint_t bint_t::operator<<(const int64_t n) const {
    bint_t result = *this;
    result <<= n;
    return result;
}

bint_t bint_t::operator>>(const int64_t n) const {
    bint_t result = *this;
    result >>= n;
    return result;
}

bint_t& bint_t::operator+=(const bint_t& other) {
    return *this = *this + other;
}

bint_t& bint_t::operator-=(const bint_t& other) {
    return *this = *this - other;
}

bint_t& bint_t::operator*=(const bint_t& other) {
    if (this == &other) {
        big_int_impl::square(*this);
        return *this;
    }

    return *this = *this * other;
}

bint_t& bint_t::operator/=(const bint_t& other) {
    bint_t rem;
    big_int::divide_abs(*this, other, rem);
    this->sign = this->sign ^ other.sign;
    if (this->sign && rem != bint_t(0ll)) {
        big_int_impl::add_abs_inplace(*this, bint_t(1ll));
    }
    return *this;
}

bint_t& bint_t::operator%=(const bint_t& other) {
    bint_t rem;
    big_int::divide_abs(*this, other, rem);
    if (rem == bint_t(0ll)) {
        *this = bint_t(0);
        return *this;
    }

    if (this->sign != other.sign) {
        *this = other;
        big_int_impl::sub_abs_inplace(*this, rem);
    } else {
        *this = rem;
    }
    this->sign = other.sign;
    return *this;
}

bint_t& bint_t::operator<<=(const int64_t n) {
    big_int::shift_left_inplace(*this, n);
    return *this;
}

bint_t& bint_t::operator>>=(const int64_t n) {
    big_int::shift_right_inplace(*this, n);
    return *this;
}

void small_to_string(const bint_t& a, std::string& buffer, int digits) {
    std::string result;
    bint_t current = a;
    uint64_t rem;
    if (current.data.empty() || current.data.size() == 1 && current.data[0] == 0ull) { // if current == 0
        result.push_back('0');
    }
    while (current.data.size() > 1 || (!current.data.empty() && current.data[0] != 0ull)) { // while current != 0
        constexpr uint64_t divisor = 10000000000000000000ull;
        big_int_impl::div_abs_inplace(current, divisor, rem);
        for (int _ = 0; _ < 19; _++) {
            result.push_back("0123456789"[rem % 10]);
            rem /= 10;
        }
    }
    while (result.size() > 1 && result.back() == '0') {
        result.pop_back();
    }
    if (digits > 0) {
        if (result.size() > digits) {
            throw std::runtime_error("small_to_string: result.size()[" + std::to_string(result.size()) +
                "] > digits[" + std::to_string(digits) + "]: " + result);
        }
        for (int i = 0; i < digits - result.size(); i++) {
            buffer.push_back('0');
        }
    }
    std::ranges::reverse(result);
    buffer.append(result);
}

void to_string(const bint_t& a, std::string& buffer, int digits) { // NOLINT(*-no-recursion)
    if (a.data.size() < 10)
        return small_to_string(a, buffer, digits);

    const int half_length10 = static_cast<int>(a.data.size() * 32 / std::log2(10.0));
    bint_t big10(10ll);
    big_int::fast_pow_inplace(big10, half_length10);

    bint_t high = a;
    bint_t low;
    big_int_impl::normalize(high);
    big_int::divide_abs(high, big10, low);
    to_string(high, buffer, digits - half_length10);
    to_string(low, buffer, half_length10);
}

std::string bint_t::to_string() const {
    std::string result;
    const int expected_length = this->data.size() * 64 / std::log2(10.0) + 1;
    result.reserve(expected_length + (this->sign ? 1 : 0));

    if (this->sign) result.push_back('-');
    ::to_string(*this, result, 0);
    return result;
}


