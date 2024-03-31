//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int/big_int.h"
#include "big_int/big_int_ops.h"

std::vector<bint_t> bint_t::power10_conversion_cache = {bint_t(10)};

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

void small_to_string(const bint_t& a, std::string& buffer, const int digits) {
    std::vector<uint64_t> groups;
    bint_t current = a;
    uint64_t rem;

    groups.reserve(std::ceil(static_cast<double>(current.data.size()) * 64 / std::log2(10) / 19));
    while (current.data.size() > 1 || (!current.data.empty() && current.data[0] != 0ull)) { // while current != 0
        constexpr uint64_t divisor = 10000000000000000000ull;
        big_int_impl::div_abs_inplace(current, divisor, rem);
        groups.push_back(rem);
    }

    constexpr int digits_per_group = 19;
    if (digits > 0) {
        const auto& last_group = std::to_string(groups.back());
        const int total_digits = (groups.size() - 1) * digits_per_group + last_group.length();
        if (total_digits > digits) {
            throw std::runtime_error("small_to_string: result.size()[" + std::to_string(total_digits) +
                "] > digits[" + std::to_string(digits) + "]");
        }
        buffer.append(digits - total_digits, '0');
        buffer.append(last_group);
    } else {
        buffer.append(std::to_string(groups.back()));
    }

    for (int i = groups.size() - 2; i >= 0; i--) {
        auto sgroup = std::to_string(groups[i]);
        buffer.append(digits_per_group - sgroup.length(), '0');
        buffer.append(sgroup);
    }
}

const bint_t& bint_t::compute_power10_with_cache(const int n) {
    if (n < power10_conversion_cache.size())
        return power10_conversion_cache[n];

    const int old_size = power10_conversion_cache.size();
    power10_conversion_cache.resize(n + 1);
    for (int i = old_size; i <= n; i++) {
        bint_t next_value = power10_conversion_cache[i - 1];
        next_value *= next_value;
        power10_conversion_cache[i] = next_value;
    }
    return power10_conversion_cache[n];
}

void bint_t::to_string(const bint_t& a, std::string& buffer, const int digits) { // NOLINT(*-no-recursion)
    if (a.data.size() < TO_STRING_THRESHOLD)
        return small_to_string(a, buffer, digits);

    // Compute n, such as 10^(2^(n + 1)) is approximately a
    const int n = static_cast<int>(std::round(
        std::log2(static_cast<double>(big_int_impl::count_bits(a)) / std::log2(10)) - 1
    ));
    const int expected_digits = 1 << n;
    const bint_t& big10 = compute_power10_with_cache(n);

    bint_t high = a;
    bint_t low;
    big_int_impl::normalize(high);
    big_int::divide_abs(high, big10, low);
    to_string(high, buffer, digits - expected_digits);
    to_string(low, buffer, expected_digits);
}

std::string bint_t::to_string() const {
    std::string result;
    const int expected_length = this->data.size() * 64 / std::log2(10.0) + 1;
    result.reserve(expected_length + (this->sign ? 1 : 0));

    if (this->sign) result.push_back('-');
    to_string(*this, result, 0);
    return result;
}


