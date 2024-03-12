//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#include "big_int.h"
#include "big_int_ops.h"
#include "karatsuba.h"

bint_t::bint_t() : sign(false) {}

bint_t::bint_t(const uint64_t value) : sign(false) {
    data.push_back(value);
}

bint_t::bint_t(const int64_t value) : sign(value < 0) {
    data.push_back(abs(value));
}

void bint_t::normalize() {
    while (data.size() > 1 && data.back() == 0) {
        data.pop_back();
    }
    if (data.empty() || data.size() == 1 && data[0] == 0) {
        sign = false;
    }
}

std::strong_ordering bint_t::operator<=>(const bint_t& other) const {
    return compare(*this, other);
}

bool bint_t::operator==(const bint_t& other) const {
    return (*this <=> other) == 0;
}

bint_t bint_t::operator+(const bint_t& other) const {
    return add(*this, other);
}

bint_t bint_t::operator-(const bint_t& other) const {
    return sub(*this, other);
}

bint_t bint_t::operator*(const bint_t& other) const {
    return karatsuba(*this, other);
}

void small_to_string(const bint_t& a, std::string& buffer, int digits) {
    std::string result;
    bint_t current = a;
    uint64_t rem;
    while (current.data.size() > 1 || (!current.data.empty() && current.data[0] != 0ull)) { // while current != 0
        constexpr uint64_t divisor = 10000000000000000000ull;
        div_abs_inplace(current, divisor, rem);
        for (int _ = 0; _ < 19; _++) {
            result.push_back("0123456789"[rem % 10]);
            rem /= 10;
        }
    }
    while (result.size() > 1 && result.back() == '0') {
        result.pop_back();
    }
    if (digits > 0) {
        if (result.size() > digits) throw std::runtime_error("small_to_string: result.size() > digits");
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
    const bint_t big10 = fast_pow(bint_t(10ll), half_length10);

    bint_t high = a;
    bint_t low;
    div_abs_inplace(high, big10, low);
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


