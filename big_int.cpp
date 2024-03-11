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

std::string bint_t::to_string() const {
    std::string result;
    bint_t current = *this;
    uint64_t rem;
    while (current.data.size() > 1 || current.data[0] != 0ull) { // while current != 0
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
    if (this->sign) result.push_back('-');
    std::ranges::reverse(result);
    return result;
}


