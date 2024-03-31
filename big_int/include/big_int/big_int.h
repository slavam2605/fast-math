//
// Created by Vyacheslav.Moklev on 09/03/2024.
//

#ifndef BIG_INT_H
#define BIG_INT_H

#include <vector>
#include <string>

struct bint_t {
    bool sign;
    std::vector<uint64_t> data;

    bint_t();
    bint_t(const bint_t& other) = default;
    bint_t(bint_t&& other) noexcept;
    explicit bint_t(uint64_t value);
    explicit bint_t(int64_t value);
    explicit bint_t(int value);

    bint_t& pow(int n);
    [[nodiscard]] bint_t pow(int n) const;

    bint_t& operator=(const bint_t& other) = default;
    bint_t& operator=(bint_t&& other) = default;
    std::strong_ordering operator<=>(const bint_t& other) const;
    bool operator==(const bint_t& other) const;
    bint_t operator-() const;

    bint_t operator+(const bint_t& other) const;
    bint_t operator-(const bint_t& other) const;
    bint_t operator*(const bint_t& other) const;
    bint_t operator/(const bint_t& other) const;
    bint_t operator%(const bint_t& other) const;
    bint_t operator<<(int64_t n) const;
    bint_t operator>>(int64_t n) const;

    bint_t& operator+=(const bint_t& other);
    bint_t& operator-=(const bint_t& other);
    bint_t& operator*=(const bint_t& other);
    bint_t& operator/=(const bint_t& other);
    bint_t& operator%=(const bint_t& other);
    bint_t& operator<<=(int64_t n);
    bint_t& operator>>=(int64_t n);

    [[nodiscard]] std::string to_string() const;
    [[nodiscard]] std::string to_string_old() const;

private:
    static constexpr int TO_STRING_THRESHOLD = 20;
    static std::vector<bint_t> power10_conversion_cache;

    static const bint_t& compute_power10_with_cache(int n);
    static void to_string(const bint_t& a, std::string& buffer, int digits);
};

#endif //BIG_INT_H
