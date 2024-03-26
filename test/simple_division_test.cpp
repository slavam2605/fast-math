#include <gtest/gtest.h>
#include <big_int/big_int.h>
#include <big_int/big_int_ops.h>

std::pair<bint_t, bint_t> get_test_div_values(const bool first_sign, const bool second_sign) {
    bint_t a(3);
    big_int::fast_pow_inplace(a, 100);
    a.sign = first_sign;
    return std::make_pair(a, bint_t(1000000000000000000ll * (second_sign ? -1 : 1)));
}

TEST(SimpleDivisionTest, TestPosPos) {
    auto [a, b] = get_test_div_values(false, false);
    const auto mod = a % b;
    const auto div = a / b;
    EXPECT_EQ(mod, bint_t(621272702107522001ll));
    EXPECT_STREQ(div.to_string().c_str(), "515377520732011331036461129765");
}

TEST(SimpleDivisionTest, TestPosNeg) {
    auto [a, b] = get_test_div_values(false, true);
    const auto mod = a % b;
    const auto div = a / b;
    EXPECT_EQ(mod, bint_t(-378727297892477999ll));
    EXPECT_STREQ(div.to_string().c_str(), "-515377520732011331036461129766");
}

TEST(SimpleDivisionTest, TestNegPos) {
    auto [a, b] = get_test_div_values(true, false);
    const auto mod = a % b;
    const auto div = a / b;
    EXPECT_EQ(mod, bint_t(378727297892477999ll));
    EXPECT_STREQ(div.to_string().c_str(), "-515377520732011331036461129766");
}

TEST(SimpleDivisionTest, TestNegNeg) {
    auto [a, b] = get_test_div_values(true, true);
    const auto mod = a % b;
    const auto div = a / b;
    EXPECT_EQ(mod, bint_t(-621272702107522001ll));
    EXPECT_STREQ(div.to_string().c_str(), "515377520732011331036461129765");
}

TEST(SimpleDivisionTest, TestPosPosAssign) {
    auto [a, b] = get_test_div_values(false, false);
    auto c = a;
    a %= b;
    c /= b;
    EXPECT_EQ(a, bint_t(621272702107522001ll));
    EXPECT_STREQ(c.to_string().c_str(), "515377520732011331036461129765");
}

TEST(SimpleDivisionTest, TestPosNegAssign) {
    auto [a, b] = get_test_div_values(false, true);
    auto c = a;
    a %= b;
    c /= b;
    EXPECT_EQ(a, bint_t(-378727297892477999ll));
    EXPECT_STREQ(c.to_string().c_str(), "-515377520732011331036461129766");
}

TEST(SimpleDivisionTest, TestNegPosAssign) {
    auto [a, b] = get_test_div_values(true, false);
    auto c = a;
    a %= b;
    c /= b;
    EXPECT_EQ(a, bint_t(378727297892477999ll));
    EXPECT_STREQ(c.to_string().c_str(), "-515377520732011331036461129766");
}

TEST(SimpleDivisionTest, TestNegNegAssign) {
    auto [a, b] = get_test_div_values(true, true);
    auto c = a;
    a %= b;
    c /= b;
    EXPECT_EQ(a, bint_t(-621272702107522001ll));
    EXPECT_STREQ(c.to_string().c_str(), "515377520732011331036461129765");
}