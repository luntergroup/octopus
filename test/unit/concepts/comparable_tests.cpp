// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "concepts/equitable.hpp"
#include "concepts/comparable.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(concepts)
BOOST_AUTO_TEST_SUITE(comparable)

namespace {

struct Foo : public Equitable<Foo>
{
    explicit Foo(int x) : x {x} {}
    int x;
};

bool operator==(const Foo& lhs, const Foo& rhs) noexcept { return lhs.x == rhs.x; }

std::ostream& operator<<(std::ostream& os, const Foo& f)
{
    os << f.x;
    return os;
}

struct Bar : public Comparable<Bar>
{
    explicit Bar(int x) : x {x} {}
    int x;
};

bool operator==(const Bar& lhs, const Bar& rhs) noexcept { return lhs.x == rhs.x; }
bool operator<(const Bar& lhs, const Bar& rhs) noexcept { return lhs.x < rhs.x; }

std::ostream& operator<<(std::ostream& os, const Bar& b)
{
    os << b.x;
    return os;
}

} // namespace

BOOST_AUTO_TEST_CASE(equitable_types_automatically_get_the_not_equals_operator)
{
    Foo f1 {0}, f2 {1}, f3 {2};
    
    BOOST_CHECK_NE(f1, f2);
    BOOST_CHECK_NE(f1, f3);
    BOOST_CHECK_NE(f2, f3);
}

BOOST_AUTO_TEST_CASE(comparable_types_automatically_get_all_comparison_operators)
{
    Bar b1 {0}, b2 {1}, b3 {2};
    
    BOOST_CHECK_LT(b1, b2);
    BOOST_CHECK_LT(b1, b3);
    BOOST_CHECK_LT(b2, b3);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
    
} // namespace test
} // namespace octopus
