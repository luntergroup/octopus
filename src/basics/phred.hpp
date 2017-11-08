// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef phred_hpp
#define phred_hpp

#include <type_traits>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <iostream>
#include <functional>

#include "concepts/comparable.hpp"

namespace octopus {

template <typename Q = double, typename = std::enable_if_t<std::is_arithmetic<Q>::value>>
class Phred : public Comparable<Phred<Q>>
{
public:
    struct Probability
    {
        Q value;
        operator Q() const noexcept { return value; }
    };
    
    Phred() = default;
    
    explicit Phred(const Q score) : score_ {score} // to convert -0.0 to +0.0
    {
        if (score_ < Q {0}) {
            throw std::domain_error {"Phred: negative score " + std::to_string(score)};
        }
    }
    
    explicit Phred(const Probability error)
    {
        if (error < Q {0}) {
            throw std::domain_error {"Phred: negative error probability " + std::to_string(error)};
        }
        static const Q minProbability = std::nextafter(Q {0}, Q {1});
        score_ = std::abs(Q {-10} * std::log10(std::max(std::min(error.value, Q {1}), minProbability)));
    }
    
    Phred(const Phred&)            = default;
    Phred& operator=(const Phred&) = default;
    Phred(Phred&&)                 = default;
    Phred& operator=(Phred&&)      = default;
    
    ~Phred() = default;
    
    Q score() const noexcept
    {
        return score_;
    }
    
    explicit operator Q() const noexcept
    {
        return score_;
    }
    
    Probability probability_true() const noexcept
    {
        static constexpr Q one {1};
        return Probability {one - probability_false()};
    }
    
    Probability probability_false() const noexcept
    {
        static constexpr Q ten {10};
        return Probability {std::pow(ten, -score_ / ten)};
    }
private:
    Q score_;
};

template <typename Q>
auto probability_to_phred(const Q p)
{
    return Phred<Q> {typename Phred<Q>::Probability {p}};
}

template <typename Q>
bool operator==(const Phred<Q>& lhs, const Phred<Q>& rhs) noexcept
{
    return lhs.score() == rhs.score();
}

template <typename Q>
bool operator<(const Phred<Q>& lhs, const Phred<Q>& rhs) noexcept
{
    return lhs.score() < rhs.score();
}

template <typename Q>
std::ostream& operator<<(std::ostream& os, const Phred<Q>& p)
{
    os << p.score();
    return os;
}

template <typename Q>
std::istream& operator>>(std::istream& is, Phred<Q>& p)
{
    Q value;
    is >> value;
    p = Phred<Q> {value};
    return is;
}

template <typename Q>
std::string to_string(const Phred<Q>& phred)
{
    return std::to_string(phred.score());
}

} // namespace octopus

namespace std {
    template <typename T> struct hash<octopus::Phred<T>>
    {
        size_t operator()(const octopus::Phred<T>& phred) const noexcept
        {
            return hash<T>(phred.score());
        }
    };
}

#endif
