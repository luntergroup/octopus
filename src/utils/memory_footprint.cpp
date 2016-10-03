#include "memory_footprint.hpp"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <cctype>

namespace octopus {
    
MemoryFootprint::MemoryFootprint(std::size_t num_bytes) noexcept
: num_bytes_ {num_bytes}
{}
    
std::size_t MemoryFootprint::num_bytes() const noexcept
{
    return num_bytes_;
}

std::ostream& operator<<(std::ostream& os, MemoryFootprint footprint)
{
    os << footprint.num_bytes();
    return os;
}

std::istream& operator>>(std::istream& is, MemoryFootprint& result)
{
    if (is.good()) {
        std::string input;
        std::getline(is, input, ' ');
        auto footprint = parse_footprint(input);
        if (footprint) result = *footprint;
    }
    return is;
}

enum class MemoryUnit { kB, KiB, MB, MiB, GB, GiB, TB, TiB, PB, PiB, EB, EiB, ZB, ZiB, YB, YiB, };

boost::optional<MemoryUnit> parse_units(const std::string& str)
{
    static const std::unordered_map<std::string, MemoryUnit> units {
        {"kB", MemoryUnit::kB},
        {"MB", MemoryUnit::MB},
        {"GB", MemoryUnit::GB},
        {"TB", MemoryUnit::TB},
        {"PB", MemoryUnit::PB},
        {"EB", MemoryUnit::EB},
        {"ZB", MemoryUnit::ZB},
        {"YB", MemoryUnit::YB},
        {"KiB", MemoryUnit::KiB},
        {"MiB", MemoryUnit::MiB},
        {"GiB", MemoryUnit::GiB},
        {"TiB", MemoryUnit::TiB},
        {"PiB", MemoryUnit::PiB},
        {"EiB", MemoryUnit::EiB},
        {"ZiB", MemoryUnit::ZiB},
        {"YiB", MemoryUnit::YiB}
    };
    const auto iter = units.find(str);
    if (iter != std::cend(units)) {
        return iter->second;
    } else {
        return boost::none;
    }
}

constexpr std::size_t ipow(std::size_t base, unsigned exp) noexcept
{
    auto result = 1;
    while (exp) {
        if (exp & 1) result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

struct MemoryUnitHash
{
    std::size_t operator()(const MemoryUnit& mu) const noexcept
    {
        return static_cast<std::size_t>(mu);
    }
};

std::size_t get_multiplier(const MemoryUnit units)
{
    static const std::unordered_map<MemoryUnit, std::size_t, MemoryUnitHash> multiplier {
        {MemoryUnit::kB, ipow(1000, 1)},
        {MemoryUnit::MB, ipow(1000, 2)},
        {MemoryUnit::GB, ipow(1000, 3)},
        {MemoryUnit::TB, ipow(1000, 4)},
        {MemoryUnit::PB, ipow(1000, 5)},
        {MemoryUnit::EB, ipow(1000, 6)},
        {MemoryUnit::ZB, ipow(1000, 7)},
        {MemoryUnit::YB, ipow(1000, 8)},
        {MemoryUnit::KiB, ipow(1024, 1)},
        {MemoryUnit::MiB, ipow(1024, 2)},
        {MemoryUnit::GiB, ipow(1024, 3)},
        {MemoryUnit::TiB, ipow(1024, 4)},
        {MemoryUnit::PiB, ipow(1024, 5)},
        {MemoryUnit::EiB, ipow(1024, 6)},
        {MemoryUnit::ZiB, ipow(1024, 7)},
        {MemoryUnit::YiB, ipow(1024, 8)}
    };
    return multiplier.at(units);
}

boost::optional<MemoryFootprint> parse_footprint(std::string str)
{
    const auto iter = std::find_if_not(std::cbegin(str), std::cend(str),
                                       [] (char c) { return std::isdigit(c); });
    if (iter == std::cbegin(str)) return boost::none;
    const auto unit_begin = std::find_if_not(iter, std::cend(str), [] (char c) { return c == ' '; });
    const std::string unit_part {unit_begin, std::cend(str)};
    std::size_t multiplier {1};
    if (!unit_part.empty()) {
        str.erase(iter, std::cend(str));
        const auto units = parse_units(unit_part);
        if (!units) return boost::none;
        multiplier = get_multiplier(*units);
    }
    return MemoryFootprint {multiplier * std::stoll(str)};
}

} // namespace octopus
