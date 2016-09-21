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

constexpr std::size_t cpow(const std::size_t a, unsigned n) noexcept
{
    auto result = a;
    for (; n > 0; --n) result *= a;
    return result;
}

std::size_t get_multiplier(const MemoryUnit units)
{
    static const std::unordered_map<MemoryUnit, std::size_t> multiplier {
        {MemoryUnit::kB, cpow(1000, 1)},
        {MemoryUnit::MB, cpow(1000, 2)},
        {MemoryUnit::GB, cpow(1000, 3)},
        {MemoryUnit::TB, cpow(1000, 4)},
        {MemoryUnit::PB, cpow(1000, 5)},
        {MemoryUnit::EB, cpow(1000, 6)},
        {MemoryUnit::ZB, cpow(1000, 7)},
        {MemoryUnit::YB, cpow(1000, 8)},
        {MemoryUnit::KiB, cpow(1024, 1)},
        {MemoryUnit::MiB, cpow(1024, 2)},
        {MemoryUnit::GiB, cpow(1024, 3)},
        {MemoryUnit::TiB, cpow(1024, 4)},
        {MemoryUnit::PiB, cpow(1024, 5)},
        {MemoryUnit::EiB, cpow(1024, 6)},
        {MemoryUnit::ZiB, cpow(1024, 7)},
        {MemoryUnit::YiB, cpow(1024, 8)}
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
