// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "memory_footprint.hpp"

#include <array>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "string_utils.hpp"

namespace octopus {

bool operator==(const MemoryFootprint& lhs, const MemoryFootprint& rhs) noexcept
{
    return lhs.bytes() == rhs.bytes();
}

bool operator<(const MemoryFootprint& lhs, const MemoryFootprint& rhs) noexcept
{
    return lhs.bytes() < rhs.bytes();
}

namespace {

enum class MemoryUnit { B, kB, KiB, MB, MiB, GB, GiB, TB, TiB, PB, PiB, EB, EiB, ZB, ZiB, YB, YiB, };

boost::optional<MemoryUnit> parse_units(std::string& units_str)
{
    static const std::unordered_map<std::string, MemoryUnit> units
    {
        {"B",   MemoryUnit::B},
        {"K",   MemoryUnit::kB},
        {"KB",  MemoryUnit::kB},
        {"M",   MemoryUnit::MB},
        {"MB",  MemoryUnit::MB},
        {"G",   MemoryUnit::GB},
        {"GB",  MemoryUnit::GB},
        {"T",   MemoryUnit::TB},
        {"TB",  MemoryUnit::TB},
        {"P",   MemoryUnit::PB},
        {"PB",  MemoryUnit::PB},
        {"E",   MemoryUnit::EB},
        {"EB",  MemoryUnit::EB},
        {"Z",   MemoryUnit::ZB},
        {"ZB",  MemoryUnit::ZB},
        {"Y",   MemoryUnit::YB},
        {"YB",  MemoryUnit::YB},
        {"KIB", MemoryUnit::KiB},
        {"MIB", MemoryUnit::MiB},
        {"GIB", MemoryUnit::GiB},
        {"TIB", MemoryUnit::TiB},
        {"PIB", MemoryUnit::PiB},
        {"EIB", MemoryUnit::EiB},
        {"ZIB", MemoryUnit::ZiB},
        {"YIB", MemoryUnit::YiB}
    };
    utils::capitalise(units_str);
    const auto itr = units.find(units_str);
    if (itr != std::cend(units)) {
        return itr->second;
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

constexpr std::size_t get_multiplier(const MemoryUnit unit) noexcept
{
    switch (unit) {
        case MemoryUnit::B:   return ipow(1000, 0);
        case MemoryUnit::kB:  return ipow(1000, 1);
        case MemoryUnit::MB:  return ipow(1000, 2);
        case MemoryUnit::GB:  return ipow(1000, 3);
        case MemoryUnit::TB:  return ipow(1000, 4);
        case MemoryUnit::PB:  return ipow(1000, 5);
        case MemoryUnit::EB:  return ipow(1000, 6);
        case MemoryUnit::ZB:  return ipow(1000, 7);
        case MemoryUnit::YB:  return ipow(1000, 8);
        case MemoryUnit::KiB: return ipow(1024, 1);
        case MemoryUnit::MiB: return ipow(1024, 2);
        case MemoryUnit::GiB: return ipow(1024, 3);
        case MemoryUnit::TiB: return ipow(1024, 4);
        case MemoryUnit::PiB: return ipow(1024, 5);
        case MemoryUnit::EiB: return ipow(1024, 6);
        case MemoryUnit::ZiB: return ipow(1024, 7);
        case MemoryUnit::YiB: return ipow(1024, 8);
        default: return ipow(1000, 0);
    }
}

std::string to_string(const MemoryUnit unit)
{
    switch (unit) {
        case MemoryUnit::B:   return "B";
        case MemoryUnit::kB:  return "kB";
        case MemoryUnit::MB:  return "MB";
        case MemoryUnit::GB:  return "GB";
        case MemoryUnit::TB:  return "TB";
        case MemoryUnit::PB:  return "PB";
        case MemoryUnit::EB:  return "EB";
        case MemoryUnit::ZB:  return "ZB";
        case MemoryUnit::YB:  return "YB";
        case MemoryUnit::KiB: return "KiB";
        case MemoryUnit::MiB: return "MiB";
        case MemoryUnit::GiB: return "GiB";
        case MemoryUnit::TiB: return "TiB";
        case MemoryUnit::PiB: return "PiB";
        case MemoryUnit::EiB: return "EiB";
        case MemoryUnit::ZiB: return "ZiB";
        case MemoryUnit::YiB: return "YiB";
        default: return "B";
    }
}

auto get_human_format_units(std::size_t bytes) noexcept
{
    using MU = MemoryUnit;
    if (bytes >= 1000) {
        if (bytes % 10 == 0) {
            static constexpr std::array<MU, 8> units { MU::kB, MU::MB, MU::GB, MU::TB, MU::PB, MU::EB, MU::ZB, MU::YB };
            std::size_t unit_idx {0};
            bytes /= 1000;
            while (bytes >= 1000 && bytes % 10 == 0) {
                bytes /= 1000;
                ++unit_idx;
            }
            return units[std::min(unit_idx, units.size() - 1)];
        } else {
            static constexpr std::array<MU, 9> units { MU::B, MU::KiB, MU::MiB, MU::GiB, MU::TiB, MU::PiB, MU::EiB, MU::ZiB, MU::YiB };
            const auto unit_idx = static_cast<std::size_t>(std::floor(std::log(bytes) / std::log(1024)));
            const auto dv = std::lldiv(bytes, std::pow(1024, unit_idx));
            if (dv.rem == 0 || (dv.rem % 2 == 0 && std::log2(dv.rem) > 6)) {
                return units[std::min(unit_idx, units.size() - 1)];
            }
        }
    }
    return MU::B;
}

auto get_human_format_units(const MemoryFootprint& footprint) noexcept
{
    return get_human_format_units(footprint.num_bytes());
}

} // namespace

boost::optional<MemoryFootprint> parse_footprint(std::string footprint_str)
{
    using std::cbegin; using std::cend;
    static const auto is_digit = [] (char c) { return std::isdigit(c); };
    auto last_digit_itr = std::find_if_not(cbegin(footprint_str), cend(footprint_str), is_digit);
    if (last_digit_itr == cend(footprint_str)) {
        MemoryFootprint {static_cast<std::size_t>(std::stoll(footprint_str))};
    }
    bool is_float {false};
    if (*last_digit_itr == '.') {
        last_digit_itr = std::find_if_not(std::next(last_digit_itr), cend(footprint_str), is_digit);
        is_float = true;
    }
    if (last_digit_itr == cbegin(footprint_str)) return boost::none;
    const auto unit_begin_itr = std::find_if_not(last_digit_itr, cend(footprint_str), [] (char c) { return c == ' '; });
    std::string unit_part {unit_begin_itr, cend(footprint_str)};
    std::size_t multiplier {1};
    if (!unit_part.empty()) {
        footprint_str.erase(last_digit_itr, cend(footprint_str));
        const auto units = parse_units(unit_part);
        if (!units) return boost::none;
        multiplier = get_multiplier(*units);
    }
    std::size_t bytes {};
    if (is_float) {
        bytes = multiplier * std::stod(footprint_str);
    } else {
        bytes = multiplier * std::stoll(footprint_str);
    }
    return MemoryFootprint {bytes};
}

std::ostream& operator<<(std::ostream& os, MemoryFootprint footprint)
{
    const auto units = get_human_format_units(footprint);
    const auto multiplier = get_multiplier(units);
    if (footprint.bytes() % multiplier == 0) {
        os << footprint.bytes() / multiplier;
    } else {
        os << static_cast<double>(footprint.bytes()) / multiplier;
    }
    os << to_string(units);
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

} // namespace octopus
