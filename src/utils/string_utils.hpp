// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef string_utils_hpp
#define string_utils_hpp

#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <cstddef>
#include <cctype>
#include <array>
#include <sstream>
#include <iomanip>
#include <locale>

#include "maths.hpp"

namespace octopus { namespace utils {

std::vector<std::string> split(const std::string& str, const char delim);
std::vector<std::string> split(const std::string& str, const std::string delims);

std::string join(const std::vector<std::string>& strings, const std::string delim = "");
std::string join(const std::vector<std::string>& strings, const char delim);

bool is_prefix(const std::string& prefix, const std::string& text) noexcept;
bool is_suffix(const std::string& suffix, const std::string& text) noexcept;

std::size_t length(const char* str);
std::size_t length(const std::string& str);

bool find(const std::string& lhs, const std::string& rhs);

std::string& capitalise(std::string& str) noexcept;
std::string capitalise(const std::string& str);
std::string& capitalise_front(std::string& str) noexcept;
std::string capitalise_front(const std::string& str);
std::string& to_lower(std::string& str) noexcept;
std::string to_lower(const std::string& str);
std::string& strip_leading_zeroes(std::string& str);
std::string strip_leading_zeroes(const std::string& str);

bool is_vowel(const char c);
bool begins_with_vowel(const std::string& str);

enum class PrecisionRule { dp, sf };

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
std::string to_string(const T val, const unsigned precision = 2, const PrecisionRule rule = PrecisionRule::dp)
{
    if (rule == PrecisionRule::dp) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(precision) << val;
        return out.str();
    } else {
        return to_string(val, maths::count_leading_zeros(val) + precision, PrecisionRule::dp);
    }
}

template <typename T>
std::vector<std::string> to_strings(const std::vector<T>& values)
{
    std::vector<std::string> result {};
    result.reserve(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::back_inserter(result),
                   [] (auto value) { return std::to_string(value); });
    return result;
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
std::vector<std::string> to_strings(const std::vector<T>& values, const unsigned precision = 2)
{
    std::vector<std::string> result {};
    result.reserve(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::back_inserter(result),
                   [precision] (auto value) { return to_string(value, precision); });
    return result;
}

namespace detail {
    class CommaNumPunct : public std::numpunct<char>
    {
    protected:
        virtual char do_thousands_sep() const noexcept override
        {
            return ',';
        }
        
        virtual std::string do_grouping() const override
        {
            return "\03";
        }
    };
}

template <typename T>
std::string format_with_commas(const T value)
{
    static_assert(std::is_arithmetic<T>::value, "T must be numeric");
    // std::locale is responsable for calling ~detail::CommaNumPunct
    std::locale comma_locale {std::locale(), new detail::CommaNumPunct {}};
    std::stringstream ss;
    ss.imbue(comma_locale);
    ss << std::fixed << value;
    return ss.str();
}

} // namespace utils
} // namespace octopus

#endif
