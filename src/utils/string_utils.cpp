// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "string_utils.hpp"

#include <boost/algorithm/string/join.hpp>

namespace octopus { namespace utils {

std::vector<std::string> split(const std::string& str, const char delim) {
    std::vector<std::string> elems;
    elems.reserve(std::count(std::cbegin(str), std::cend(str), delim) + 1);
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::string join(const std::vector<std::string>& strings, const std::string delim)
{
    return boost::algorithm::join(strings, delim);
}

std::string join(const std::vector<std::string>& strings, const char delim)
{
    const std::array<char, 2> Delim {delim, '\0'};
    return join(strings, Delim.data());
}

bool is_prefix(const std::string& lhs, const std::string& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

bool is_suffix(const std::string& lhs, const std::string& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::next(std::cbegin(rhs)));
}

std::size_t length(const char* str)
{
    return std::strlen(str);
}

std::size_t length(const std::string& str)
{
    return str.length();
}

bool find(const std::string& lhs, const std::string& rhs)
{
    return lhs.find(rhs) != std::string::npos;
}

std::string& capitalise_front(std::string& str) noexcept
{
    if (!str.empty()) str.front() = std::toupper(str.front());
        return str;
}

std::string capitalise_front(const std::string& str)
{
    auto result = str;
    return capitalise_front(result);
}

bool is_vowel(const char c)
{
    static constexpr std::array<char, 5> vowels {'a', 'e', 'i', 'o', 'u'};
    return std::find(std::cbegin(vowels), std::cend(vowels), std::tolower(c)) != std::cend(vowels);
}

bool begins_with_vowel(const std::string& str)
{
    return !str.empty() && is_vowel(str.front());
}

} // namespace utils
} // namespace octopus
