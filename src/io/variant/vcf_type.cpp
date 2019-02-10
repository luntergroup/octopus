// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "vcf_type.hpp"

#include <unordered_map>
#include <functional>
#include <stdexcept>

namespace octopus {

VcfType& VcfType::operator+=(const VcfType& rhs)
{
    boost::apply_visitor(detail::inplace_add(), *this, rhs);
    return *this;
}

VcfType& VcfType::operator-=(const VcfType& rhs)
{
    boost::apply_visitor(detail::inplace_subtract(), *this, rhs);
    return *this;
}

VcfType& VcfType::operator*=(const VcfType& rhs)
{
    boost::apply_visitor(detail::inplace_multiply(), *this, rhs);
    return *this;
}

VcfType& VcfType::operator/=(const VcfType& rhs)
{
    boost::apply_visitor(detail::inplace_divide(), *this, rhs);
    return *this;
}

std::string VcfType::type_name() const
{
    return boost::apply_visitor(detail::reflect(), *this);
}

VcfType operator+(VcfType lhs, const VcfType& rhs)
{
    lhs += rhs;
    return lhs;
}

VcfType operator-(VcfType lhs, const VcfType& rhs)
{
    lhs -= rhs;
    return lhs;
}

VcfType operator*(VcfType lhs, const VcfType& rhs)
{
    lhs *= rhs;
    return lhs;
}

VcfType operator/(VcfType lhs, const VcfType& rhs)
{
    lhs /= rhs;
    return lhs;
}

bool operator==(const VcfType& lhs, const VcfType& rhs)
{
    return boost::apply_visitor(detail::is_equal(), lhs, rhs);
}

bool operator!=(const VcfType& lhs, const VcfType& rhs)
{
    return !operator==(lhs, rhs);
}

bool operator<(const VcfType& lhs, const VcfType& rhs)
{
    return boost::apply_visitor(detail::is_less_than(), lhs, rhs);
}

bool operator>(const VcfType& lhs, const VcfType& rhs)
{
    return operator<(rhs, lhs);
}

bool operator<=(const VcfType& lhs, const VcfType& rhs)
{
    return !operator>(lhs, rhs);
}

bool operator>=(const VcfType& lhs, const VcfType& rhs)
{
    return !operator<(lhs, rhs);
}

class UnknownVcfType : std::runtime_error
{
public:
    UnknownVcfType(std::string type)
    : runtime_error {"Invalid VcfType"}
    , type_ {std::move(type)}
    {}
    
    virtual ~UnknownVcfType() noexcept = default;
    
    virtual const char* what() const noexcept
    {
        return (std::string{runtime_error::what()} + ": " + type_).c_str();
    }
    
private:
    std::string type_;
};

class BadVcfType : std::invalid_argument
{
public:
    BadVcfType(std::string type, std::string value)
    : invalid_argument {"invalid VcfType"}
    , type_ {std::move(type)}
    , value_ {std::move(value)}
    {}
    
    virtual ~BadVcfType() noexcept = default;
    
    virtual const char* what() const noexcept
    {
        return (std::string{invalid_argument::what()} + ": " + value_ + " to VcfType " + type_).c_str();
    }
    
private:
    std::string type_, value_;
};

template <typename T>
VcfType make_vcf_type(T&& val)
{
    return VcfType(std::forward<T>(val));
}

VcfType make_vcf_type(const std::string& type, const std::string& value)
{
    static const std::unordered_map<std::string, std::function<VcfType(std::string)>> typeMap {
        {"String",    [] (const auto& value) { return make_vcf_type(value); }},
        {"Integer",   [] (const auto& value) { return make_vcf_type(std::stoi(value)); }},
        {"Float",     [] (const auto& value) { return make_vcf_type(std::stod(value)); }},
        {"Character", [] (const auto& value) { return make_vcf_type(value.front()); }},
        {"Flag",      [] (const auto& value) { return make_vcf_type(value == "1"); }}
    };
    
    if (typeMap.count(type) == 0) throw UnknownVcfType {type};
    
    try {
        return typeMap.at(type)(value);
    } catch (std::invalid_argument& e) {
        throw BadVcfType {type, value};
    } catch (...) {
        throw;
    }
}
    
} // namespace octopus
