//
//  vcf_type.cpp
//  Octopus
//
//  Created by Daniel Cooke on 06/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_type.h"

#include <unordered_map>
#include <functional>
#include <stdexcept>

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

std::string VcfType::name() const
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

class UnknownVcfType : std::runtime_error {
public:
    UnknownVcfType(std::string type)
    :
    runtime_error {"Invalid VcfType"},
    type_ {std::move(type)}
    {}
    
    const char* what() const noexcept
    {
        return (std::string{runtime_error::what()} + ": " + type_).c_str();
    }
    
private:
    std::string type_;
};

class BadVcfType : std::invalid_argument {
public:
    BadVcfType(std::string type, std::string value)
    :
    invalid_argument {"Incompatible VcfType"},
    type_ {std::move(type)},
    value_ {std::move(value)}
    {}
    
    const char* what() const noexcept
    {
        return (std::string{invalid_argument::what()} + ": " + value_ + " to VcfType " + type_).c_str();
    }
    
private:
    std::string type_, value_;
};

VcfType make_vcf_type(const std::string& type, const std::string& value)
{
    static const std::unordered_map<std::string, std::function<VcfType(std::string)>> type_map {
        {"String", [] (const auto& value) { return detail::make_vcf_type(value); }},
        {"Integer",   [] (const auto& value) { return detail::make_vcf_type(std::stoi(value)); }},
        {"Float",     [] (const auto& value) { return detail::make_vcf_type(std::stod(value)); }},
        {"Character", [] (const auto& value) { return detail::make_vcf_type(value.front()); }},
        {"Flag",      [] (const auto& value) { return detail::make_vcf_type(value == "1"); }}
    };
    
    if (type_map.count(type) == 0) throw UnknownVcfType {type};
    
    try {
        return type_map.at(type)(value);
    } catch (std::invalid_argument& e) {
        throw BadVcfType {type, value};
    } catch (...) { // e.g. std::out_of_range
        throw;
    }
    
}
