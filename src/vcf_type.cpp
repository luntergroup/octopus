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

VcfType operator+(const VcfType& lhs, const VcfType& rhs)
{
    return boost::apply_visitor(detail::add(), lhs, rhs);
}

VcfType operator-(const VcfType& lhs, const VcfType& rhs)
{
    return boost::apply_visitor(detail::subtract(), lhs, rhs);
}

VcfType operator*(const VcfType& lhs, const VcfType& rhs)
{
    return boost::apply_visitor(detail::multiply(), lhs, rhs);
}

VcfType operator/(const VcfType& lhs, const VcfType& rhs)
{
    return boost::apply_visitor(detail::divide(), lhs, rhs);
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

VcfType vcf_type_factory(const std::string& type, const std::string& value)
{
    static std::unordered_map<std::string, std::function<VcfType(std::string)>> type_map;
    type_map.emplace("String",    [] (const auto& value) { return detail::make_vcf_type(value); });
    type_map.emplace("Integer",   [] (const auto& value) { return detail::make_vcf_type(std::stoi(value)); });
    type_map.emplace("Float",     [] (const auto& value) { return detail::make_vcf_type(std::stod(value)); });
    type_map.emplace("Character", [] (const auto& value) { return detail::make_vcf_type(value.front()); });
    type_map.emplace("Flag",      [] (const auto& value) { return detail::make_vcf_type(value == "1"); });
    
    if (type_map.count(type) == 0) throw UnknownVcfType {type};
    
    try {
        return type_map[type](value);
    } catch (std::invalid_argument& e) {
        throw BadVcfType {type, value};
    } catch (...) { // e.g. std::out_of_range
        throw;
    }
    
}
