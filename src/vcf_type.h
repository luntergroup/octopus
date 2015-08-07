//
//  vcf_type.h
//  Octopus
//
//  Created by Daniel Cooke on 06/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_vcf_type_h
#define Octopus_vcf_type_h

#include <string>
#include <type_traits>
#include <stdexcept>
#include <boost/variant.hpp>

#include <iostream> // TEST

namespace detail
{
    using VcfTypeBase = boost::variant<int, double, char, std::string, bool>;
}

// subclass from boost::variant so we can define our own comparison operators
class VcfType : public detail::VcfTypeBase
{
public:
    template <typename T> VcfType(T&& v) : detail::VcfTypeBase {std::forward<T>(v)} {}
};

class VcfTypeBinaryOpError : std::runtime_error {
public:
    VcfTypeBinaryOpError(std::string op, std::string type1, std::string type2)
    :
    runtime_error {"Cannot apply binary operation to given types"},
    op_ {std::move(op)},
    type1_ {std::move(type1)},
    type2_ {std::move(type2)}
    {}
    
    const char* what() const noexcept
    {
        return (std::string{runtime_error::what()} + ": operation=" + op_ + " type1= " + type1_ + " type2=" + type2_).c_str();
    }
    
private:
    std::string op_, type1_, type2_;
};

namespace detail
{
    template <typename T, typename U, typename R = VcfType>
    using enable_if_both_arithmetic = std::enable_if_t<std::is_arithmetic<T>::value && std::is_arithmetic<U>::value, R>;
    
    template <typename T, typename U, typename R = VcfType>
    using enable_if_not_both_arithmetic = std::enable_if_t<!(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value), R>;
    
    class add : public boost::static_visitor<VcfType>
    {
    public:
        std::string operator()(const std::string& lhs, const std::string& rhs) const
        {
            return lhs + rhs;
        }
        
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            return lhs + rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            throw std::logic_error {"cannot add these types"};
        }
    };
    
    class subtract : public boost::static_visitor<VcfType>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            return lhs - rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            throw std::logic_error {"cannot subtract these types"};
        }
    };
    
    class multiply : public boost::static_visitor<VcfType>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            return lhs * rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            throw std::logic_error {"cannot multiply these types"};
        }
    };
    
    class divide : public boost::static_visitor<VcfType>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            return lhs / rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U> operator()(const T& lhs, const U& rhs) const
        {
            throw std::logic_error {"cannot divide these types"};
        }
    };
    
    class is_equal : public boost::static_visitor<bool>
    {
    public:
        template <typename T>
        bool operator()(const T& lhs, const T& rhs) const
        {
            return lhs == rhs;
        }
        
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, bool> operator()(const T& lhs, const U& rhs) const
        {
            return lhs == rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, bool> operator()(const T& lhs, const U& rhs) const
        {
            return false;
        }
    };
    
    class is_less_than : public boost::static_visitor<bool>
    {
    public:
        template <typename T>
        bool operator()(const T& lhs, const T& rhs) const
        {
            return lhs < rhs;
        }
        
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, bool> operator()(const T& lhs, const U& rhs) const
        {
            return lhs < rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, bool> operator()(const T& lhs, const U& rhs) const
        {
            return false;
        }
    };
}

VcfType operator+(const VcfType& lhs, const VcfType& rhs);
VcfType operator-(const VcfType& lhs, const VcfType& rhs);
VcfType operator*(const VcfType& lhs, const VcfType& rhs);
VcfType operator/(const VcfType& lhs, const VcfType& rhs);
bool operator==(const VcfType& lhs, const VcfType& rhs);
bool operator!=(const VcfType& lhs, const VcfType& rhs);
bool operator<(const VcfType& lhs, const VcfType& rhs);
bool operator>(const VcfType& lhs, const VcfType& rhs);
bool operator<=(const VcfType& lhs, const VcfType& rhs);
bool operator>=(const VcfType& lhs, const VcfType& rhs);

namespace detail
{
    template <typename T>
    VcfType make_vcf_type(T&& val)
    {
        return VcfType(std::forward<T>(val));
    }
}

VcfType vcf_type_factory(const std::string& type, const std::string& value);

#endif
