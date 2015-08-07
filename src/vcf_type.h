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

using VcfType = boost::variant<int, double, char, std::string, bool>;

namespace detail
{
    template <typename T, typename U>
    using enable_if_both_arithmetic = std::enable_if_t<std::is_arithmetic<T>::value && std::is_arithmetic<U>::value, VcfType>;
    
    template <typename T, typename U>
    using enable_if_not_both_arithmetic = std::enable_if_t<!(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value), VcfType>;
    
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
}

VcfType operator+(const VcfType& lhs, const VcfType& rhs);
VcfType operator-(const VcfType& lhs, const VcfType& rhs);
VcfType operator*(const VcfType& lhs, const VcfType& rhs);
VcfType operator/(const VcfType& lhs, const VcfType& rhs);

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
