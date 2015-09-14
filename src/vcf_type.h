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
#include <boost/type_index.hpp>

namespace detail
{
    using VcfTypeBase = boost::variant<int, double, char, std::string, bool>;
    
    template <typename T>
    std::string type_name()
    {
        return boost::typeindex::type_id<T>().pretty_name();
    }
    
    template <typename T>
    std::string type_name(const T& v)
    {
        return boost::typeindex::type_id<T>().pretty_name();
    }
    
    template <typename T>
    class convert : public boost::static_visitor<T>
    {
    public:
        template <typename U>
        T operator()(const U& v) const
        {
            return static_cast<T>(v);
        }
        
        T operator()(const std::string& v) const
        {
            return do_convert(v, T {});
        }
        
    private:
        int do_convert(const std::string& v, int) const
        {
            return std::stoi(v);
        }
        
        double do_convert(const std::string& v, double) const
        {
            return std::stod(v);
        }
        
        template <typename U>
        U do_convert(const std::string& v, U) const
        {
            throw std::runtime_error {"Cannot convert std::string to type " + type_name<U>()};
        }
    };
    
    template <>
    class convert<std::string> : public boost::static_visitor<std::string>
    {
    public:
        template <typename U>
        std::string operator()(const U& v) const
        {
            return std::to_string(v);
        }
        
        std::string operator()(const std::string& v) const
        {
            return v;
        }
    };
    
    template <typename T, typename U, typename R>
    using enable_if_both_arithmetic = std::enable_if_t<std::is_arithmetic<T>::value && std::is_arithmetic<U>::value, R>;
    
    template <typename T, typename U, typename R>
    using enable_if_not_both_arithmetic = std::enable_if_t<!(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value), R>;
    
    class VcfTypeError : std::runtime_error {
    public:
        VcfTypeError(std::string op, std::string type1, std::string type2)
        :
        runtime_error {"invalid operation on types"},
        op_ {std::move(op)},
        type1_ {std::move(type1)},
        type2_ {std::move(type2)}
        {}
        
        const char* what() const noexcept
        {
            return (std::string {runtime_error::what()} + ": operation=" + op_ + " lhs= " + type1_ + " rhs=" + type2_).c_str();
        }
        
    private:
        std::string op_, type1_, type2_;
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
            throw VcfTypeError {"==", type_name<T>(), type_name<U>()};
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
            throw VcfTypeError {"<", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_add : public boost::static_visitor<void>
    {
    public:
        template <typename T>
        void operator()(T& lhs, const T& rhs) const // this case is really just for std::string
        {
            lhs += rhs;
        }
        
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs += rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"+=", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_subtract : public boost::static_visitor<void>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs -= rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"-=", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_multiply : public boost::static_visitor<void>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs *= rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"*=", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_divide : public boost::static_visitor<void>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs /= rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"/=", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_mod : public boost::static_visitor<void>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs %= rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"%=", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_bit_and : public boost::static_visitor<void>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs &= rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"&=", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_bit_or : public boost::static_visitor<void>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs |= rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"|=", type_name<T>(), type_name<U>()};
        }
    };
    
    class inplace_bit_xor : public boost::static_visitor<void>
    {
    public:
        template <typename T, typename U>
        enable_if_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            lhs ^= rhs;
        }
        
        template <typename T, typename U>
        enable_if_not_both_arithmetic<T, U, void> operator()(T& lhs, const U& rhs) const
        {
            throw VcfTypeError {"^=", type_name<T>(), type_name<U>()};
        }
    };
    
    class reflect : public boost::static_visitor<std::string>
    {
    public:
        template <typename T>
        std::string operator()(const T& v) const
        {
            return type_name(v);
        }
    };
}

// subclass from boost::variant so we can define our own comparison operators
class VcfType : public detail::VcfTypeBase
{
public:
    // need first constructor to avoid const char* converting to bool
    explicit VcfType(const char* str) : detail::VcfTypeBase {std::string {str}} {}
    template <typename T, typename = std::enable_if_t<!std::is_base_of<VcfType, std::decay_t<T>>::value>>
    explicit VcfType(T&& v) : detail::VcfTypeBase {std::forward<T>(v)} {}
    
    VcfType(const VcfType&)            = default;
    VcfType& operator=(const VcfType&) = default;
    VcfType(VcfType&&)                 = default;
    VcfType& operator=(VcfType&&)      = default;
    
    template <typename T>
    operator T() const
    {
        return boost::apply_visitor(detail::convert<T>(), *this);
    }
    
    VcfType& operator+=(const VcfType& rhs);
    VcfType& operator-=(const VcfType& rhs);
    VcfType& operator*=(const VcfType& rhs);
    VcfType& operator/=(const VcfType& rhs);
    
    std::string type_name() const;
};

VcfType operator+(VcfType lhs, const VcfType& rhs);
VcfType operator-(VcfType lhs, const VcfType& rhs);
VcfType operator*(VcfType lhs, const VcfType& rhs);
VcfType operator/(VcfType lhs, const VcfType& rhs);
bool operator==(const VcfType& lhs, const VcfType& rhs);
bool operator!=(const VcfType& lhs, const VcfType& rhs);
bool operator<(const VcfType& lhs, const VcfType& rhs);
bool operator>(const VcfType& lhs, const VcfType& rhs);
bool operator<=(const VcfType& lhs, const VcfType& rhs);
bool operator>=(const VcfType& lhs, const VcfType& rhs);

VcfType make_vcf_type(const std::string& type, const std::string& value);

#endif
