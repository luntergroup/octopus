// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

// Based of https://stackoverflow.com/a/15166623/2970186

#ifndef beta_distribution_hpp
#define beta_distribution_hpp

#include <random>

namespace std {

template <typename RealType = double>
class beta_distribution
{
public:
    using result_type = RealType;
    
    class param_type
    {
    public:
        using distribution_type = beta_distribution;
        
        explicit param_type(RealType a, RealType) : a_param(a), b_param(b) {}
        
        RealType a() const noexcept { return a_param; }
        RealType b() const noexcept { return b_param; }
        
        bool operator==(const param_type& other) const noexcept
        {
            return (a_param == other.a_param && b_param == other.b_param);
        }
        
        bool operator!=(const param_type& other) const noexcept
        {
            return !(*this == other);
        }
    
    private:
        RealType a_param, b_param;
    };
    
    explicit beta_distribution(RealType a = 2.0, RealType b = 2.0) : a_gamma(a), b_gamma(b) {}
    
    explicit beta_distribution(const param_type& param) : a_gamma(param.a()), b_gamma(param.b()) {}
    
    void reset() noexcept {}
    
    param_type param() const noexcept
    {
        return param_type(a(), b());
    }
    
    void param(const param_type& param)
    {
        a_gamma = gamma_dist_type(param.a());
        b_gamma = gamma_dist_type(param.b());
    }
    
    template <typename URNG>
    result_type operator()(URNG& engine)
    {
        return generate(engine, a_gamma, b_gamma);
    }
    
    template <typename URNG>
    result_type operator()(URNG& engine, const param_type& param)
    {
        gamma_dist_type a_param_gamma(param.a()),
        b_param_gamma(param.b());
        return generate(engine, a_param_gamma, b_param_gamma);
    }
    
    result_type min() const noexcept { return 0.0; }
    result_type max() const noexcept { return 1.0; }
    result_type a() const noexcept { return a_gamma.alpha(); }
    result_type b() const noexcept { return b_gamma.alpha(); }
    
    bool operator==(const beta_distribution<result_type>& other) const noexcept
    {
        return (param() == other.param() &&
                a_gamma == other.a_gamma &&
                b_gamma == other.b_gamma);
    }
    
    bool operator!=(const beta_distribution<result_type>& other) const noexcept
    {
        return !(*this == other);
    }

private:
    using gamma_dist_type = std::gamma_distribution<result_type>;
    
    gamma_dist_type a_gamma, b_gamma;
    
    template <typename URNG>
    result_type generate(URNG& engine,
                         gamma_dist_type& x_gamma,
                         gamma_dist_type& y_gamma)
    {
        result_type x = x_gamma(engine);
        return x / (x + y_gamma(engine));
    }
};

} // namespace std

#endif
