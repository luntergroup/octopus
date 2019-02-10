// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef csr_pedigree_hpp
#define csr_pedigree_hpp

#include <string>
#include <vector>
#include <functional>

#include "basics/pedigree.hpp"
#include "facet.hpp"

namespace octopus { namespace csr {

class Pedigree : public Facet
{
public:
    using ResultType = std::reference_wrapper<const octopus::Pedigree>;
    
    Pedigree() = default;
    Pedigree(octopus::Pedigree pedigree);

private:
    static const std::string name_;
    
    octopus::Pedigree pedigree_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
