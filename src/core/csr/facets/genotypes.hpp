// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genotypes_hpp
#define genotypes_hpp

#include <string>
#include <functional>

#include "facet.hpp"

namespace octopus { namespace csr {

class Genotypes : public Facet
{
public:
    using ResultType = std::reference_wrapper<const GenotypeMap>;
    
    Genotypes() = default;
    Genotypes(GenotypeMap genotypes);

private:
    static const std::string name_;
    
    GenotypeMap genotypes_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
