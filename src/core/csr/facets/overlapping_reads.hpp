// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef overlapping_reads_hpp
#define overlapping_reads_hpp

#include <string>

#include "facet.hpp"
#include "config/common.hpp"
#include "basics/aligned_read.hpp"

namespace octopus { namespace csr {

class OverlappingReads : public Facet
{
public:
    using ResultType = ReadMap;
    
    OverlappingReads(ReadMap reads);

private:
    static const std::string name_;
    
    ReadMap reads_;
    
    const std::string& do_name() const noexcept { return name_; }
    Facet::ResultType do_get() const;
};

} // namespace csr
} // namespace octopus

#endif
