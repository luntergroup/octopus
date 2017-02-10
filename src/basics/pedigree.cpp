// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "pedigree.hpp"

namespace octopus {

// public methods

void Pedigree::clear()
{
    tree_.clear();
}

std::size_t Pedigree::size() const
{
    return boost::num_vertices(tree_);
}
    
} // namespace octopus