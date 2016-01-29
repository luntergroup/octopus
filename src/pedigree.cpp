//
//  pedigree.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "pedigree.hpp"

// public methods

void Pedigree::clear()
{
    tree_.clear();
}

std::size_t Pedigree::size() const
{
    return boost::num_vertices(tree_);
}
