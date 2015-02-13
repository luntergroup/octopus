//
//  contig_name.cpp
//  Octopus
//
//  Created by Daniel Cooke on 13/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "contig_name.h"

// non-const static variables must be initialised in source file
std::unordered_map<uint_fast32_t, std::string> ContigName::contig_id_to_name_ {};
std::unordered_map<std::string, uint_fast32_t> ContigName::contig_name_to_id_ {};

void
ContigName::give_all_reference_contig_names(const std::unordered_set<std::string>& the_contig_names)
{
    ContigName::contig_id_to_name_.reserve(the_contig_names.size());
    uint_fast32_t contig_id {0};
    for (const auto& contig_name : the_contig_names) {
        ContigName::contig_id_to_name_.emplace(std::make_pair(contig_id, contig_name));
        ContigName::contig_name_to_id_.emplace(std::make_pair(contig_name, contig_id));
    }
}
