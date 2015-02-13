//
//  contig_name.h
//  Octopus
//
//  Created by Daniel Cooke on 13/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_contig_name_h
#define Octopus_contig_name_h

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>

#include "equitable.h"

using std::uint_fast32_t;

class ContigName : Equitable<ContigName>
{
public:
    ContigName() = delete;
    ContigName(const std::string& contig_name);
    
    static void give_all_reference_contig_names(const std::unordered_set<std::string>& the_contig_names);
    
    std::string get_name() const;
private:
    static std::unordered_map<uint_fast32_t, std::string> contig_id_to_name_;
    static std::unordered_map<std::string, uint_fast32_t> contig_name_to_id_;
    uint_fast32_t contig_id_;
};

inline ContigName::ContigName(const std::string& contig_name)
:   contig_id_ {contig_name_to_id_.at(contig_name)}
{}

inline std::string ContigName::get_name() const{
    return contig_id_to_name_.at(contig_id_);
}

inline std::string to_string(const ContigName& a_contig)
{
    return a_contig.get_name();
}

inline bool operator==(const ContigName& lhs, const ContigName& rhs)
{
    return lhs.get_name() == rhs.get_name();
}

namespace std {
    template <> struct hash<ContigName>
    {
        size_t operator()(const ContigName& a_contig) const
        {
            return hash<std::string>()(to_string(a_contig));
        }
    };
}

inline
std::ostream& operator<<(std::ostream& os, const ContigName& a_contig)
{
    os << to_string(a_contig);
    return os;
}

#endif
