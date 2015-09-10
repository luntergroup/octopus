//
//  hash_functions.h
//  Octopus
//
//  Created by Daniel Cooke on 09/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_hash_functions_h
#define Octopus_hash_functions_h

#include <boost/utility/string_ref.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem/path.hpp>

namespace std
{
    template<>
    struct hash<boost::string_ref> {
        size_t operator()(const boost::string_ref& sr) const {
            return boost::hash_range(sr.begin(), sr.end());
        }
    };
}

namespace std {
    template <> struct hash<boost::filesystem::path>
    {
        size_t operator()(const boost::filesystem::path& path) const
        {
            return hash<string>()(path.string());
        }
    };
}

#endif
