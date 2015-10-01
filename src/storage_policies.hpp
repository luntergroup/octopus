//
//  storage_policies.hpp
//  Octopus
//
//  Created by Daniel Cooke on 24/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_storage_policies_hpp
#define Octopus_storage_policies_hpp

#include <string>
#include <unordered_set>
#include <boost/utility/string_ref.hpp>

namespace policies {
    
    class StoreStringCopy
    {
    protected:
        using InputType     = const std::string&;
        using ValueType     = std::string;
        using ReferenceType = boost::string_ref;
        
        StoreStringCopy() = default;
        
        ReferenceType store(InputType value)
        {
            auto copy_it = copies_.emplace(value);
            return *copy_it.first;
        }
        
    private:
        std::unordered_set<ValueType> copies_;
    };
    
    class StoreStringReference
    {
    protected:
        using InputType     = boost::string_ref;
        using ValueType     = boost::string_ref;
        using ReferenceType = boost::string_ref;
        
        StoreStringReference() = default;
        
        ReferenceType store(InputType value)
        {
            return value;
        }
    };
    
} // end namespace policies

#endif
