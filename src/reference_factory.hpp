//
//  reference_factory.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef reference_factory_hpp
#define reference_factory_hpp

#include <memory>

#include <boost/filesystem/path.hpp>

#include "reference_genome_impl.hpp"
#include "reference_genome.hpp"

class ReferenceFactory
{
public:
    ReferenceFactory(boost::filesystem::path path);
    
    ~ReferenceFactory() = default;
    
    ReferenceFactory(const ReferenceFactory&)            = default;
    ReferenceFactory& operator=(const ReferenceFactory&) = default;
    ReferenceFactory(ReferenceFactory&&)                 = default;
    ReferenceFactory& operator=(ReferenceFactory&&)      = default;
    
    ReferenceGenome make() const;
    
private:
    std::unique_ptr<ReferenceGenomeImpl> impl_;
};

#endif /* reference_factory_hpp */
