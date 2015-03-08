//
//  read_filter_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_filter_factory_h
#define Octopus_read_filter_factory_h

class ReadFilter;

class ReadFilterFactory
{
public:
    ReadFilterFactory() = default;
    ~ReadFilterFactory() = default;
    
    ReadFilterFactory(const ReadFilterFactory&)            = default;
    ReadFilterFactory& operator=(const ReadFilterFactory&) = default;
    ReadFilterFactory(ReadFilterFactory&&)                 = default;
    ReadFilterFactory& operator=(ReadFilterFactory&&)      = default;
    
    ReadFilter make(/* options in here */) const;
};

#endif
