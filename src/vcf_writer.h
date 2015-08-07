//
//  vcf_writer.h
//  Octopus
//
//  Created by Daniel Cooke on 29/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_writer__
#define __Octopus__vcf_writer__

class VcfWriter
{
public:
    VcfWriter()  = delete;
    ~VcfWriter() = default;
    
    VcfWriter(const VcfWriter&)            = default;
    VcfWriter& operator=(const VcfWriter&) = default;
    VcfWriter(VcfWriter&&)                 = default;
    VcfWriter& operator=(VcfWriter&&)      = default;
    
private:
    
};

#endif /* defined(__Octopus__vcf_writer__) */
