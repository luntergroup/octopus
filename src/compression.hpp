//
//  compression.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_compression_hpp
#define Octopus_compression_hpp

#include <string>
#include <sstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/copy.hpp>

namespace Octopus
{
    inline std::string compress(const std::string& data)
    {
        std::stringstream decompressed {data};
        boost::iostreams::filtering_streambuf<boost::iostreams::input> stream;
        stream.push(boost::iostreams::zlib_compressor());
        stream.push(decompressed);
        std::stringstream compressed {};
        boost::iostreams::copy(stream, compressed);
        return compressed.str();
    }
    
    inline std::string decompress(const std::string& data)
    {
        std::stringstream compressed {data};
        boost::iostreams::filtering_streambuf<boost::iostreams::input> stream;
        stream.push(boost::iostreams::zlib_decompressor());
        stream.push(compressed);
        std::stringstream decompressed;
        boost::iostreams::copy(stream, decompressed);
        return decompressed.str();
    }
} // namespace Octopus

#endif
