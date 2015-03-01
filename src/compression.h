//
//  compression.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_compression_h
#define Octopus_compression_h

#include <string>
#include <sstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

namespace policies {
    
    class ZlibCompression
    {
    public:
        template <typename T>
        std::string compress(const T& data)
        {
            std::stringstream decompressed {data};
            boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
            out.push(boost::iostreams::zlib_compressor());
            out.push(decompressed);
            std::stringstream compressed {};
            boost::iostreams::copy(out, compressed);
            return compressed.str();
        }
        
        template <typename T>
        std::string decompress(const T& data)
        {
            std::stringstream compressed {data};
            boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
            in.push(boost::iostreams::zlib_decompressor());
            in.push(compressed);
            std::stringstream decompressed;
            boost::iostreams::copy(in, decompressed);
            return decompressed.str();
        }
    };
}

#endif
