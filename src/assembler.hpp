//
//  assembler.hpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef assembler_hpp
#define assembler_hpp

#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
#include <iterator>
#include <cstddef>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/functional/hash.hpp>

class Assembler
{
public:
    using SequenceType = std::string;
    
    Assembler() = delete;
    
    Assembler(std::size_t kmer_size);
    
    ~Assembler() = default;
    
    Assembler(const Assembler&)            = default;
    Assembler& operator=(const Assembler&) = default;
    Assembler(Assembler&&)                 = default;
    Assembler& operator=(Assembler&&)      = default;
    
    void insert(const SequenceType& sequence);
    void insert(const SequenceType& sequence, std::size_t position);
    
    
    
    bool empty() const noexcept;
    
    void clear();
    
private:
    using SequenceIterator = SequenceType::const_iterator;
    
    struct KmerNode
    {
        SequenceIterator sequence;
        std::vector<std::size_t> positions;
    };
    
    using KmerGraph = boost::adjacency_list<boost::listS, boost::listS, boost::directedS, KmerNode>;
    
    std::size_t k_;
};

#endif /* assembler_hpp */
