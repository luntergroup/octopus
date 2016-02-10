//
//  genotype_probability_map.h
//  Octopus
//
//  Created by Daniel Cooke on 04/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef genotype_probability_map_h
#define genotype_probability_map_h

#include <vector>
#include <unordered_map>
#include <functional>
#include <cstddef>

template <typename SampleTp, typename GenotypeTp>
class GenotypeProbabilityMap
{
public:
    GenotypeProbabilityMap() = default;
    ~GenotypeProbabilityMap() = default;
    
    GenotypeProbabilityMap(const GenotypeProbabilityMap&)            = default;
    GenotypeProbabilityMap& operator=(const GenotypeProbabilityMap&) = default;
    GenotypeProbabilityMap(GenotypeProbabilityMap&&)                 = default;
    GenotypeProbabilityMap& operator=(GenotypeProbabilityMap&&)      = default;
    
    double at(const SampleTp& sample, const GenotypeTp& genotype) const;
    
    bool empty() const noexcept;
    std::size_t size() const noexcept;
    
    void reserve(std::size_t n);
    
    void clear() noexcept;
    
    
private:
    std::vector<GenotypeTp> genotypes_;
    std::unordered_map<SampleTp, std::vector<double>> sample_genotype_probabilities_;
    std::unordered_map<std::reference_wrapper<const GenotypeTp>, std::size_t> genotype_indicies_;
};

#endif /* genotype_probability_map_h */
