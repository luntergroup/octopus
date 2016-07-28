//
//  i_candidate_variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_i_candidate_variant_generator__
#define Octopus_i_candidate_variant_generator__

#include <vector>
#include <cstddef>
#include <type_traits>
#include <iostream>
#include <string>

#include "common.hpp"
#include "variant.hpp"
#include "mappable_flat_multi_set.hpp"

class AlignedRead;
class GenomicRegion;

namespace octopus
{
    class ICandidateVariantGenerator
    {
    public:
        virtual ~ICandidateVariantGenerator() = default;
        
        virtual std::vector<Variant> generate_candidates(const GenomicRegion& region) = 0;
        
        virtual bool requires_reads() const noexcept { return false; };
        
        virtual void add_read(const AlignedRead&) {};
        
        // add_reads is not strictly necessary as the effect of calling add_reads must be the same as
        // calling add_read for each read. However, there may be significant performance benifits
        // to having an add_reads method to avoid many virtual dispatches.
        // Ideally add_reads would be templated to accept any InputIterator, but it is not possible
        // to have template virtual methods. The best solution is therefore to just overload add_reads
        // for common container iterators, more can easily be added if needed.
        virtual void add_reads(std::vector<AlignedRead>::const_iterator first,
                               std::vector<AlignedRead>::const_iterator last) {};
        virtual void add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                               MappableFlatMultiSet<AlignedRead>::const_iterator last) {};
        
        virtual void reserve(size_t n) {};
        virtual void clear() {};
    };
    
    namespace detail
    {
        template <typename Container, typename G>
        void add_reads(const Container& reads, G& generator, std::true_type)
        {
            generator.add_reads(std::cbegin(reads), std::cend(reads));
        }
        
        template <typename ReadMap, typename G>
        void add_reads(const ReadMap& reads, G& generator, std::false_type)
        {
            for (const auto& p : reads) {
                generator.add_reads(std::cbegin(p.second), std::cend(p.second));
            }
        }
    } // namespace detail
    
    template <typename Container, typename G,
              typename = std::enable_if_t<std::is_base_of<ICandidateVariantGenerator, G>::value>>
    void add_reads(const Container& reads, G& generator)
    {
        using ValueType = typename std::decay_t<typename Container::value_type>;
        detail::add_reads(reads, generator, std::is_same<ValueType, AlignedRead> {});
    }
    
    namespace debug
    {
        template <typename S, typename Container>
        void print_generated_candidates(S&& stream, const Container& candidates,
                                        const std::string& generator_name)
        {
            if (candidates.empty()) {
                stream << "No candidates generated from " << generator_name << '\n';
            } else {
                stream << "Generated " << candidates.size();
                stream << " candidate";
                if (candidates.size() > 1) {
                    stream << "s";
                }
                stream << " from " << generator_name << ":\n";
                for (const auto& c : candidates) stream << c << '\n';
            }
        }
        
        template <typename Container>
        void print_generated_candidates(const Container& candidates,
                                        const std::string& generator_name)
        {
            print_generated_candidates(std::cout, candidates, generator_name);
        }
        
    } // namespace debug
} // namespace octopus

#endif
