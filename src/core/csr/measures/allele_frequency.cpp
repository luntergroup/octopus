// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_frequency.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> AlleleFrequency::do_clone() const
{
    return std::make_unique<AlleleFrequency>(*this);
}

Measure::ResultType AlleleFrequency::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    boost::optional<double> result {};
    for (const auto& p : assignments) {
        std::vector<Allele> alleles; bool has_ref;
        std::tie(alleles, has_ref) = get_called_alleles(call, p.first, true);
        std::size_t read_count {0};
        std::vector<unsigned> allele_counts(alleles.size());
        for (const auto& h : p.second) {
            const auto& haplotype = h.first;
            const auto& reads = h.second;
            const auto haplotype_support_depth = count_overlapped(reads, call);
            if (haplotype_support_depth > 0) {
                std::transform(std::cbegin(alleles), std::cend(alleles), std::cbegin(allele_counts), std::begin(allele_counts),
                               [&] (const auto& allele, auto count) {
                                   if (haplotype.includes(allele)) {
                                       count += haplotype_support_depth;
                                   }
                                   return count;
                               });
                read_count += haplotype_support_depth;
            }
        }
        if (read_count > 0) {
            auto first_called_count_itr = std::cbegin(allele_counts);
            if (has_ref) ++first_called_count_itr;
            assert(first_called_count_itr != std::cend(allele_counts));
            const auto min_count_itr = std::min_element(first_called_count_itr, std::cend(allele_counts));
            const auto maf = static_cast<double>(*min_count_itr) / read_count;
            if (result) {
                result = std::min(*result, maf);
            } else {
                result = maf;
            }
        }
    }
    return result;
}

std::string AlleleFrequency::do_name() const
{
    return "AF";
}

std::vector<std::string> AlleleFrequency::do_requirements() const
{
    return {"ReadAssignments"};
}

} // namespace csr
} // namespace octopus
