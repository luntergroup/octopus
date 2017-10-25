// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_frequency.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "core/tools/read_assigner.hpp"
#include "../facets/read_assignments.hpp"

#include <algorithm>
#include <iterator>

#include "core/types/allele.hpp"

namespace octopus { namespace csr {

auto get_called_alleles(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    const auto& call_region = mapped_region(call);
    auto genotype = get_genotype(call, sample);
    std::sort(std::begin(genotype), std::end(genotype));
    genotype.erase(std::unique(std::begin(genotype), std::end(genotype)), std::end(genotype));
    std::vector<Allele> result {};
    result.reserve(genotype.size());
    std::transform(std::cbegin(genotype), std::cend(genotype), std::back_inserter(result),
                   [&] (const auto& alt_seq) { return Allele {call_region, alt_seq}; });
    return result;
}

Measure::ResultType AlleleFrequency::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto assignments = boost::get<ReadAssignments::ResultType>(facets.at("ReadAssignments").get());
    double min_freq {1.0}, max_freq {0.0};
    for (const auto& p : assignments) {
        if (call.is_heterozygous(p.first)) {
            auto alleles = get_called_alleles(call, p.first);
            auto allele_support = compute_allele_support(alleles, p.second);
            std::vector<double> allele_counts(alleles.size());
            std::size_t read_count {0};
            std::transform(std::cbegin(allele_support), std::cend(allele_support), std::begin(allele_counts),
                           [&] (const auto& p) {
                               read_count += p.second.size();
                               return p.second.size();
                           });
            std::sort(std::begin(allele_counts), std::end(allele_counts), std::greater<> {});
            auto major_af = allele_counts.front() / read_count;
            auto minor_af = allele_counts.back() / read_count;
            if (major_af > max_freq) {
                max_freq = major_af;
            }
            if (minor_af < min_freq) {
                min_freq = minor_af;
            }
        }
    }
    return std::max(min_freq, 1.0 - max_freq);
}

std::string AlleleFrequency::do_name() const
{
    return "allele frequency";
}

std::vector<std::string> AlleleFrequency::do_requirements() const
{
    return {"ReadAssignments"};
}

} // namespace csr
} // namespace octopus
