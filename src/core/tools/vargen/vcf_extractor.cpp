// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "vcf_extractor.hpp"

#include <deque>
#include <algorithm>
#include <iterator>
#include <utility>

#include "io/variant/vcf_spec.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/sequence_utils.hpp"

namespace octopus { namespace coretools {

VcfExtractor::VcfExtractor(std::unique_ptr<const VcfReader> reader)
: VcfExtractor {std::move(reader), Options {}}
{}

VcfExtractor::VcfExtractor(std::unique_ptr<const VcfReader> reader, Options options)
: reader_ {std::move(reader)}
, options_ {options}
{}

std::unique_ptr<VariantGenerator> VcfExtractor::do_clone() const
{
    return std::make_unique<VcfExtractor>(*this);
}

namespace {

static bool is_canonical(const VcfRecord::NucleotideSequence& allele)
{
    return allele != vcfspec::missingValue
           && std::none_of(std::cbegin(allele), std::cend(allele),
                           [](const auto base) { return base == vcfspec::deletedBase; });
}

bool is_good_quality(const VcfRecord& record, boost::optional<VcfRecord::QualityType> min_quality) noexcept
{
    return !min_quality || (record.qual() && *record.qual() >= *min_quality);
}

template <typename Iterator>
auto make_allele(const Iterator first_base, const Iterator last_base)
{
    Variant::NucleotideSequence result{first_base, last_base};
    utils::capitalise(result);
    return result;
}

template <typename Container>
void extract_variants(const VcfRecord& record, Container& result)
{
    for (const auto& alt_allele : record.alt()) {
        if (is_canonical(alt_allele)) {
            const auto& ref_allele = record.ref();
            if (ref_allele.size() != alt_allele.size()) {
                auto begin = record.pos();
                const auto p = std::mismatch(std::cbegin(ref_allele), std::cend(ref_allele),
                                             std::cbegin(alt_allele), std::cend(alt_allele));
                if (p.first != std::cend(ref_allele) && alt_allele.size() > ref_allele.size()) {
                    // Split non-reference padded insertions into snv (or mnv) and insertion with empty
                    // reference (e.g. A -> TT makes two variants A -> T and -> T).
                    const auto ref_pad_size = std::distance(std::cbegin(ref_allele), p.first);
                    begin += ref_pad_size;
                    const auto remaining_ref_size = ref_allele.size() - ref_pad_size;
                    const auto first_alt_end = std::next(p.second, remaining_ref_size);
                    result.emplace_back(record.chrom(), begin - 1,
                                        make_allele(p.first, std::cend(ref_allele)),
                                        make_allele(p.second, first_alt_end));
                    begin += remaining_ref_size;
                    result.emplace_back(record.chrom(), begin - 1, "",
                                        make_allele(first_alt_end, std::cend(alt_allele)));
                } else {
                    begin += std::distance(std::cbegin(ref_allele), p.first);
                    result.emplace_back(record.chrom(), begin - 1,
                                        make_allele(p.first, std::cend(ref_allele)),
                                        make_allele(p.second, std::cend(alt_allele)));
                }
            } else {
                using utils::capitalise_copy;
                result.emplace_back(record.chrom(), record.pos() - 1,
                                    capitalise_copy(record.ref()),
                                    capitalise_copy(alt_allele));
            }
        }
    }
}

std::vector<Variant> fetch_variants(const GenomicRegion& region, const VcfReader& reader,
                                    const boost::optional<VcfRecord::QualityType> min_quality)
{
    std::deque<Variant> variants{}; // Use deque to prevent reallocating
    auto p = reader.iterate(region, VcfReader::UnpackPolicy::sites);
    std::for_each(std::move(p.first), std::move(p.second),
                  [&variants, min_quality](const auto& record) {
                      if (is_good_quality(record, min_quality)) {
                          extract_variants(record, variants);
                      }
                  });
    std::vector<Variant> result{std::make_move_iterator(std::begin(variants)),
                                std::make_move_iterator(std::end(variants))};
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

} // namespace

std::vector<Variant> VcfExtractor::do_generate_variants(const GenomicRegion& region)
{
    return fetch_variants(region, *reader_, options_.min_quality);
}

std::string VcfExtractor::name() const
{
    return "VCF extraction";
}

} // namespace coretools
} // namespace octopus
