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

} // namespace

std::vector<Variant> VcfExtractor::do_generate_variants(const GenomicRegion& region)
{
    std::deque<Variant> variants {};
    for (auto p = reader_->iterate(region, VcfReader::UnpackPolicy::sites); p.first != p.second; ++p.first) {
        if (is_good(*p.first)) {
            extract_variants(*p.first, variants);
        }
    }
    std::vector<Variant> result {std::make_move_iterator(std::begin(variants)),
                                 std::make_move_iterator(std::end(variants))};
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

std::string VcfExtractor::name() const
{
    return "VCF extraction";
}

bool VcfExtractor::is_good(const VcfRecord& record)
{
    if (!options_.extract_filtered && is_filtered(record)) return false;
    return !options_.min_quality || (record.qual() && *record.qual() >= *options_.min_quality);
}

} // namespace coretools
} // namespace octopus
