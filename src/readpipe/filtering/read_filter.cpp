// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_filter.hpp"

namespace octopus { namespace readpipe
{

HasWellFormedCigar::HasWellFormedCigar() : BasicReadFilter {"HasWellFormedCigar"} {}

HasWellFormedCigar::HasWellFormedCigar(std::string name) : BasicReadFilter {std::move(name)} {}

bool is_all_clipped(const CigarString& cigar) noexcept
{
    return cigar.size() == 1 && is_clipping(cigar.front());
}

bool HasWellFormedCigar::passes(const AlignedRead& read) const noexcept
{
    return is_valid(read.cigar()) && is_minimal(read.cigar()) && !is_all_clipped(read.cigar());
}

HasValidBaseQualities::HasValidBaseQualities() : BasicReadFilter {"HasValidBaseQualities"} {}

HasValidBaseQualities::HasValidBaseQualities(std::string name) : BasicReadFilter {std::move(name)} {}

bool HasValidBaseQualities::passes(const AlignedRead& read) const noexcept
{
    return read.sequence().size() == read.base_qualities().size();
}

IsNotSecondaryAlignment::IsNotSecondaryAlignment()
: BasicReadFilter {"IsNotSecondaryAlignment"} {}

IsNotSecondaryAlignment::IsNotSecondaryAlignment(std::string name)
: BasicReadFilter {std::move(name)} {}

bool IsNotSecondaryAlignment::passes(const AlignedRead& read) const noexcept
{
    return !read.is_marked_secondary_alignment();
}

IsNotSupplementaryAlignment::IsNotSupplementaryAlignment()
: BasicReadFilter {"IsNotSupplementaryAlignment"} {}

IsNotSupplementaryAlignment::IsNotSupplementaryAlignment(std::string name)
: BasicReadFilter {std::move(name)} {}

bool IsNotSupplementaryAlignment::passes(const AlignedRead& read) const noexcept
{
    return !read.is_marked_supplementary_alignment();
}

IsGoodMappingQuality::IsGoodMappingQuality(MappingQuality good_mapping_quality)
:
BasicReadFilter {"IsGoodMappingQuality"}
, good_mapping_quality_ {good_mapping_quality} {}

IsGoodMappingQuality::IsGoodMappingQuality(std::string name, MappingQuality good_mapping_quality)
: BasicReadFilter {std::move(name)}
, good_mapping_quality_ {good_mapping_quality} {}

bool IsGoodMappingQuality::passes(const AlignedRead& read) const noexcept
{
    return read.mapping_quality() >= good_mapping_quality_;
}

HasSufficientGoodBaseFraction::HasSufficientGoodBaseFraction(BaseQuality good_base_quality,
                                                             double min_good_base_fraction)
: BasicReadFilter {"HasSufficientGoodBaseFraction"}
, good_base_quality_ {good_base_quality}
, min_good_base_fraction_ {min_good_base_fraction}
{}

HasSufficientGoodBaseFraction::HasSufficientGoodBaseFraction(std::string name,
                                                             BaseQuality good_base_quality,
                                                             double min_good_base_fraction)
: BasicReadFilter {std::move(name)}
, good_base_quality_ {good_base_quality}
, min_good_base_fraction_ {min_good_base_fraction}
{}

bool HasSufficientGoodBaseFraction::passes(const AlignedRead& read) const noexcept
{
    const auto& qualities = read.base_qualities();
    auto num_good_bases = std::count_if(std::cbegin(qualities), std::cend(qualities),
                                        [this] (const auto quality) {
                                            return quality >= good_base_quality_;
                                        });
    auto good_base_fraction = static_cast<double>(num_good_bases) / static_cast<double>(sequence_size(read));
    return good_base_fraction >= min_good_base_fraction_;
}

HasSufficientGoodQualityBases::HasSufficientGoodQualityBases(BaseQuality good_base_quality,
                                                             unsigned min_good_bases)
: BasicReadFilter {"HasSufficientGoodQualityBases"}
, good_base_quality_ {good_base_quality}, min_good_bases_ {min_good_bases} {}

HasSufficientGoodQualityBases::HasSufficientGoodQualityBases(std::string name,
                                                             BaseQuality good_base_quality,
                                                             unsigned min_good_bases)
: BasicReadFilter {std::move(name)}
, good_base_quality_ {good_base_quality}, min_good_bases_ {min_good_bases} {}

bool HasSufficientGoodQualityBases::passes(const AlignedRead& read) const noexcept
{
    const auto& qualities = read.base_qualities();
    return std::count_if(std::cbegin(qualities), std::cend(qualities), [this]
                         (const auto quality) {
                             return quality >= good_base_quality_;
                         }) >= min_good_bases_;
}

IsMapped::IsMapped() : BasicReadFilter {"IsMapped"} {}
IsMapped::IsMapped(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsMapped::passes(const AlignedRead& read) const noexcept
{
    return !read.is_marked_unmapped();
}

IsNotChimeric::IsNotChimeric() : BasicReadFilter {"IsNotChimeric"} {}
IsNotChimeric::IsNotChimeric(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsNotChimeric::passes(const AlignedRead& read) const noexcept
{
    return !read.has_other_segment();
}

IsNextSegmentMapped::IsNextSegmentMapped() : BasicReadFilter {"IsNextSegmentMapped"} {}
IsNextSegmentMapped::IsNextSegmentMapped(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsNextSegmentMapped::passes(const AlignedRead& read) const noexcept
{
    return !read.has_other_segment() || !read.next_segment().is_marked_unmapped();
}

IsNotMarkedDuplicate::IsNotMarkedDuplicate() : BasicReadFilter {"IsNotMarkedDuplicate"} {}
IsNotMarkedDuplicate::IsNotMarkedDuplicate(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsNotMarkedDuplicate::passes(const AlignedRead& read) const noexcept
{
    return !read.is_marked_duplicate();
}

IsShort::IsShort(Length max_length)
: BasicReadFilter {"IsShort"}
, max_length_ {max_length} {}

IsShort::IsShort(std::string name, Length max_length)
: BasicReadFilter {std::move(name)}
, max_length_ {max_length} {}

bool IsShort::passes(const AlignedRead& read) const noexcept
{
    return sequence_size(read) <= max_length_;
}

IsLong::IsLong(Length min_length)
: BasicReadFilter {"IsLong"}
, min_length_ {min_length} {}

IsLong::IsLong(std::string name, Length min_length)
: BasicReadFilter {std::move(name)}
, min_length_ {min_length} {}

bool IsLong::passes(const AlignedRead& read) const noexcept
{
    return sequence_size(read) >= min_length_;
}

IsNotContaminated::IsNotContaminated() : BasicReadFilter {"IsNotContaminated"} {}
IsNotContaminated::IsNotContaminated(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsNotContaminated::passes(const AlignedRead& read) const noexcept
{
    if (!read.has_other_segment() || read.next_segment().is_marked_unmapped()) {
        return true;
    }
    const auto template_length = read.next_segment().inferred_template_length();
    return template_length > region_size(read);
}

IsNotMarkedQcFail::IsNotMarkedQcFail() : BasicReadFilter {"IsNotMarkedQcFail"} {}
IsNotMarkedQcFail::IsNotMarkedQcFail(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsNotMarkedQcFail::passes(const AlignedRead& read) const noexcept
{
    return !read.is_marked_qc_fail();
}

IsProperTemplate::IsProperTemplate() : BasicReadFilter {"IsProperTemplate"} {}
IsProperTemplate::IsProperTemplate(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsProperTemplate::passes(const AlignedRead& read) const noexcept
{
    return !read.has_other_segment() || read.is_marked_all_segments_in_read_aligned();
}

IsLocalTemplate::IsLocalTemplate() : BasicReadFilter {"IsLocalTemplate"} {}
IsLocalTemplate::IsLocalTemplate(std::string name) :  BasicReadFilter {std::move(name)} {}

bool IsLocalTemplate::passes(const AlignedRead& read) const noexcept
{
    return !read.has_other_segment() || read.next_segment().contig_name() == contig_name(read);
}

NoUnlocalizedSupplementaryAlignments::NoUnlocalizedSupplementaryAlignments(boost::optional<MappingQuality> min_mapping_quality)
: NoUnlocalizedSupplementaryAlignments {"NoUnlocalizedSupplementaryAlignments", min_mapping_quality} {}
NoUnlocalizedSupplementaryAlignments::NoUnlocalizedSupplementaryAlignments(std::string name, boost::optional<MappingQuality> min_mapping_quality)
:  BasicReadFilter {std::move(name)}
, min_mapping_quality_ {min_mapping_quality}
{}

bool NoUnlocalizedSupplementaryAlignments::passes(const AlignedRead& read) const noexcept
{
    const auto supplementary_alignments = read.supplementary_alignments();
    if (min_mapping_quality_) {
        const auto is_good_unlocalized = [mq = *min_mapping_quality_] (const auto& alignment) noexcept {
            return alignment.mapping_quality() >= mq && is_unlocalized(alignment);};
        return std::none_of(std::cbegin(supplementary_alignments), std::cend(supplementary_alignments), is_good_unlocalized);
    } else {
        return std::none_of(std::cbegin(supplementary_alignments), std::cend(supplementary_alignments), is_unlocalized);
    }
}

NoUnplacedSupplementaryAlignments::NoUnplacedSupplementaryAlignments(boost::optional<MappingQuality> min_mapping_quality)
: NoUnplacedSupplementaryAlignments {"NoUnplacedSupplementaryAlignments", min_mapping_quality} {}
NoUnplacedSupplementaryAlignments::NoUnplacedSupplementaryAlignments(std::string name, boost::optional<MappingQuality> min_mapping_quality)
:  BasicReadFilter {std::move(name)}
, min_mapping_quality_ {min_mapping_quality}
{}

bool NoUnplacedSupplementaryAlignments::passes(const AlignedRead& read) const noexcept
{
    const auto supplementary_alignments = read.supplementary_alignments();
    if (min_mapping_quality_) {
        const auto is_good_unplaced = [mq = *min_mapping_quality_] (const auto& alignment) noexcept {
            return alignment.mapping_quality() >= mq && is_unplaced(alignment);};
        return std::none_of(std::cbegin(supplementary_alignments), std::cend(supplementary_alignments), is_good_unplaced);
    } else {
        return std::none_of(std::cbegin(supplementary_alignments), std::cend(supplementary_alignments), is_unplaced);
    }
}

NoDecoySupplementaryAlignments::NoDecoySupplementaryAlignments(boost::optional<MappingQuality> min_mapping_quality)
: NoDecoySupplementaryAlignments {"NoDecoySupplementaryAlignments", min_mapping_quality} {}
NoDecoySupplementaryAlignments::NoDecoySupplementaryAlignments(std::string name, boost::optional<MappingQuality> min_mapping_quality)
:  BasicReadFilter {std::move(name)}
, min_mapping_quality_ {min_mapping_quality}
{}

bool NoDecoySupplementaryAlignments::passes(const AlignedRead& read) const noexcept
{
    const auto supplementary_alignments = read.supplementary_alignments();
    if (min_mapping_quality_) {
        const auto is_good_decoy = [mq = *min_mapping_quality_] (const auto& alignment) noexcept {
            return alignment.mapping_quality() >= mq && is_decoy(alignment);};
        return std::none_of(std::cbegin(supplementary_alignments), std::cend(supplementary_alignments), is_good_decoy);
    } else {
        return std::none_of(std::cbegin(supplementary_alignments), std::cend(supplementary_alignments), is_decoy);
    }
}

} // namespace readpipe
} // namespace octopus
