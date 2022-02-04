// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "vcf_record.hpp"

#include <algorithm>
#include <iterator>

#include <boost/lexical_cast.hpp>

#include "io/reference/reference_genome.hpp"
#include "vcf_spec.hpp"

namespace octopus {

// public methods

const GenomicRegion& VcfRecord::mapped_region() const noexcept
{
    return region_;
}

const GenomicRegion::ContigName& VcfRecord::chrom() const noexcept
{
    return region_.contig_name();
}

GenomicRegion::Position VcfRecord::pos() const noexcept
{
    return region_.begin() + 1;
}

const std::string& VcfRecord::id() const noexcept
{
    return id_;
}

const VcfRecord::NucleotideSequence& VcfRecord::ref() const noexcept
{
    return ref_;
}

unsigned VcfRecord::num_alt() const noexcept
{
    return static_cast<unsigned>(alt_.size());
}

const std::vector<VcfRecord::NucleotideSequence>& VcfRecord::alt() const noexcept
{
    return alt_;
}

boost::optional<VcfRecord::QualityType> VcfRecord::qual() const noexcept
{
    return qual_;
}

bool VcfRecord::has_filter(const KeyType& filter) const noexcept
{
    return std::find(std::cbegin(filter_), std::cend(filter_), filter) != std::cend(filter_);
}

const std::vector<VcfRecord::KeyType>& VcfRecord::filter() const noexcept
{
    return filter_;
}

bool VcfRecord::has_info(const KeyType& key) const noexcept
{
    return info_.count(key) == 1;
}

std::vector<VcfRecord::KeyType> VcfRecord::info_keys() const
{
    std::vector<KeyType> result {};
    result.reserve(info_.size());
    std::transform(info_.cbegin(), info_.cend(), std::back_inserter(result), [] (const auto& p) {
        return p.first;
    });
    return result;
}

const std::vector<VcfRecord::ValueType>& VcfRecord::info_value(const KeyType& key) const
{
    return info_.at(key);
}

bool VcfRecord::has_format(const KeyType& key) const noexcept
{
    return std::find(std::cbegin(format_), std::cend(format_), key) != std::cend(format_);
}

boost::optional<unsigned> VcfRecord::format_cardinality(const KeyType& key) const noexcept
{
    boost::optional<unsigned> result {};
    if (has_format(key)) {
        for (const auto& p : samples_) {
            const auto sample_format_cardinality = p.second.other.at(key).size();
            if (result) {
                if (*result != sample_format_cardinality) return boost::none;
            } else {
                result = sample_format_cardinality;
            }
        }
    }
    return result;
}

const std::vector<VcfRecord::KeyType>& VcfRecord::format() const noexcept
{
    return format_;
}

unsigned VcfRecord::num_samples() const noexcept
{
    return samples_.size();
}

bool VcfRecord::has_genotypes() const noexcept
{
    return std::find(std::cbegin(format_), std::cend(format_), vcfspec::format::genotype) != std::cend(format_);
}

unsigned VcfRecord::ploidy(const SampleName& sample) const
{
    return get_genotype(sample).indices.size();
}

bool VcfRecord::is_sample_phased(const SampleName& sample) const
{
    return get_genotype(sample).phased;
}

bool VcfRecord::is_homozygous(const SampleName& sample) const
{
    const auto& genotype = get_genotype(sample).indices;
    return std::adjacent_find(std::cbegin(genotype), std::cend(genotype), std::not_equal_to<>{}) == std::cend(genotype);
}

bool VcfRecord::is_heterozygous(const SampleName& sample) const
{
    return !is_homozygous(sample);
}

bool VcfRecord::is_homozygous_ref(const SampleName& sample) const
{
    const auto& genotype = get_genotype(sample).indices;
    return std::all_of(std::cbegin(genotype), std::cend(genotype),
                       [] (const auto& allele) { return allele == 0; });
}

bool VcfRecord::is_refcall() const
{
    return alt_.size() == 1 && alt_.front() == vcfspec::allele::nonref;
}

bool VcfRecord::is_homozygous_non_ref(const SampleName& sample) const
{
    const auto& genotype = get_genotype(sample).indices;
    return genotype.front() > 0 && is_homozygous(sample);
}

bool VcfRecord::has_ref_allele(const SampleName& sample) const
{
    const auto& genotype = get_genotype(sample).indices;
    return std::find(std::cbegin(genotype), std::cend(genotype), 0) != std::cend(genotype);
}

bool VcfRecord::has_alt_allele(const SampleName& sample) const
{
    const auto& genotype = get_genotype(sample).indices;
    return std::find_if_not(std::cbegin(genotype), std::cend(genotype),
                            [] (const auto& allele) {
                                return allele == 0;
                            }) != std::cend(genotype);
}

const std::vector<VcfRecord::AlleleIndex>& VcfRecord::genotype(const SampleName& sample) const
{
    return get_genotype(sample).indices;
}

const std::vector<VcfRecord::ValueType>& VcfRecord::get_sample_value(const SampleName& sample, const KeyType& key) const
{
    return samples_.at(sample).other.at(key);
}

// helper non-members needed for printing

namespace {

template <typename T>
std::ostream& print(std::ostream& os, const std::vector<T>& v, const std::string& delim = ",",
                    const std::string& empty_value = ".")
{
    if (v.empty()) {
        os << empty_value;
    } else {
        std::copy(std::cbegin(v), std::prev(std::cend(v)),
                  std::ostream_iterator<T>(os, delim.c_str()));
        os << v.back();
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    return print(os, v);
}

} // namespace

// private methods

std::vector<VcfRecord::SampleName> VcfRecord::samples() const
{
    std::vector<SampleName> result {};
    result.reserve(samples_.size());
    for (const auto& p : samples_) {
        result.push_back(p.first);
    }
    return result;
}

const VcfRecord::Genotype& VcfRecord::get_genotype(const SampleName& sample) const
{
    return *samples_.at(sample).genotype;
}

std::string VcfRecord::get_allele_number(const NucleotideSequence& allele) const
{
    if (allele == vcfspec::missingValue) {
        return vcfspec::missingValue;
    } else if (allele == ref_) {
        return "0";
    } else {
        const auto it = std::find(std::cbegin(alt_), std::cend(alt_), allele);
        return std::to_string(std::distance(std::cbegin(alt_), it) + 1);
    }
}

void VcfRecord::print_info(std::ostream& os) const
{
    if (info_.empty()) {
        os << vcfspec::missingValue;
    } else {
        auto last = std::next(std::cbegin(info_), info_.size() - 1);
        std::for_each(std::cbegin(info_), last,
                      [&os] (const auto& p) {
                          os << p.first;
                          if (!p.second.empty()) {
                              os << "=" << p.second;
                          }
                          os << ';';
                      });
        os << last->first;
        if (!last->second.empty()) {
            os << "=" << last->second;
        }
    }
}

void VcfRecord::print_genotype_allele_numbers(std::ostream& os, const SampleName& sample) const
{
    std::vector<std::string> allele_numbers(ploidy(sample));
    const auto& genotype = get_genotype(sample);
    std::transform(std::cbegin(genotype.indices), std::cend(genotype.indices), std::begin(allele_numbers),
                   [] (auto number) -> std::string { 
                       if (number < 0) {
                           return std::to_string(number); 
                        } else {
                            return vcfspec::missingValue;
                   }});
    print(os, allele_numbers, (genotype.phased) ? "|" : "/");
}

void VcfRecord::print_other_sample_data(std::ostream& os, const SampleName& sample) const
{
    if (!samples_.empty()) {
        if (samples_.at(sample).other.empty()) {
            os << vcfspec::missingValue;
        } else {
            const auto& data = samples_.at(sample).other;
            auto last = std::next(cbegin(data), data.size() - 1);
            std::for_each(std::cbegin(data), last, [&os] (const auto& p) {
                print(os, p.second, ",");
                os << ":";
            });
            print(os, last->second, ",");
        }
    }
}

void VcfRecord::print_sample_data(std::ostream& os) const
{
    if (num_samples() > 0) {
        print(os, format_, ":");
        os << '\t';
        auto samples = this->samples();
        std::for_each(std::cbegin(samples), std::prev(std::cend(samples)),
                      [this, &os] (const SampleName& sample) {
                          auto it = std::cbegin(format_);
                          if (*it == vcfspec::format::genotype) {
                              print_genotype_allele_numbers(os, sample);
                              ++it;
                          }
                          std::for_each(it, std::cend(format_),
                                        [this, &os, &sample] (const KeyType& key) {
                                            os << ':';
                                            print(os, get_sample_value(sample, key), ",");
                                        });
                          os << '\t';
        });
        auto it = std::cbegin(format_);
        if (*it == vcfspec::format::genotype) {
            print_genotype_allele_numbers(os, samples.back());
            ++it;
        }
        std::for_each(it, std::cend(format_),
                      [this, &os, &samples] (const KeyType& key) {
                          os << ':';
                          print(os, get_sample_value(samples.back(), key), ",");
                      });
    }
}

// non-member functions

const VcfRecord::NucleotideSequence& get_allele(const VcfRecord& record, const VcfRecord::AlleleIndex index)
{
    const static std::string missing_value {vcfspec::missingValue};
    if (index < 0) return missing_value;
    if (index == 0) return record.ref();
    return record.alt()[index - 1];
}

std::vector<VcfRecord::NucleotideSequence> get_genotype(const VcfRecord& record, const VcfRecord::SampleName& sample)
{
    const auto& gt = record.genotype(sample);
    std::vector<VcfRecord::NucleotideSequence> result(gt.size());
    std::transform(std::cbegin(gt), std::cend(gt), std::begin(result),
                   [&] (auto index) { return get_allele(record, index); });
    return result;
}

bool is_missing(const std::vector<VcfRecord::ValueType>& values) noexcept
{
    return values.size() < 2 && values.front() == vcfspec::missingValue;
}

bool is_info_missing(const VcfRecord::KeyType& key, const VcfRecord& record)
{
    return !record.has_info(key) || is_missing(record.info_value(key));
}

bool is_refcall(const VcfRecord& record)
{
    return record.is_refcall();
}

bool is_filtered(const VcfRecord& record) noexcept
{
    const auto& filters = record.filter();
    return !filters.empty() && !(filters[0] == vcfspec::filter::pass || filters[0] == vcfspec::missingValue);
}

bool is_dbsnp_member(const VcfRecord& record) noexcept
{
    return record.has_info(vcfspec::info::dbSNPMember);
}

bool is_hapmap2_member(const VcfRecord& record) noexcept
{
    return record.has_info(vcfspec::info::hapmap2Member);
}

bool is_hapmap3_member(const VcfRecord& record) noexcept
{
    return record.has_info(vcfspec::info::hapmap3Member);
}

bool is_1000g_member(const VcfRecord& record) noexcept
{
    return record.has_info(vcfspec::info::thousandGenomes);
}

bool is_somatic(const VcfRecord& record) noexcept
{
    return record.has_info(vcfspec::info::somatic);
}

bool is_validated(const VcfRecord& record) noexcept
{
    return record.has_info(vcfspec::info::validated);
}

boost::optional<GenomicRegion> get_phase_region(const VcfRecord& record, const VcfRecord::SampleName& sample)
{
    if (record.is_sample_phased(sample) && record.has_format(vcfspec::format::phaseSet)) {
        return GenomicRegion {
        record.chrom(),
        boost::lexical_cast<ContigRegion::Position>(record.get_sample_value(sample, vcfspec::format::phaseSet).front()) - 1,
        static_cast<ContigRegion::Position>(record.pos() + record.ref().size()) - 1
        };
    } else {
        return boost::none;
    }
}

bool operator==(const VcfRecord& lhs, const VcfRecord& rhs)
{
    // TODO: this should really compare other fields
    return mapped_region(lhs) == mapped_region(rhs) && lhs.ref() == rhs.ref() && lhs.alt() == rhs.alt();
}

bool operator<(const VcfRecord& lhs, const VcfRecord& rhs)
{
    if (mapped_region(lhs) == mapped_region(rhs)) {
        if (lhs.ref() == rhs.ref()) {
            return lhs.alt() < rhs.alt();
        } else {
            return lhs.ref() < rhs.ref();
        }
    } else {
        return mapped_region(lhs) < mapped_region(rhs);
    }
}

std::ostream& operator<<(std::ostream& os, const VcfRecord& record)
{
    os << record.chrom() << "\t";
    os << record.pos() << "\t";
    os << record.id_ << "\t";
    os << record.ref_ << "\t";
    os << record.alt_ << "\t";
    if (record.qual_) {
        os << static_cast<float>(*record.qual_) << "\t";
    } else {
        os << vcfspec::missingValue << "\t";
    }
    os << record.filter_ << "\t";
    record.print_info(os);
    os << "\t";
    record.print_sample_data(os);
    return os;
}

// VcfRecord::Builder

VcfRecord::Builder::Builder(const ReferenceGenome& reference)
: chrom_ {}
, pos_ {}
, id_ {}
, ref_ {}
, alt_ {}
, qual_ {}
, filter_ {}
, info_ {}
, format_ {}
, samples_ {}
, reference_ {std::addressof(reference)}
{}

VcfRecord::Builder::Builder(const VcfRecord& call)
: chrom_ {call.chrom()}
, pos_ {call.pos()}
, id_ {call.id()}
, ref_ {call.ref()}
, alt_ {call.alt()}
, qual_ {call.qual()}
, filter_ {call.filter()}
, info_ {call.info_}
, format_ {call.format()}
, samples_ {call.samples_}
, reference_ {}
{}

VcfRecord::Builder& VcfRecord::Builder::set_chrom(std::string name)
{
    chrom_ = std::move(name);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_pos(GenomicRegion::Position pos)
{
    pos_ = pos;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_id(std::string id)
{
    id_ = std::move(id);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_ref(const char allele)
{
    ref_ = allele;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_ref(NucleotideSequence allele)
{
    ref_ = std::move(allele);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_alt(const char allele)
{
    alt_.resize(1);
    alt_.front() = allele;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_alt(NucleotideSequence allele)
{
    alt_.resize(1);
    alt_.front() = std::move(allele);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_alt(std::vector<NucleotideSequence> alleles)
{
    alt_ = std::move(alleles);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_qual(QualityType quality)
{
    qual_ = quality;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_passed()
{
    filter_.assign({vcfspec::filter::pass});
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_filter(std::vector<KeyType> filter)
{
    filter_ = std::move(filter);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_filter(std::initializer_list<KeyType> filter)
{
    filter_ = filter;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::add_filter(KeyType filter)
{
    filter_.push_back(std::move(filter));
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::clear_filter() noexcept
{
    filter_.clear();
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::reserve_info(unsigned n)
{
    info_.reserve(n);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::add_info(const KeyType& key)
{
    info_.emplace(key, std::vector<ValueType> {});
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_info(const KeyType& key, const ValueType& value)
{
    return this->set_info(key, {value});
}

VcfRecord::Builder& VcfRecord::Builder::set_info(const KeyType& key, std::vector<ValueType> values)
{
    if (key == "END") {
        if (values.size() != 1)
            throw std::runtime_error {"VcfRecord::Builder INFO key END requires 1 value"};
        end_ = std::stoll(values.front());
    }
    info_[key] = std::move(values);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_info(const KeyType& key, std::initializer_list<ValueType> values)
{
    return this->set_info(key, std::vector<ValueType> {values});
}

VcfRecord::Builder& VcfRecord::Builder::set_info_flag(KeyType key)
{
    return this->set_info(std::move(key), {});
}

VcfRecord::Builder& VcfRecord::Builder::set_info_missing(const KeyType& key)
{
    return this->set_info(key, {vcfspec::missingValue});
}

VcfRecord::Builder& VcfRecord::Builder::clear_info() noexcept
{
    info_.clear();
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::clear_info(const KeyType& key)
{
    info_.erase(key);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_format(std::vector<KeyType> format)
{
    format_ = std::move(format);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_format(std::initializer_list<KeyType> format)
{
    format_ = std::move(format);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::add_format(KeyType key)
{
    format_.push_back(std::move(key));
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::add_format(std::initializer_list<KeyType> keys)
{
    format_.insert(std::cend(format_), keys);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::reserve_samples(unsigned n)
{
    samples_.reserve(n);
    return *this;
}

VcfRecord::Builder&VcfRecord::Builder:: set_homozygous_ref_genotype(const SampleName& sample, const unsigned ploidy)
{
    auto& genotype = samples_[sample].genotype;
    genotype = VcfRecord::Genotype {};
    genotype->indices.resize(ploidy, 0);
    genotype->phased = true;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_genotype(const SampleName& sample, const std::vector<NucleotideSequence>& alleles,
                                                     const Phasing phasing)
{
    auto& genotype = samples_[sample].genotype;
    genotype = VcfRecord::Genotype {};
    genotype->indices.resize(alleles.size());
    std::transform(std::cbegin(alleles), std::cend(alleles), std::begin(genotype->indices),
                   [this] (const auto& allele) -> VcfRecord::AlleleIndex {
                       if (allele == vcfspec::missingValue) {
                           return -1;
                       } else if (allele == ref_) {
                           return 0;
                       } else {
                           const auto itr = std::find(std::cbegin(alt_), std::cend(alt_), allele);
                           if (itr != std::cend(alt_)) {
                               return std::distance(std::cbegin(alt_), itr) + 1; // + 1 for ref
                           } else {
                               return -1;
                           }
                       }
                   });
    genotype->phased = (phasing == Phasing::phased);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_genotype(const SampleName& sample, const std::vector<boost::optional<unsigned>>& alleles,
                                                     const Phasing phasing)
{
    auto& genotype = samples_[sample].genotype;
    genotype = VcfRecord::Genotype {};
    genotype->indices.resize(alleles.size());
    std::transform(std::cbegin(alleles), std::cend(alleles), std::begin(genotype->indices),
                   [] (const auto& allele) -> VcfRecord::AlleleIndex {
                       if (allele) {
                           return *allele; 
                       } else {
                           return -1;
                       }
                   });
    genotype->phased = (phasing == Phasing::phased);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::clear_genotype(const SampleName& sample) noexcept
{
    samples_.at(sample).genotype = boost::none;
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_format(const SampleName& sample, const KeyType& key, const ValueType& value)
{
    return this->set_format(sample, key, std::vector<ValueType> {value});
}

VcfRecord::Builder& VcfRecord::Builder::set_format(const SampleName& sample, const KeyType& key, std::vector<ValueType> values)
{
    samples_[sample].other[key] = std::move(values);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_format(const SampleName& sample, const KeyType& key, std::initializer_list<ValueType> values)
{
    return this->set_format(sample, key, std::vector<ValueType> {values});
}

VcfRecord::Builder& VcfRecord::Builder::set_format_missing(const SampleName& sample, const KeyType& key)
{
    return this->set_format(sample, key, std::string {vcfspec::missingValue});
}

VcfRecord::Builder& VcfRecord::Builder::clear_format() noexcept
{
    format_.clear();
    samples_.clear();
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::clear_format(const SampleName& sample) noexcept
{
    samples_.erase(sample);
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::clear_format(const SampleName& sample, const KeyType& key) noexcept
{
    if (key == vcfspec::format::genotype) {
        clear_genotype(sample);
    } else {
        const auto sample_itr = samples_.find(sample);
        if (sample_itr != std::cend(samples_)) {
            sample_itr->second.other.erase(key);
        }
    }
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_passed(const SampleName& sample)
{
    return this->set_format(sample, vcfspec::format::filter, {vcfspec::filter::pass});
}

VcfRecord::Builder& VcfRecord::Builder::set_filter(const SampleName& sample, std::vector<KeyType> filter)
{
    return this->set_format(sample, vcfspec::format::filter, std::move(filter));
}

VcfRecord::Builder& VcfRecord::Builder::set_filter(const SampleName& sample, std::initializer_list<KeyType> filter)
{
    return this->set_format(sample, vcfspec::format::filter, std::move(filter));
}

VcfRecord::Builder& VcfRecord::Builder::add_filter(const SampleName& sample, KeyType filter)
{
    auto& filters = samples_[sample].other[vcfspec::format::filter];
    // FORMAT/FT has Number=1 in the spec, as FILTER
    if (filters.empty()) {
        filters.push_back(std::move(filter));
    } else {
        filters.back() += vcfspec::filter::seperator + filter;
    }
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::clear_filter(const SampleName& sample) noexcept
{
    return this->clear_format(sample, vcfspec::format::filter);
}

VcfRecord::Builder& VcfRecord::Builder::clear_all_sample_filters() noexcept
{
    for (const auto& p : samples_) {
        this->clear_filter(p.first);
    }
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::set_somatic()
{
    return this->set_info_flag(vcfspec::info::somatic);
}

VcfRecord::Builder& VcfRecord::Builder::set_denovo()
{
    return this->set_info_flag("DENOVO");
}

VcfRecord::Builder& VcfRecord::Builder::set_reference_reversion()
{
    return this->set_info_flag("REVERSION");
}

VcfRecord::Builder& VcfRecord::Builder::set_blocked_reference()
{
    const static std::vector<VcfRecord::NucleotideSequence> refcall_alts {vcfspec::allele::nonref};
    if (alt_ != refcall_alts)
        throw std::runtime_error {"Cannot block a non reference call"};
    if (ref_.size() > 1) {
        set_info("END", pos_ + ref_.size() - 1);
        ref_.resize(1);
    }
    return *this;
}

VcfRecord::Builder& VcfRecord::Builder::squash_reference_if_blocked()
{
    if (info_.count("END") > 0) {
        ref_.resize(1);
    }
    return *this;
}

GenomicRegion::Position VcfRecord::Builder::pos() const noexcept
{
    return pos_;
}

void VcfRecord::Builder::collapse_spanning_deletions()
{
    const static std::vector<VcfRecord::NucleotideSequence> refcall_alts {vcfspec::allele::nonref};
    if (alt_ != refcall_alts) {
        for (auto& alt : alt_) {
            if (alt.size() > 1 && std::find(std::cbegin(alt), std::cend(alt), vcfspec::deleteMaskAllele[0]) != std::cend(alt)) {
                alt = vcfspec::deleteMaskAllele;
            }
        }
    }
}

VcfRecord VcfRecord::Builder::build() const
{
    if (format_.empty()) {
        if (end_) {
            GenomicRegion region {chrom_, pos_ - 1, *end_ - 1};
            return VcfRecord {std::move(region), id_, ref_, alt_, qual_, filter_, info_};
        } else {
            return VcfRecord {chrom_, pos_, id_, ref_, alt_, qual_, filter_, info_};
        }
    } else {
        if (end_) {
            GenomicRegion region {chrom_, pos_ - 1, *end_ - 1};
            return VcfRecord {std::move(region), id_, ref_, alt_, qual_, filter_,
                              info_, format_, samples_};
        } else {
            return VcfRecord {chrom_, pos_, id_, ref_, alt_, qual_, filter_,
                             info_, format_, samples_};
        }
    }
}

VcfRecord VcfRecord::Builder::build_once() noexcept
{
    if (format_.empty()) {
        if (end_) {
            GenomicRegion region {std::move(chrom_), pos_ - 1, *end_ - 1};
            if (reference_) {
                ref_ = reference_->fetch_sequence(region);
            }
            return VcfRecord {std::move(region), std::move(id_), std::move(ref_),
                              std::move(alt_), qual_, std::move(filter_), std::move(info_)};
        } else {
            return VcfRecord {std::move(chrom_), pos_, std::move(id_), std::move(ref_),
                              std::move(alt_), qual_, std::move(filter_), std::move(info_)};
        }
    } else {
        if (end_) {
            GenomicRegion region {std::move(chrom_), pos_ - 1, *end_ - 1};
            if (reference_) {
                ref_ = reference_->fetch_sequence(region);
            }
            return VcfRecord {std::move(region), std::move(id_), std::move(ref_),
                              std::move(alt_), qual_, std::move(filter_), std::move(info_),
                              std::move(format_), std::move(samples_)};
        } else {
            return VcfRecord {std::move(chrom_), pos_, std::move(id_), std::move(ref_),
                              std::move(alt_), qual_, std::move(filter_), std::move(info_),
                              std::move(format_), std::move(samples_)};
        }
    }
}

} // namespace octopus
