// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "vcf_utils.hpp"

#include <unordered_map>
#include <deque>
#include <queue>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <numeric>

#include "htslib/vcf.h"
#include "htslib/tbx.h"

#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "vcf_spec.hpp"

namespace octopus {

std::vector<std::string> get_contigs(const VcfHeader& header)
{
    using namespace vcfspec::header::meta;
    std::vector<std::string> result {};
    const auto& contigs_fields = header.structured_fields(tag::contig);
    result.reserve(contigs_fields.size());
    std::transform(std::cbegin(contigs_fields), std::cend(contigs_fields), std::back_inserter(result),
                   [] (const auto& field) { return field.at(struc::id); });
    return result;
}

unsigned get_field_cardinality(const VcfHeader::StructuredKey& key, const VcfRecord& record)
{
    return 0;
}

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                                           const VcfHeader::StructuredKey& key)
{
    return get_typed_info_values(header, key, record.info_value(key.value));
}

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                             const VcfRecord::SampleName sample,
                                             const VcfHeader::StructuredKey& key)
{
    return get_typed_format_values(header, key, record.get_sample_value(sample, key.value));
}

bool is_indexable(const boost::filesystem::path& vcf_path)
{
    const auto extension = vcf_path.extension().string();
    return extension == ".bcf" || extension == ".gz";
}

void index_vcf(const boost::filesystem::path& vcf_path)
{
    auto* const fp = hts_open(vcf_path.c_str(), "r");
    if (fp == nullptr) {
        throw std::ios::failure {vcf_path.c_str()};
    }
    const auto type = *hts_get_format(fp);
    hts_close(fp);
    if (type.format == bcf) {
        bcf_index_build(vcf_path.c_str(), 14);
    } else {
        tbx_index_build(vcf_path.c_str(), 0, &tbx_conf_vcf);
    }
}

void index_vcf(const VcfReader& reader)
{
    index_vcf(reader.path());
}

void index_vcfs(const std::vector<VcfReader>& readers)
{
    for (const auto& reader : readers) index_vcf(reader);
}

std::vector<VcfReader> writers_to_readers(std::vector<VcfWriter>&& writers, const bool keep_open)
{
    std::vector<VcfReader> result {};
    result.reserve(writers.size());
    for (auto&& writer : writers) {
        auto path = writer.path();
        writer.close();
        if (path) {
            result.emplace_back(std::move(*path));
            if (!keep_open) result.back().close();
        }
    }
    writers.clear();
    return result;
}

void copy(VcfReader& src, VcfWriter& dst)
{
    const bool is_closed {!src.is_open()};
    if (is_closed) src.open();
    if (!dst.is_header_written()) {
        dst << src.fetch_header();
    }
    constexpr std::size_t maxBufferSize {1000000};
    if (src.count_records() <= maxBufferSize) {
        dst << src.fetch_records();
    } else {
        auto p = src.iterate();
        std::copy(std::move(p.first), std::move(p.second), VcfWriterIterator {dst});
    }
    if (is_closed) src.close();
}

void copy(const VcfReader& src, VcfWriter& dst)
{
    if (!dst.is_header_written()) {
        dst << src.fetch_header();
    }
    constexpr std::size_t maxBufferSize {1000000};
    if (src.count_records() <= maxBufferSize) {
        dst << src.fetch_records();
    } else {
        auto p = src.iterate();
        std::copy(std::move(p.first), std::move(p.second), VcfWriterIterator {dst});
    }
}

void sort(const VcfReader& src, VcfWriter& dst)
{
    if (!dst.is_header_written()) {
        dst << src.fetch_header();
    }
    auto records = src.fetch_records();
    std::sort(std::begin(records), std::end(records));
    dst << records;
}

bool all_same_format(const std::vector<VcfHeader>& headers)
{
    return std::adjacent_find(std::cbegin(headers), std::cend(headers),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.file_format() != rhs.file_format();
                              }) == std::cend(headers);
}

bool contain_same_samples(const std::vector<VcfHeader>& headers)
{
    return std::adjacent_find(std::cbegin(headers), std::cend(headers),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.samples() != rhs.samples();
                              }) == std::cend(headers);
}

bool all_equal(const std::vector<VcfHeader>& headers)
{
    return std::adjacent_find(std::cbegin(headers), std::cend(headers), std::not_equal_to<VcfHeader>()) == std::cend(headers);
}

void add_contigs(const std::vector<VcfHeader>& headers, VcfHeader::Builder& hb)
{
    using namespace vcfspec::header::meta;
    std::unordered_set<std::string> added {};
    for (const auto& header : headers) {
        for (auto contig_fields : header.structured_fields(tag::contig)) {
            const auto& contig = contig_fields.at(struc::id);
            if (added.count(contig) == 0) {
                contig_fields.erase(contig);
                std::unordered_map<std::string, std::string> tmp {std::begin(contig_fields), std::end(contig_fields)};
                hb.add_contig(contig, std::move(tmp));
                added.insert(contig);
            }
        }
    }
}

VcfHeader merge(const std::vector<VcfHeader>& headers)
{
    if (headers.empty()) {
        return VcfHeader {};
    }
    if (!all_same_format(headers)) {
        throw std::logic_error {"cannot merge VcfHeader's with different formats"};
    }
    if (!contain_same_samples(headers)) {
        throw std::logic_error {"cannot merge VcfHeader's that do not contain the same samples"};
    }
    if (all_equal(headers)) {
        return headers.front();
    }
    VcfHeader::Builder hb {headers.front()};
    add_contigs(headers, hb);
    return hb.build_once();
}

VcfHeader fetch_header(VcfReader& reader)
{
    const bool is_closed {!reader.is_open()};
    if (is_closed) reader.open();
    const auto result = reader.fetch_header();
    if (is_closed) reader.close();
    return result;
}

std::vector<VcfHeader> get_headers(std::vector<VcfReader>& readers)
{
    std::vector<VcfHeader> result {};
    result.reserve(readers.size());
    std::transform(std::begin(readers), std::end(readers), std::back_inserter(result),
                   [] (VcfReader& reader) { return fetch_header(reader); });
    return result;
}

using VcfReaderRef               = std::reference_wrapper<VcfReader>;
using ContigRecordCountMap       = std::unordered_map<std::string, std::size_t>;
using ReaderContigRecordCountMap = std::unordered_map<VcfReaderRef, ContigRecordCountMap>;

ReaderContigRecordCountMap
get_contig_count_map(std::vector<VcfReader>& readers, const std::vector<std::string>& contigs)
{
    ReaderContigRecordCountMap result {};
    result.reserve(readers.size());
    for (auto& reader : readers) {
        const bool is_closed {!reader.is_open()};
        if (is_closed) reader.open();
        ContigRecordCountMap contig_counts {};
        contig_counts.reserve(contigs.size());
        for (const auto& contig : contigs) {
            contig_counts.emplace(contig, reader.count_records(contig));
        }
        if (is_closed) reader.close();
        result.emplace(reader, std::move(contig_counts));
    }
    return result;
}

auto count_active_contigs(const ContigRecordCountMap& counts)
{
    return std::count_if(std::cbegin(counts), std::cend(counts),
                         [] (const auto& p) { return p.second > 0; });
}

bool is_one_contig_per_reader(const ReaderContigRecordCountMap& counts)
{
    return std::all_of(std::cbegin(counts), std::cend(counts),
                       [] (const auto& p) { return count_active_contigs(p.second) == 1; });
}

auto find_active_contig(ContigRecordCountMap::const_iterator first, ContigRecordCountMap::const_iterator last)
{
    return std::find_if(first, last, [] (const auto& p) { return p.second > 0; });
}

bool is_unique_contig_per_reader(const ReaderContigRecordCountMap& counts)
{
    std::deque<std::reference_wrapper<const std::string>> contigs {};
    for (const auto& p : counts) {
        const auto it = find_active_contig(std::cbegin(p.second), std::cend(p.second));
        if (it != std::cend(p.second)) {
            if (find_active_contig(std::next(it), std::cend(p.second)) != std::cend(p.second)) {
                return false;
            }
            contigs.emplace_back(it->first);
        }
    }
    std::sort(std::begin(contigs), std::end(contigs),
              [] (const auto& lhs, const auto& rhs) { return lhs.get() < rhs.get(); });
    return std::adjacent_find(std::cbegin(contigs), std::cend(contigs),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.get() == rhs.get();
                              }) == std::cend(contigs);
}

std::size_t count_records(const ReaderContigRecordCountMap& counts)
{
    return std::accumulate(std::cbegin(counts), std::cend(counts), std::size_t {0},
                           [] (const auto curr, const auto& reader_counts) {
                               return curr + std::accumulate(std::cbegin(reader_counts.second),
                                                             std::cend(reader_counts.second),
                                                             std::size_t {0},
                                                             [] (const auto r_curr, const auto& p) {
                                                                 return r_curr + p.second;
                                                             }); });
}

auto extract_unique_readers(const ReaderContigRecordCountMap& reader_contig_counts)
{
    std::unordered_map<std::string, VcfReaderRef> result {};
    result.reserve(reader_contig_counts.size());
    for (const auto& p : reader_contig_counts) {
        const auto itr = find_active_contig(std::cbegin(p.second), std::cend(p.second));
        if (itr != std::cend(p.second)) {
            result.emplace(itr->first, p.first);
        }
    }
    return result;
}

void merge_contig_unique(const std::vector<VcfReader>& sources, VcfWriter& dst,
                         const std::vector<std::string>& contigs,
                         const ReaderContigRecordCountMap& reader_contig_counts)
{
    const auto contig_readers = extract_unique_readers(reader_contig_counts);
    for (const auto& contig : contigs) {
        if (contig_readers.count(contig) == 1) {
            copy(contig_readers.at(contig), dst);
        }
    }
}

using VcfRecordQueue = std::priority_queue<VcfRecord, std::deque<VcfRecord>, std::greater<VcfRecord>>;

void write(VcfRecordQueue& records, VcfWriter& dst)
{
    while (!records.empty()) {
        dst << records.top();
        records.pop();
    }
}

void one_step_merge(std::vector<VcfReader>& sources, VcfWriter& dst,
                    const std::vector<std::string>& contigs,
                    ReaderContigRecordCountMap& reader_contig_counts)
{
    VcfRecordQueue record_queue {};
    for (const auto& contig : contigs) {
        for (auto& reader : sources) {
            if (reader_contig_counts[reader].count(contig) == 1) {
                auto records = reader.fetch_records(contig, VcfReader::UnpackPolicy::all);
                for (auto&& record : records) {
                    record_queue.emplace(std::move(record));
                }
            }
        }
        write(record_queue, dst);
    }
}

void merge_pair(VcfReader& first, VcfReader& second, VcfWriter& dst,
                const std::vector<std::string>& contigs,
                ReaderContigRecordCountMap& reader_contig_counts)
{
    VcfWriterIterator out {dst};
    for (const auto& contig : contigs) {
        using std::move;
        if (reader_contig_counts[first].count(contig) == 1) {
            auto p1 = first.iterate(contig);
            if (reader_contig_counts[second].count(contig) == 1) {
                auto p2 = second.iterate(contig);
                std::merge(std::move(p1.first), std::move(p1.second), std::move(p2.first), std::move(p2.second), out);
            } else {
                std::copy(std::move(p1.first), std::move(p1.second), out);
            }
        } else if (reader_contig_counts[second].count(contig) == 1) {
            auto p2 = second.iterate(contig);
            std::copy(std::move(p2.first), std::move(p2.second), out);
        }
    }
}
    
bool is_empty(const VcfReader::RecordIteratorPair& p)
{
    return p.first == p.second;
}

const VcfRecord& front(const VcfReader::RecordIteratorPair& p)
{
    return *p.first;
}

struct RecordIteratorCompare
{
    using ItrPair = VcfReader::RecordIteratorPair;
    bool operator()(const ItrPair& lhs, const ItrPair& rhs) const
    {
        return front(lhs) > front(rhs);
    }
};

using RecordIteratorQueue = std::priority_queue<VcfReader::RecordIteratorPair,
                                                std::deque<VcfReader::RecordIteratorPair>,
                                                RecordIteratorCompare>;

auto make_record_iterator_queue(const std::vector<VcfReader>& sources, const std::string& contig)
{
    RecordIteratorQueue result {};
    for (const auto& reader : sources) {
        result.push(reader.iterate(contig));
    }
    return result;
}

void merge(std::vector<VcfReader>& sources, VcfWriter& dst, const std::vector<std::string>& contigs)
{
    if (sources.empty()) return;
    if (sources.size() == 1) {
        copy(sources.front(), dst);
        return;
    }
    if (!dst.is_header_written()) {
        dst << merge(get_headers(sources));
    }
    auto reader_contig_counts = get_contig_count_map(sources, contigs);
    if (is_unique_contig_per_reader(reader_contig_counts)) {
        merge_contig_unique(sources, dst, contigs, reader_contig_counts);
    } else {
        static constexpr std::size_t maxBufferSize {100000};
        if (count_records(reader_contig_counts) <= maxBufferSize) {
            one_step_merge(sources, dst, contigs, reader_contig_counts);
        } else if (sources.size() == 2) {
            merge_pair(sources.front(), sources.back(), dst, contigs, reader_contig_counts);
        } else {
            for (const auto& contig : contigs) {
                auto record_queue = make_record_iterator_queue(sources, contig);
                while (!record_queue.empty()) {
                    auto p = record_queue.top();
                    dst << front(p);
                    record_queue.pop();
                    ++p.first;
                    if (!is_empty(p)) {
                        record_queue.push(std::move(p));
                    }
                }
            }
        }
    }
}

void merge(std::vector<VcfReader>& sources, VcfWriter& dst)
{
    if (sources.empty()) return;
    if (sources.size() == 1) {
        copy(sources.front(), dst);
        return;
    }
    const auto header = merge(get_headers(sources));
    if (!dst.is_header_written()) {
        dst << header;
    }
    return merge(sources, dst, get_contigs(header));
}

namespace {

VcfHeader to_legacy(const VcfHeader& native)
{
    return VcfHeader::Builder(native).set_file_format("VCFv4.2").build_once();
}

bool is_deleted(const VcfRecord::NucleotideSequence& allele) noexcept
{
    return allele == vcfspec::deleteMaskAllele;
}

bool is_missing(const VcfRecord::NucleotideSequence& allele) noexcept
{
    static const std::string missing {vcfspec::missingValue};
    return allele == missing;
}

bool is_missing_or_deleted(const VcfRecord::NucleotideSequence& allele) noexcept
{
    return is_missing(allele) || is_deleted(allele);
}

VcfRecord convert_to_legacy(const VcfRecord& record, const std::vector<std::string>& samples)
{
    VcfRecord::Builder cb {record};
    const auto& alt = record.alt();
    const auto first_deleted_itr = std::find_if(std::cbegin(alt), std::cend(alt), is_deleted);
    if (first_deleted_itr != std::cend(alt)) {
        auto new_alt = alt;
        new_alt.erase(std::next(std::begin(new_alt), std::distance(std::cbegin(alt), first_deleted_itr)));
        cb.set_alt(std::move(new_alt));
    }
    for (const auto& sample : samples) {
        const auto& gt = record.get_sample_value(sample, vcfspec::format::genotype);
        const auto first_non_legacy = std::find_if(std::cbegin(gt), std::cend(gt), is_missing_or_deleted);
        if (first_non_legacy != std::cend(gt)) {
            const auto& ref = record.ref();
            auto new_gt = gt;
            auto ref_itr = std::next(std::begin(new_gt), std::distance(std::cbegin(gt), first_non_legacy));
            *ref_itr = ref;
            std::replace_if(std::next(ref_itr), std::end(new_gt), is_missing_or_deleted, ref);
            auto phasing = VcfRecord::Builder::Phasing::phased;
            if (!record.is_sample_phased(sample)) {
                phasing = VcfRecord::Builder::Phasing::unphased;
            }
            cb.set_genotype(sample, std::move(new_gt), phasing);
        }
    }
    return cb.build_once();
}

void convert_to_legacy_keep_all(const VcfReader& src, VcfWriter& dst)
{
    const auto samples = src.fetch_header().samples();
    for (auto p = src.iterate(); p.first != p.second; ++p.first) {
        dst << convert_to_legacy(*p.first, samples);
    }
}

bool is_snv(const VcfRecord::NucleotideSequence& ref, const VcfRecord::NucleotideSequence& alt) noexcept
{
    return ref.size() == 1 && alt.size() == 1 && ref != alt;
}

bool only_snvs(const VcfRecord& record) noexcept
{
    return std::all_of(std::cbegin(record.alt()), std::cend(record.alt()),
                       [&] (const auto& allele) { return is_snv(record.ref(), allele); });
}

bool is_indel(const VcfRecord::NucleotideSequence& ref, const VcfRecord::NucleotideSequence& alt) noexcept
{
    return ref.size() != alt.size();
}

bool only_indels(const VcfRecord& record) noexcept
{
    return std::any_of(std::cbegin(record.alt()), std::cend(record.alt()),
                       [&] (const auto& allele) { return is_indel(record.ref(), allele); });
}

bool contains(const std::vector<VcfRecord::NucleotideSequence>& genotype, const VcfRecord::NucleotideSequence& allele) noexcept
{
    return std::find(std::cbegin(genotype), std::cend(genotype), allele) != std::cend(genotype);
}

bool both_alts_genotyped(const VcfRecord& lhs, const VcfRecord::NucleotideSequence& lhs_alt,
                         const VcfRecord& rhs, const VcfRecord::NucleotideSequence& rhs_alt,
                         const std::string& sample) noexcept
{
    const auto& lhs_genotype = get_genotype(lhs, sample);
    if (contains(lhs_genotype, lhs_alt)) {
        const auto& rhs_genotype = get_genotype(rhs, sample);
        return contains(rhs_genotype, rhs_alt);
    }
    return false;
}

bool genotyped_in_same_samples(const VcfRecord& lhs, const VcfRecord::NucleotideSequence& lhs_alt,
                               const VcfRecord& rhs, const VcfRecord::NucleotideSequence& rhs_alt,
                               const std::vector<std::string>& samples) noexcept
{
    return std::all_of(std::cbegin(samples), std::cend(samples),
                       [&] (const auto& sample) { return both_alts_genotyped(lhs, lhs_alt, rhs, rhs_alt, sample); });
}

bool is_lhs_ref_flank_snv_duplicate(const VcfRecord& lhs, const VcfRecord& rhs,
                                    const std::vector<std::string>& samples) noexcept
{
    // Positive examples:
    //
    // N	1000	.	T	G	.	.	.	GT	1|0
    // N	1000	.	TCCTTTC	G	.	.	.	GT	1|0
    //
    // N    1000   .   C    G   .   .   .   GT  1|0
    // N    1000    .   C   GT,CTT  .   .   .   GT  1|2
    //
    // N    1000   .   C    G   .   .   .   GT  1|0 0|0
    // N    1000    .   C   GT,CTT  .   .   .   GT  1|2 0|2
    //
    // Negative examples:
    //
    // N	1000	.	T	G	.	.	.	GT	1|0
    // N	1000	.	TCCTTTC	T	.	.	.	GT	0|1
    //
    // N    1000   .   C    G   .   .   .   GT  1|0 0|1
    // N    1000    .   C   GT,CTT  .   .   .   GT  1|2 1|0
    //
    if (!is_same_contig(lhs, rhs) || !begins_equal(lhs, rhs)) return false;
    if (!only_snvs(lhs) || !only_indels(rhs)) return false;
    bool result {false};
    for (const auto& lhs_alt : lhs.alt()) {
        for (const auto& rhs_alt : rhs.alt()) {
            if (lhs_alt[0] == rhs_alt.front()) {
                if (!genotyped_in_same_samples(lhs, lhs_alt, rhs, rhs_alt, samples)) {
                    return false;
                } else {
                    result = true;
                }
            }
        }
    }
    return result;
}

void convert_to_legacy_dedup(const VcfReader& src, VcfWriter& dst)
{
    auto p = src.iterate();
    if (p.first == p.second) return;
    const auto samples = src.fetch_header().samples();
    VcfRecord prev_record {*p.first};
    ++p.first;
    for (; p.first != p.second; ++p.first) {
        const auto& record = *p.first;
        if (record != prev_record && !is_lhs_ref_flank_snv_duplicate(prev_record, record, samples)) {
            dst << convert_to_legacy(prev_record, samples);
        }
        prev_record = record;
    }
    dst << convert_to_legacy(prev_record, samples);
}

} // namespace

void convert_to_legacy(const VcfReader& src, VcfWriter& dst, const bool remove_ref_pad_duplicates)
{
    if (!dst.is_header_written()) {
        dst << to_legacy(src.fetch_header());
    }
    if (remove_ref_pad_duplicates) {
        convert_to_legacy_dedup(src, dst);
    } else {
        convert_to_legacy_keep_all(src, dst);
    }
}

} // namespace octopus
