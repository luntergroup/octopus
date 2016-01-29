/*  bioio.hpp -- FASTA/Q I/O
 
 Copyright (C) 2015 University of Oxford.
 
 Author: Daniel Cooke <dcooke@well.ox.ac.uk>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.  */

#ifndef __bioio__bioio__
#define __bioio__bioio__

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <limits>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <utility>
#include <type_traits>

namespace bioio
{   
    namespace detail
    {
        inline std::vector<std::string> split(const std::string& str, const char delim)
        {
            std::stringstream ss {str};
            std::string item;
            std::vector<std::string> result {};
            while (std::getline(ss, item, delim)) result.emplace_back(item);
            return result;
        }
    } // namespace detail
    
   /*=======================================================================================
    Types
    =======================================================================================*/
    
    struct FastaContigIndex
    {
        std::string contig_name;
        size_t length;
        size_t offset;
        size_t line_length;
        size_t line_byte_length;
        
        FastaContigIndex() = default;
        
        template <typename T>
        explicit FastaContigIndex(T&& contig_name, size_t length, size_t offset,
                                  size_t line_length, size_t line_byte_length)
        :
        contig_name {std::forward<T>(contig_name)},
        offset {offset}, 
        length {length},
        line_length {line_length}, 
        line_byte_length {line_byte_length} 
        {}
        
        template <typename T>
        explicit FastaContigIndex(const T& fasta_index_line)
        {
            const auto parts = detail::split(fasta_index_line, '\t');
            
            contig_name      = parts[0];
            length           = std::stoull(parts[1]);
            offset           = std::stoull(parts[2]);
            line_length      = std::stoull(parts[3]);
            line_byte_length = std::stoull(parts[4]);
        }
    };
    
    using FastaIndex = std::unordered_map<std::string, FastaContigIndex>;
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    struct FastaRecord
    {
        StringType name;
        SequenceType sequence;
        
        FastaRecord() = delete;
        template<typename StringType_, typename SequenceType_>
        explicit FastaRecord(StringType_&& name, SequenceType_&& sequence)
        : 
        name {std::forward<StringType_>(name)}, 
        sequence {std::forward<SequenceType_>(sequence)}
        {}
    };
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string>
    struct FastqRecord
    {
        StringType name;
        SequenceType1 seq;
        SequenceType2 qual;
        
        FastqRecord() = delete;
        template<typename StringType_, typename SequenceType1_, typename SequenceType2_>
        explicit FastqRecord(StringType_&& name, SequenceType1_&& seq, SequenceType2_&& qual)
        :
        name {std::forward<StringType_>(name)}, 
        seq {std::forward<SequenceType1_>(seq)}, 
        qual {std::forward<SequenceType2_>(qual)} 
        {}
    };
    
    template <typename StringType = std::string>
    using ReadIds = std::vector<StringType>;
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    using FastaMap = std::unordered_map<StringType, SequenceType>;
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string>
    using FastqMap = std::unordered_map<StringType, std::pair<SequenceType1, SequenceType2>>;
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    using FastaReads = std::pair<ReadIds<StringType>, FastaMap<StringType, SequenceType>>;
    
    template <typename StringType = std::string, typename SequenceType = std::string,
              typename SequenceType2 = std::string>
    using FastqReads = std::pair<ReadIds<StringType>, FastqMap<StringType, SequenceType, SequenceType2>>;
    
    namespace detail
    {
        // Counts the occurrences of record_delim that follow a newline '\n'
        inline size_t count_records(std::istream& is, const char record_delim)
        {
            const auto current_position = is.tellg();
            
            size_t result {};
            
            while (is) {
                if (is.peek() == record_delim) ++result;
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            
            is.clear();
            is.seekg(current_position, std::ios::beg);
            
            return result;
        }
        
        static constexpr char fasta_delim {'>'};
        static constexpr char fastq_delim {'@'};
        
        template<typename StringType = std::string, typename SequenceType = std::string>
        ::bioio::FastaRecord<StringType, SequenceType>
        read_fasta_record(std::istream& fasta)
        {
            StringType name;
            std::getline(fasta, name); // name is always a single line
            
            SequenceType line;
            std::getline(fasta, line);
            
            // The FASTA format is not as simple as FASTQ - the sequence
            // may be broken into multiple lines. We assume each line is the same size.
            if (!fasta.good() || fasta.peek() == fasta_delim) {
                return ::bioio::FastaRecord<StringType, SequenceType> {std::move(name), std::move(line)};
            } else {
                const auto line_size = line.size();
                
                line.resize(line_size);
                
                SequenceType seq {};
                seq.reserve(line_size * 100);
                
                const auto line_begin = line.cbegin(), line_end = line.cend();
                
                seq.insert(seq.end(), line_begin, line_end);
                
                while (fasta.good()) {
                    fasta.getline(&line[0], line_size + 1);
                    
                    if (line.front() == fasta_delim) {
                        if (fasta.good()) fasta.seekg(-fasta.gcount(), std::ios_base::cur);
                        break;
                    }
                    
                    if (fasta.gcount() < line_size) {
                        seq.insert(seq.end(), line_begin, line_begin + fasta.gcount() - 1);
                        break;
                    }
                    
                    seq.insert(seq.end(), line_begin, line_end);
                }
                
                seq.shrink_to_fit();
                
                return ::bioio::FastaRecord<StringType, SequenceType> {std::move(name), std::move(seq)};
            }
        }
        
        template<typename StringType = std::string, typename SequenceType1 = std::string,
                 typename SequenceType2 = std::string>
        ::bioio::FastqRecord<StringType, SequenceType1, SequenceType2>
        read_fastq_record(std::istream& fastq)
        {
            StringType name;
            SequenceType1 seq;
            SequenceType2 quals;
            
            // Unlike FASTA, FASTQ always use one line per field.
            std::getline(fastq, name);
            std::getline(fastq, seq);
            fastq.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::getline(fastq, quals);
            
            return ::bioio::FastqRecord<StringType, SequenceType1, SequenceType2> 
                    {std::move(name), std::move(seq), std::move(quals)};
        }
    } // namespace detail
    
    /*=======================================================================================
     INDEX: For reading FASTA index files
     =======================================================================================*/
    
    inline size_t count_contigs_in_fasta_index(std::istream& fasta_index)
    {
        return std::count(std::istreambuf_iterator<char>(fasta_index), 
                          std::istreambuf_iterator<char>(), '\n');
    }
    
    inline size_t count_contigs_in_fasta_index(const std::string& fasta_index_path)
    {
        std::ifstream fasta_index {fasta_index_path, std::ios::binary};
        return count_contigs_in_fasta_index(fasta_index);
    }
    
    inline std::vector<std::string> read_fasta_index_contig_names(std::istream& fasta_index)
    {
        std::vector<std::string> result {};
        result.reserve(100);
        
        std::string line;
        
        while (std::getline(fasta_index, line)) {
            result.emplace_back(line.substr(0, line.find('\t')));
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    inline std::vector<std::string> read_fasta_index_contig_names(const std::string& fasta_index_path)
    {
        std::ifstream fasta_index {fasta_index_path, std::ios::binary};
        return read_fasta_index_contig_names(fasta_index);
    }
    
    inline FastaIndex read_fasta_index(std::istream& fasta_index)
    {
        FastaIndex result {};
        result.reserve(100);
        
        std::string line;
        
        while (std::getline(fasta_index, line)) {
            FastaContigIndex contig_index {line};
            result.emplace(contig_index.contig_name, std::move(contig_index));
        }
        
        result.rehash(result.size());
        
        return result;
    }
    
    inline FastaIndex read_fasta_index(const std::string& fasta_index_path)
    {
        std::ifstream fasta_index {fasta_index_path, std::ios::binary};
        return read_fasta_index(fasta_index);
    }
    
    inline size_t calculate_contig_size(std::istream& fasta_index, const std::string& contig_name)
    {
        const auto current_position = fasta_index.tellg();
        
        fasta_index.seekg(0, std::ios::beg);
        
        std::string line;
        
        while (std::getline(fasta_index, line)) {
            FastaContigIndex contig_index {line};
            if (contig_index.contig_name == contig_name) return contig_index.length;
        }
        
        fasta_index.clear();
        fasta_index.seekg(current_position, std::ios::beg);
        
        return 0;
    }
    
    inline size_t calculate_contig_size(const std::string& fasta_index_path, const std::string& contig_name)
    {
        std::ifstream fasta_index {fasta_index_path, std::ios::binary};
        return calculate_contig_size(fasta_index, contig_name);
    }
    
    /*=======================================================================================
     FASTA: For reading FASTAs with index files
     =======================================================================================*/
    
    namespace detail
    {
        inline size_t line_offset(const ::bioio::FastaContigIndex& index, const size_t begin)
        {
            return begin % index.line_length;
        }
        
        inline size_t region_offset(const ::bioio::FastaContigIndex& index, const size_t begin)
        {
            return index.offset + begin / index.line_length * index.line_byte_length + line_offset(index, begin);
        }
        
        inline size_t remaining_line_length(const ::bioio::FastaContigIndex& index, const size_t begin)
        {
            return index.line_length - line_offset(index, begin);
        }
    } // namespace detail
    
    template <typename SequenceType = std::string>
    SequenceType 
    read_fasta_contig(std::istream& fasta, const FastaContigIndex& index,
                      const size_t begin, size_t length)
    {
        SequenceType result {};
        
        if (length == 0 || begin >= index.length) return result;
        
        fasta.seekg(detail::region_offset(index, begin), std::ios::beg);
        
        length = std::min(length, index.length - begin);
        
        if (index.line_length == index.line_byte_length) {
            result.resize(length);
            fasta.read(&result[0], length);
        } else {
            const auto num_line_end_bytes = index.line_byte_length - index.line_length;
            
            const auto num_remaining_curr_line_bytes = detail::remaining_line_length(index, begin);
            
            if (length <= num_remaining_curr_line_bytes) {
                result.resize(length);
                fasta.read(&result[0], length);
            } else {
                // Allocate enough space to fit a full last line so we don't need to keep
                // checking how much of the final line to read.
                result.resize(length + detail::remaining_line_length(index, begin + length) + num_line_end_bytes);
                
                fasta.read(&result[0], num_remaining_curr_line_bytes + num_line_end_bytes);
                
                for (auto pos = num_remaining_curr_line_bytes; pos < length; pos += index.line_length) {
                    fasta.read(&result[pos], index.line_byte_length);
                }
                
                result.resize(length);
            }
        }
        
        fasta.clear(); // assumes indexed queries do not need eof flag
        
        return result;
    }
    
    template <typename SequenceType = std::string>
    SequenceType
    read_fasta_contig(const std::string& fasta_path, const FastaContigIndex& index,
                      const size_t begin, const size_t length)
    {
        std::ifstream fasta {fasta_path, std::ios::binary | std::ios::beg};
        return read_fasta_contig<SequenceType>(fasta, index, begin, length);
    }
    
    template <typename SequenceType = std::string>
    SequenceType read_fasta_contig(std::istream& fasta, const FastaContigIndex& index)
    {
        return read_fasta_contig<SequenceType>(fasta, index, 0, index.length);
    }
    
    template <typename SequenceType = std::string>
    SequenceType read_fasta_contig(const std::string& fasta_path, const FastaContigIndex& index)
    {
        std::ifstream fasta {fasta_path, std::ios::binary | std::ios::beg};
        return read_fasta_contig<SequenceType>(fasta, index);
    }
    
    template <typename T, typename U>
    std::ostream& operator<<(std::ostream& os, const FastaRecord<T, U>& data)
    {
        os << ">" << data.first << "\n" << data.second;
        return os;
    }
    
   /*=======================================================================================
    FASTA: Optimised for reading a Fasta file with a single contig without an index
    =======================================================================================*/
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    FastaRecord<StringType, SequenceType> read_single_contig_fasta(const std::string& fasta_path)
    {
        std::ifstream fasta {fasta_path, std::ios::binary | std::ios::ate};
        
        const auto file_length = static_cast<size_t>(fasta.tellg());
        
        fasta.seekg(0, std::ios::beg);
        
        StringType name;
        std::getline(fasta, name);
        
        SequenceType sequence;
        sequence.resize(file_length - name.size() - 1);
        
        fasta.read(&sequence[0], sequence.size());
        
        sequence.erase(std::remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
        
        return FastaRecord<StringType, SequenceType> {std::move(name), std::move(sequence)};
    }
    
    /*=======================================================================================
     FASTA: For reading multiple record Fasta files without an index
     =======================================================================================*/
    
    inline size_t count_fasta_records(std::istream& fasta)
    {
        return detail::count_records(fasta, detail::fasta_delim);
    }
    
    inline size_t count_fasta_records(const std::string& fasta_path)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return count_fasta_records(fasta);
    }
    
    template <typename UnaryPredicate,
              typename = typename std::enable_if<
                                      !std::is_convertible<UnaryPredicate, std::string>::value
                                  >::type>
    bool seek_fasta_record(std::istream& fasta, UnaryPredicate pred)
    {
        std::string record;
        while (std::getline(fasta, record)) {
            if (pred(record)) {
                fasta.seekg(-(record.size() + 1), std::ios_base::cur);
                return true;
            }
            while (fasta) {
                if (fasta.peek() == detail::fasta_delim) break;
                fasta.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
        return false;
    }
    
    inline bool seek_fasta_record(std::istream& fasta, const std::string& name)
    {
        return seek_fasta_record(fasta,
                                 [&] (const std::string& fasta_record) {
                                     if (fasta_record.size() < name.size()) return false;
                                     return std::equal(name.cbegin(), name.cend(),
                                                       std::next(fasta_record.cbegin()));
                                 });
    }
    
    template <typename SequenceType = std::string, typename UnaryPredicate>
    SequenceType read_fasta_seq(std::istream& fasta, UnaryPredicate pred)
    {
        if (seek_fasta_record(fasta, pred)) {
            return detail::read_fasta_record<std::string, SequenceType>(fasta).sequence;
        }
        return SequenceType {};
    }
    
    template <typename SequenceType = std::string, typename UnaryPredicate>
    SequenceType read_fasta_seq(const std::string& fasta_path, UnaryPredicate pred)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta_seq(fasta, pred);
    }
    
    template <typename SequenceType = std::string>
    SequenceType read_fasta_seq(std::istream& fasta, const std::string& name)
    {
        if (seek_fasta_record(fasta, name)) {
            return detail::read_fasta_record<std::string, SequenceType>(fasta).sequence;
        }
        return SequenceType {};
    }
    
    template <typename SequenceType = std::string>
    SequenceType read_fasta_seq(const std::string& fasta_path, const std::string& name)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta_seq(fasta, name);
    }
    
    template <typename SequenceType = std::string>
    std::vector<SequenceType> read_fasta_seqs(std::istream& fasta, size_t num_records)
    {
        std::vector<SequenceType> result {};
        result.reserve(num_records);
        
        for (; num_records > 0; --num_records) {
            result.emplace_back(detail::read_fasta_record<std::string, SequenceType>(fasta).sequence);
        }
        
        return result;
    }
    
    template <typename SequenceType = std::string>
    std::vector<SequenceType> read_fasta_seqs(std::istream& fasta)
    {
        return read_fasta_seqs<SequenceType>(fasta, count_fasta_records(fasta));
    }
    
    template <typename SequenceType = std::string>
    std::vector<SequenceType> read_fasta_seqs(const std::string& fasta_path)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta_seqs<SequenceType>(fasta);
    }
    
    template <typename SequenceType = std::string>
    std::vector<SequenceType> read_fasta_seqs(const std::string& fasta_path,
                                              const size_t num_records)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta_seqs<SequenceType>(fasta, num_records);
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    std::vector<FastaRecord<StringType, SequenceType>> read_fasta(std::istream& fasta,
                                                                  size_t num_records)
    {
        std::vector<FastaRecord<StringType, SequenceType>> result {};
        result.reserve(num_records);
        
        for (; num_records > 0; --num_records) {
            result.emplace_back(detail::read_fasta_record<StringType, SequenceType>(fasta));
        }
        
        return result;
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    std::vector<FastaRecord<StringType, SequenceType>> read_fasta(std::istream& fasta)
    {
        return read_fasta<StringType, SequenceType>(fasta, count_fasta_records(fasta));
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    std::vector<FastaRecord<StringType, SequenceType>> read_fasta(const std::string& fasta_path)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta<StringType, SequenceType>(fasta);
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    std::vector<FastaRecord<StringType, SequenceType>> read_fasta(const std::string& fasta_path,
                                                                  const size_t num_records)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta<StringType, SequenceType>(fasta, num_records);
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string,
              typename UnaryOperation>
    FastaReads<StringType, SequenceType> read_fasta_map(std::istream& fasta,
                                                        size_t num_records,
                                                        UnaryOperation unary_op)
    {
        FastaMap<StringType, SequenceType> records {};
        records.reserve(num_records);
        ReadIds<StringType> contig_names {};
        contig_names.reserve(num_records);
        
        for (; num_records > 0; --num_records) {
            auto record = detail::read_fasta_record<StringType, SequenceType>(fasta);
            auto f_name = unary_op(std::move(record.name));
            records.emplace(f_name, std::move(record.sequence));
            contig_names.insert(std::move(f_name));
        }
        
        return FastaReads<StringType, SequenceType> {std::move(contig_names), std::move(records)};
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string,
              typename UnaryOperation>
    FastaReads<StringType, SequenceType> read_fasta_map(std::istream& fasta,
                                                        UnaryOperation unary_op)
    {
        return read_fasta_map<StringType, SequenceType, UnaryOperation>(fasta,
                                                                        count_fasta_records(fasta),
                                                                        unary_op);
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string,
              typename UnaryOperation>
    FastaReads<StringType, SequenceType> read_fasta_map(const std::string& fasta_path,
                                                        UnaryOperation unary_op)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta_map<StringType, SequenceType, UnaryOperation>(fasta, unary_op);
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string,
              typename UnaryOperation>
    FastaReads<StringType, SequenceType> read_fasta_map(const std::string& fasta_path,
                                                        size_t num_records,
                                                        UnaryOperation unary_op)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta_map<StringType, SequenceType, UnaryOperation>(fasta, num_records, unary_op);
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    FastaReads<StringType, SequenceType> read_fasta_map(const std::string& fasta_path)
    {
        return read_fasta_map<StringType, SequenceType>(fasta_path,
                                                        [] (StringType&& name) {
                                                            return name;
                                                        });
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string,
              typename UnaryOperation>
    FastaReads<StringType, SequenceType>
    read_fasta_map(std::istream& fasta, const ReadIds<StringType>& names,
                   UnaryOperation unary_op)
    {
        auto num_records = names.size();
        
        FastaMap<StringType, SequenceType> records {};
        records.reserve(num_records);
        ReadIds<StringType> f_names {};
        f_names.reserve(num_records);
        
        for (; num_records > 0; --num_records) {
            auto record = detail::read_fasta_record<StringType, SequenceType>(fasta);
            auto f_name = unary_op(std::move(record.name));
            if (names.find(f_name) != names.end()) {
                records.emplace(f_name, std::move(record.sequence));
                f_names.insert(std::move(f_name));
            }
        }
        
        return FastaReads<StringType, SequenceType> {std::move(f_names), std::move(records)};
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string,
              typename UnaryOperation>
    FastaReads<StringType, SequenceType>
    read_fasta_map(const std::string& fasta_path, const ReadIds<StringType>& names,
                   UnaryOperation unary_op)
    {
        std::ifstream fasta {fasta_path, std::ios::binary};
        return read_fasta_map<StringType, SequenceType, UnaryOperation>(fasta, names, unary_op);
    }
    
    template <typename StringType = std::string, typename SequenceType = std::string>
    FastaReads<StringType, SequenceType>
    read_fasta_map(const std::string& path, const ReadIds<StringType>& names)
    {
        return read_fasta_map(path, names, 
                              [] (StringType&& name) { return name; });
    }
    
    template <typename T, typename U>
    std::ostream& operator<<(std::ostream& os, const FastaReads<T, U>& records)
    {
        for (auto name : records.first) {
            os << ">" << name << "\n" << records.second.at(name);
        }
        return os;
    }
    
    template <typename T, typename U>
    void write_fasta(const std::string& path, const FastaReads<T, U>& records)
    {
        std::ofstream fasta {path, std::ios::out | std::ios::binary};
        fasta << records;
    }
    
    /*=======================================================================================
     FASTQ: For reading multiple line Fastq files.
     =======================================================================================*/
    
    inline size_t count_fastq_records(std::istream& fastq)
    {
        return detail::count_records(fastq, detail::fastq_delim);
    }
    
    inline size_t count_fastq_records(const std::string& fastq_path)
    {
        std::ifstream fastq {fastq_path, std::ios::binary};
        return count_fasta_records(fastq);
    }
    
    template <typename SequenceType = std::string>
    std::vector<SequenceType>
    read_fastq_seqs(std::istream& fastq)
    {
        auto num_records = count_fastq_records(fastq);
        
        std::vector<SequenceType> result {};
        result.reserve(num_records);
        
        for (; num_records > 0; --num_records) {
            result.push_back(detail::read_fastq_record<std::string, SequenceType, std::string>(fastq).seq);
        }
        
        return result;
    }
    
    template <typename SequenceType = std::string>
    std::vector<SequenceType>
    read_fastq_seqs(const std::string& path)
    {
        std::ifstream fastq {path, std::ios::binary};
        return read_fastq_seqs(fastq);
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string, typename F>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(std::istream& fastq, size_t num_records, F f)
    {
        std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>> result {};
        result.reserve(num_records);
        
        for (; num_records > 0; --num_records) {
            auto record = detail::read_fastq_record<StringType, SequenceType1, SequenceType2>(fastq);
            result.emplace_back(f(std::move(record.name)), std::move(record.seq), std::move(record.qual));
        }
        
        return result;
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string, typename F>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(std::istream& fastq, F f)
    {
        return read_fastq<StringType, SequenceType1, SequenceType2, F>(fastq, count_fastq_records(fastq), f);
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string, typename F>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(const std::string& path, const size_t num_records, F f)
    {
        std::ifstream fastq {path, std::ios::binary};
        return read_fastq<StringType, SequenceType1, SequenceType2, F>(fastq, num_records, f);
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string, typename F>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(const std::string& path, F f)
    {
        std::ifstream fastq {path, std::ios::binary};
        return read_fastq<StringType, SequenceType1, SequenceType2, F>(fastq, f);
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(const std::string& path)
    {
        return read_fastq<StringType, SequenceType1, SequenceType2>(path, [] (const StringType& name) { return name; });
    }    
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(std::istream& fastq)
    {
        return read_fastq<StringType, SequenceType1, SequenceType2>(fastq, [] (const StringType& name) { return name; });
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(std::istream& fastq, const size_t num_records)
    {
        return read_fastq<StringType, SequenceType1, SequenceType2>(fastq, num_records, [] (const StringType& name) { 
            return name; 
        });
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string, typename IntegerType>
    std::vector<FastqRecord<StringType, SequenceType1, SequenceType2>>
    read_fastq(const std::string& path, const size_t num_records)
    {
        std::ifstream fastq {path, std::ios::binary};
        return read_fastq<StringType, SequenceType1, SequenceType2>(fastq, num_records, [] (const StringType& name) { 
            return name; 
        });
    }
    
    template <typename StringType = std::string, typename SequenceType1 = std::string,
              typename SequenceType2 = std::string>
    FastqReads<StringType, SequenceType1, SequenceType2>
    read_fastq_map(const std::string& path)
    {
        std::ifstream fastq {path, std::ios::binary};
        
        auto num_records = count_fastq_records(fastq);
        
        ReadIds<StringType> names {};
        FastqMap<StringType, SequenceType1, SequenceType2> data {};
        data.reserve(num_records);
        
        for (; num_records > 0; --num_records) {
            auto record = detail::read_fastq_record<StringType, SequenceType1, SequenceType2>(fastq);
            data.emplace(record.name, {std::move(record.seq), std::move(record.qual)});
            names.insert(std::move(record.name));
        }
        
        return {names, data};
    }
    
} // namespace bioio

#endif /* defined(__bioio__bioio__) */
