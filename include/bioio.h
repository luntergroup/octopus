//
//  bioio.h
//
//  Created by Daniel Cooke on 13/10/2014.
//  Copyright (c) 2014 Oxford University. All rights reserved.
//

#ifndef __bioio__bioio__
#define __bioio__bioio__

#include <string>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <sstream>
#include <exception>
#include <algorithm>

namespace bioio {

using std::size_t;

/*=======================================================================================
 MISC: functions for getting file information.
 =======================================================================================*/

inline
size_t get_num_records(std::ifstream& file, char record_delimiter)
{
    size_t num_records = 0;
    std::string first_word;
    while (file >> first_word) {
        if (first_word[0] == record_delimiter) {
            ++num_records;
        }
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // Uses pre-existing filestream so, need to reset eof flags and rewind.
    file.clear();
    file.seekg(0, std::ios::beg);
    return num_records;
}

template <typename T>
std::vector<std::string> split(T&& s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.emplace_back(item);
    }
    return elems;
}

/*=======================================================================================
 TYPEDEFS
 =======================================================================================*/

struct FastaIndex
{
    std::string contig_name;
    size_t length;
    size_t offset;
    size_t line_length;
    size_t line_byte_length;
    
    FastaIndex() = default;
    FastaIndex(std::string contig_name, size_t length, size_t offset, 
              size_t line_length, size_t line_byte_length) 
    : contig_name {contig_name}, offset {offset}, length {length}, 
      line_length {line_length}, 
      line_byte_length {line_byte_length} {}
    FastaIndex(std::string fasta_index_line)
    {
        auto parts = split(std::move(fasta_index_line), '\t');
        
        contig_name      = parts[0];
        length           = std::stoll(parts[1]);
        offset           = std::stoll(parts[2]);
        line_length      = std::stoll(parts[4]);
        line_byte_length = std::stoll(parts[3]);
    }
};

struct FastaIndexHash
{
    size_t operator()(const FastaIndex& a) const
    {
        return std::hash<std::string>()(a.contig_name);
    }
};

template <typename T1=std::string, typename T2=std::string>
struct FastaRecord
{
    T1 contig_name;
    T2 seq;
    
    FastaRecord() = delete;
    template<typename T1_, typename T2_>
    FastaRecord(T1_&& contig_name, T2_&& seq) 
        : contig_name {std::forward<T1_>(contig_name)}, seq {std::forward<T2_>(seq)} {}
};

template <typename T1=std::string, typename T2=std::string, typename T3=std::string>
struct FastqRecord
{
    T2 seq;
    T3 qual;
    T1 contig_name;
    
    FastqRecord() = delete;
    template<typename T1_, typename T2_, typename T3_>
    FastqRecord(T1_&& contig_name, T2_&& seq, T3_&& qual) 
        : contig_name {std::forward<T1_>(contig_name)}, seq {std::forward<T2_>(seq)}
                , qual {std::forward<T3_>(qual)} {}
};

template <typename T=std::string>
using ReadIds = std::unordered_set<T>;

template <typename T1=std::string, typename T2=std::string>
using FastaMap = std::unordered_map<T1, T2>;
template <typename T1=std::string, typename T2=std::string, typename T3=std::string>
using FastqMap = std::unordered_map<T1, std::pair<T2, T3>>;

template <typename T1=std::string, typename T2=std::string>
using FastaReads = std::pair<ReadIds<T1>, FastaMap<T1, T2>>;
template <typename T1=std::string, typename T2=std::string, typename T3=std::string>
using FastqReads = std::pair<ReadIds<T1>, FastqMap<T1, T2, T3>>;

class BadIndexException : public std::exception
{
    virtual const char* what() const throw()
    {
        return "Invlacontig_name fasta index.";
    }
};

/*=======================================================================================
    REFERENCE: Optimised for reading a Fasta file with a single contig
 =======================================================================================*/

template <typename T=std::string>
T read_single_contig_fasta_seq(const std::string& fasta_path)
{
    std::ifstream fasta(fasta_path, std::ios::binary | std::ios::ate);
    std::string contig_name;
    T seq;
    size_t file_len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    std::getline(fasta, contig_name, '\n');
    seq.resize(file_len - contig_name.size());
    fasta.read(&seq[0], seq.size());
    auto new_end = std::remove(seq.begin(), seq.end(), '\n');
    seq.resize(new_end - seq.begin() - 1);
    return seq;
}

template <typename T1=std::string, typename T2=std::string>
FastaRecord<T1, T2>
read_ref(const std::string& fasta_path)
{
    std::ifstream fasta(fasta_path, std::ios::binary | std::ios::ate);
    T1 contig_name;
    T2 seq;
    size_t len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    std::getline(fasta, contig_name, '\n');
    seq.resize(len - contig_name.size() - 1);
    fasta.read(&seq[0], seq.size());
    std::remove(seq.begin(), seq.end(), '\n');
    return FastaRecord<T1, T2>(contig_name, seq);
}

/*=======================================================================================
    INDEX: For reading FASTA index files
 =======================================================================================*/

inline
std::vector<std::string>
get_fasta_index_contig_names(std::ifstream& fasta_index)
{
    std::vector<std::string> region_names;
    std::string line;
    while (std::getline(fasta_index, line, '\n')) {
        region_names.emplace_back(line.substr(0, line.find('\t')));
    }
    return region_names;
}

inline
std::vector<std::string>
get_fasta_index_contig_names(const std::string& fasta_index_path)
{
    std::ifstream fasta_index(fasta_index_path, std::ios::binary);
    return get_fasta_index_contig_names(fasta_index);
}

inline
std::unordered_map<std::string, FastaIndex>
read_fasta_index(const std::string& fasta_index_path)
{
    std::ifstream fasta_index(fasta_index_path, std::ios::binary);
    std::unordered_map<std::string, FastaIndex> index;
    std::string line;
    while (std::getline(fasta_index, line, '\n')) {
        FastaIndex fai(line);
        index[fai.contig_name] = std::move(fai);
    }
    return index;
}

inline
size_t get_contig_size(std::ifstream& fasta_index, const std::string& contig_name)
{
    fasta_index.seekg(0, std::ios::beg);
    std::string line;
    while (std::getline(fasta_index, line, '\n')) {
        FastaIndex contig_index(line);
        if (contig_index.contig_name == contig_name) return contig_index.length;
    }
    return 0;
}

inline
size_t get_contig_size(const std::string& fasta_index_path, const std::string& contig_name)
{
    std::ifstream fasta_index(fasta_index_path, std::ios::binary);
    return get_contig_size(fasta_index, contig_name);
}

/*=======================================================================================
    FASTA: For reading FASTAs with index files
 =======================================================================================*/

inline
size_t get_line_offset(const FastaIndex& index, size_t begin)
{
    return begin % index.line_byte_length;
}

inline
size_t get_contig_offset(const FastaIndex& index, size_t begin)
{
    return index.offset + index.line_length * (begin / index.line_byte_length) +
                get_line_offset(index, begin);
}

inline
size_t get_remaining_line_length(const FastaIndex& index, size_t begin)
{
    return index.line_length - get_line_offset(index, begin);
}

template <typename T=std::string>
T read_fasta_contig(std::ifstream& fasta, const FastaIndex& index, size_t begin, size_t length)
{
    fasta.seekg(get_contig_offset(index, begin), std::ios::beg);
    T seq {};
    seq.resize(length);
    if (index.line_length == index.line_byte_length) {
        fasta.read(&seq[0], length);
    } else {
        size_t line_length_to_read {std::min(length, get_remaining_line_length(index, begin))};
        fasta.read(&seq[0], line_length_to_read);
        size_t num_line_end_bytes {index.line_byte_length - index.line_length};
        fasta.ignore(num_line_end_bytes);
        size_t num_chars_read {line_length_to_read};
        while (num_chars_read < length) {
            fasta.read(&seq[num_chars_read], index.line_length);
            fasta.ignore(num_line_end_bytes);
            num_chars_read += index.line_length;
        }
    }
    return seq;
}

template <typename T=std::string>
T read_fasta_contig(std::ifstream& fasta, const FastaIndex& index)
{
    return read_fasta_contig(fasta, index, 0, index.length);
}

template <typename T=std::string>
T read_fasta_contig(const std::string& fasta_path, const FastaIndex& index)
{
    std::ifstream fasta(fasta_path, std::ios::binary | std::ios::beg);
    return read_fasta_contig(fasta, index);
}

template <typename T1=std::string, typename T2=std::string>
void write_fasta(const std::string& fasta_path, FastaRecord<T1, T2> data)
{
    std::ofstream fasta(fasta_path, std::ios::out | std::ios::binary);
    fasta << ">" << std::move(data.first) << "\n" << std::move(data.second);
}

/*=======================================================================================
 FASTA: For reading multiple record Fasta files.
 =======================================================================================*/

template<typename T1=std::string, typename T2=std::string>
FastaRecord<T1, T2>
read_fasta_record(std::ifstream& fasta)
{
    T1 contig_name;
    std::getline(fasta, contig_name); // contig_name always a single line
    
    T2 line;
    std::getline(fasta, line);
    
    // The FASTA format is not as simple as FASTQ - the sequence
    // may be broken into multiple lines.
    if (fasta.peek() == '>' || !fasta.good()) {
        return FastaRecord<T1, T2>(contig_name, line);
    } else {
        line.shrink_to_fit();
        size_t line_size = line.size();
        
        T2 seq;
        auto it = std::back_inserter(seq);
        std::copy(line.begin(), line.end(), it);
        
        while (fasta.peek() != '>' && fasta.good()) {
            fasta.getline(&line[0], line_size + 1);
            std::copy(line.begin(), line.begin() + fasta.gcount(), it);
        }
        
        return FastaRecord<T1, T2>(contig_name, seq);
    }
}

template<typename T=std::string>
std::vector<T>
read_fasta_seqs(std::string path)
{
    std::ifstream fasta(path, std::ios::binary);
    std::vector<T> seqs;
    size_t num_records = get_num_records(fasta, '>');
    seqs.reserve(num_records);
    while (num_records) {
        seqs.push_back(read_fasta_record<std::string, T>(fasta).seq);
        --num_records;
    }
    return seqs;
}

template<typename T1=std::string, typename T2=std::string>
std::vector<FastaRecord<T1, T2>>
read_fasta(std::string path)
{
    std::ifstream fasta(path, std::ios::binary);
    std::vector<FastaRecord<T1, T2>> seqs;
    size_t num_records = get_num_records(fasta, '>');
    seqs.reserve(num_records);
    while (num_records) {
        seqs.push_back(read_fasta_record<T1, T2>(fasta));
        --num_records;
    }
    return seqs;
}

template <typename T1=std::string, typename T2=std::string, typename F>
FastaReads<T1, T2>
read_fasta_map(std::string path, F f)
{
   std::ifstream fasta(path, std::ios::binary);
   FastaMap<T1, T2> records;
   ReadIds<T1> contig_names;
   size_t num_records = get_num_records(fasta, '>');
   records.reserve(num_records);
   contig_names.reserve(num_records);
   while (num_records) {
       auto record = read_fasta_record<T1, T2>(fasta);
       auto f_contig_name = f(record.contig_name);
       records.emplace(f_contig_name, std::move(record.seq));
       contig_names.insert(std::move(f_contig_name));
       --num_records;
   }
   return std::make_pair(contig_names, records);
}

template<typename T1=std::string, typename T2=std::string>
FastaReads<T1, T2>
read_fasta_map(std::string path)
{
    return read_fasta_map(path, [] (T1 contig_name) { return contig_name; });
}

template <typename T1=std::string, typename T2=std::string, typename F>
FastaReads<T1, T2>
read_fasta_map(std::string path, const ReadIds<T1>& contig_names, F f)
{
   std::ifstream fasta(path, std::ios::binary);
   FastaMap<T1, T2> records;
   ReadIds<T1> f_contig_names;
   size_t num_records = contig_names.size();
   records.reserve(num_records);
   f_contig_names.reserve(num_records);
   while (num_records) {
       auto record = read_fasta_record<T1, T2>(fasta);
       auto f_contig_name = f(record.contig_name);
       if (contig_names.find(f_contig_name) != contig_names.end()) {
           records.emplace(f_contig_name, std::move(record.seq));
           f_contig_names.insert(std::move(f_contig_name));
           --num_records;
       }
   }
   return {f_contig_names, records};
}

template<typename T1=std::string, typename T2=std::string>
FastaReads<T1, T2>
read_fasta_map(std::string path, const ReadIds<T1>& contig_names)
{
    return read_fasta_map(path, contig_names, [] (T1 contig_name) { return contig_name; });
}

template<typename T1=std::string, typename T2=std::string>
int write_fasta(std::string path, const FastaReads<T1, T2>& data)
{
    std::ofstream fasta(path, std::ios::out | std::ios::binary);
    if (fasta) {
        for (auto contig_name : data.first) {
            fasta << ">" << contig_name << "\n";
            fasta << data.second.at(contig_name);
        }
        return 0;
    }
    return 1;
}

/*=======================================================================================
 FASTQ: For reading multiple line Fastq files.
 =======================================================================================*/

template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
FastqRecord<T1, T2, T3>
read_fastq_record(std::ifstream& fastq)
{
    T1 contig_name;
    T2 seq;
    T3 quals;
    // Unlike FASTA, FASTQ always use one line per field.
    std::getline(fastq, contig_name);
    std::getline(fastq, seq);
    fastq.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::getline(fastq, quals);
    return FastqRecord<T1, T2, T3>(std::move(contig_name), std::move(seq), std::move(quals));
}

template<typename T=std::string>
std::vector<T>
read_fastq_seqs(std::string path)
{
    std::ifstream fastq(path, std::ios::binary);
    std::vector<T> seqs;
    size_t num_records = get_num_records(fastq, '@');
    seqs.reserve(num_records);
    while (num_records) {
        seqs.push_back(read_fastq_record<std::string, T, std::string>(fastq).seq);
        --num_records;
    }
    return seqs;
}

template<typename T1=std::string, typename T2=std::string, typename T3=std::string, typename F>
std::vector<FastqRecord<T1, T2, T3>>
read_fastq(std::string path, F f)
{
    std::ifstream fastq(path, std::ios::binary);
    std::vector<FastqRecord<T1, T2, T3>> data;
    size_t num_records = get_num_records(fastq, '@');
    data.reserve(num_records);
    while (num_records) {
        auto record = read_fastq_record<T1, T2, T3>(fastq);
        data.emplace_back(f(record.contig_name), std::move(record.seq), std::move(record.qual));
        --num_records;
    }
    return data;
}

template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
std::vector<FastqRecord<T1, T2, T3>>
read_fastq(std::string path)
{
    return read_fastq(path, [] (T1 contig_name) { return contig_name; });
}

template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
FastqReads<T1, T2, T3>
read_fastq_map(std::string path)
{
   std::ifstream fastq(path, std::ios::binary);
   ReadIds<T1> contig_names;
   FastqMap<T1, T2, T3> data;
   size_t num_records = get_num_records(fastq, '@');
   data.reserve(num_records);
   while (num_records) {
       auto record = read_fastq_record<T1, T2, T3>(fastq);
       data.emplace(record.contig_name, {std::move(record.seq), std::move(record.qual)});
       contig_names.insert(std::move(record.contig_name));
       --num_records;
   }
   return {contig_names, data};
}

}

#endif /* defined(__bioio__bioio__) */
