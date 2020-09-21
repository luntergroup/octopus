// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef vcf_writer_hpp
#define vcf_writer_hpp

#include <memory>
#include <mutex>
#include <type_traits>
#include <functional>
#include <iterator>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "htslib_bcf_facade.hpp"

namespace octopus {

class VcfHeader;
class VcfRecord;
class GenomicRegion;

class VcfWriter
{
public:
    using Path = boost::filesystem::path;
    
    VcfWriter();
    VcfWriter(Path file_path);
    VcfWriter(const VcfHeader& header);
    VcfWriter(Path file_path, const VcfHeader& header);
    
    VcfWriter(const VcfWriter&)            = delete;
    VcfWriter& operator=(const VcfWriter&) = delete;
    VcfWriter(VcfWriter&&);
    VcfWriter& operator=(VcfWriter&&);
    
    ~VcfWriter();
    
    friend void swap(VcfWriter& lhs, VcfWriter& rhs) noexcept;
    
    bool is_open() const noexcept;
    void open();
    void open(Path file_path);
    void close() noexcept;
    
    bool is_header_written() const noexcept;
    
    boost::optional<Path> path() const;
    
    void write(const VcfHeader& header);
    void write(const VcfRecord& record);
    
private:
    boost::optional<Path> file_path_;
    std::unique_ptr<HtslibBcfFacade> writer_;
    bool is_header_written_;
    mutable std::mutex mutex_;
    
    bool can_write_index() const noexcept;
};

VcfWriter& operator<<(VcfWriter& dst, const VcfHeader& header);
VcfWriter& operator<<(VcfWriter& dst, const VcfRecord& record);

template <typename Container>
void write(const Container& records, VcfWriter& dst)
{
    static_assert(std::is_same<typename Container::value_type, VcfRecord>::value, "");
    for (const auto& record : records) {
        dst << record;
    }
}

template <typename Container>
VcfWriter& operator<<(VcfWriter& dst, const Container& records)
{
    write(records, dst);
    return dst;
}

class VcfWriterIterator
{
public:
    using iterator_category = std::output_iterator_tag;
    using value_type        = void;
    using difference_type   = void;
    using pointer           = void;
    using reference         = void;
    
    VcfWriterIterator() = delete;
    
    VcfWriterIterator(VcfWriter& writer) : writer_ {writer} {}
    
    VcfWriterIterator& operator=(const VcfRecord& record)
    {
        writer_ << record;
        return *this;
    }
    
    VcfWriterIterator& operator*() { return *this; }
    VcfWriterIterator& operator++() { return *this; }
    VcfWriterIterator& operator++(int) { return *this; }
    
private:
    std::reference_wrapper<VcfWriter> writer_;
};

bool operator==(const VcfWriter& lhs, const VcfWriter& rhs);

} // namespace octopus

namespace std {
    template <> struct hash<octopus::VcfWriter>
    {
        size_t operator()(const octopus::VcfWriter& writer) const
        {
            const auto output_path = writer.path();
            return output_path ? hash<string>()(output_path->string()) : 0;
        }
    };
} // namespace std

#endif
