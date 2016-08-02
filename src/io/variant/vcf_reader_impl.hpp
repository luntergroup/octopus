//
//  vcf_reader_impl.hpp
//  Octopus
//
//  Created by Daniel Cooke on 14/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_vcf_reader_impl_hpp
#define Octopus_vcf_reader_impl_hpp

#include <vector>
#include <string>
#include <cstddef>
#include <memory>
#include <utility>

namespace octopus {

class GenomicRegion;
class VcfHeader;
class VcfRecord;

class IVcfReaderImpl
{
public:
    enum class UnpackPolicy { All, Sites };
    
    using RecordContainer = std::vector<VcfRecord>;
    
    class RecordIterator
    {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = VcfRecord;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const VcfRecord*;
        using reference         = const VcfRecord&;
        
        virtual reference operator*() const = 0;
        virtual pointer operator->() const = 0;
        
        virtual void next() = 0;
        
        virtual ~RecordIterator() = default;
    };
    
    using RecordIteratorPtr     = std::unique_ptr<RecordIterator>;
    using RecordIteratorPtrPair = std::pair<RecordIteratorPtr, RecordIteratorPtr>;
    
    virtual bool is_header_written() const noexcept = 0;
    
    virtual VcfHeader fetch_header() const = 0;
    
    virtual std::size_t count_records() const = 0;
    virtual std::size_t count_records(const std::string& contig) const = 0;
    virtual std::size_t count_records(const GenomicRegion& region) const = 0;
    
    virtual RecordContainer fetch_records(UnpackPolicy level = UnpackPolicy::All) const = 0; // fetches all records
    virtual RecordContainer fetch_records(const std::string& contig, UnpackPolicy level = UnpackPolicy::All) const = 0;
    virtual RecordContainer fetch_records(const GenomicRegion& region, UnpackPolicy level = UnpackPolicy::All) const = 0;
    
    virtual RecordIteratorPtrPair iterate(UnpackPolicy level = UnpackPolicy::All) const = 0;
    
    virtual ~IVcfReaderImpl() noexcept = default;
};

} // namespace octopus    

#endif
