//
//  htslib_bcf_facade.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__htslib_bcf_facade__
#define __Octopus__htslib_bcf_facade__

#include <string>
#include <set>
#include <memory> // std::unique_ptr
#include <boost/filesystem/path.hpp>

#include "variant_file_impl.h"
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

class GenomicRegion;
class Variant;

namespace fs = boost::filesystem;

auto htslib_file_deleter       = [] (bcf_srs_t* the_file) { bcf_sr_destroy(the_file); };
auto htslib_bcf_header_deleter = [] (bcf_hdr_t* the_header) { bcf_hdr_destroy(the_header); };
auto htslib_bcf1_deleter       = [] (bcf1_t* the_bcf1) { bcf_destroy(the_bcf1); };

class HtslibBcfFacade : public IVariantFileImpl
{
public:
    HtslibBcfFacade() = delete;
    HtslibBcfFacade(const fs::path& the_variant_file_path);
    ~HtslibBcfFacade() = default;
    
    HtslibBcfFacade(const HtslibBcfFacade&)            = default;
    HtslibBcfFacade& operator=(const HtslibBcfFacade&) = default;
    HtslibBcfFacade(HtslibBcfFacade&&)                 = default;
    HtslibBcfFacade& operator=(HtslibBcfFacade&&)      = default;
    
    std::vector<Variant> fetch_variants(const GenomicRegion& a_region) override;
    void write_variants(const std::vector<Variant>& some_variants) override;
    
private:
    fs::path the_file_path_;
    std::unique_ptr<bcf_srs_t, decltype(htslib_file_deleter)> the_file_;
    std::unique_ptr<bcf_hdr_t, decltype(htslib_bcf_header_deleter)> the_header_;
    
    void set_region(const GenomicRegion& a_region);
};

#endif /* defined(__Octopus__htslib_bcf_facade__) */
