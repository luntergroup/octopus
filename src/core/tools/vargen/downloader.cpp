// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "downloader.hpp"

//#include <boost/network/protocol/http.hpp>

#include "basics/aligned_read.hpp"

namespace octopus { namespace coretools {

Downloader::Downloader(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
{}

std::unique_ptr<VariantGenerator> Downloader::do_clone() const
{
    return std::make_unique<Downloader>(*this);
}

std::vector<Variant> Downloader::do_generate(const RegionSet& regions, OptionalThreadPool workers) const
{
    //namespace http = boost::network::http;
    
//    <?xml version="1.0" encoding="UTF-8"?>
//    <!DOCTYPE Query>
//    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
//    
//    <Dataset name = "hsapiens_snp" interface = "default" >
//    <Filter name = "chr_name" value = "X"/>
//    <Filter name = "end" value = "10500"/>
//    <Filter name = "start" value = "10000"/>
//    <Attribute name = "refsnp_id" />
//    <Attribute name = "refsnp_source" />
//    <Attribute name = "chr_name" />
//    <Attribute name = "chrom_start" />
//    <Attribute name = "chrom_end" />
//    <Attribute name = "allele" />
//    </Dataset>
//    </Query>
    
    try {
//        http::client::request request("http://www.boost.org/");
//        http::client client;
//        http::client::response response = client.get(request);
//        std::cout << body(response) << std::endl;
        
        std::vector<Variant> result {};
        return result;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw;
    }
}
    
std::string Downloader::name() const
{
    return "Download";
}

} // namespace coretools
} // namespace octopus
