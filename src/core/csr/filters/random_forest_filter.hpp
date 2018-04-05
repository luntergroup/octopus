// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef random_forest_filter_hpp
#define random_forest_filter_hpp

#include <vector>
#include <cstddef>
#include <memory>
#include <fstream>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "ranger/Forest.h"

#include "double_pass_variant_call_filter.hpp"

namespace octopus { namespace csr {

class RandomForestFilter : public DoublePassVariantCallFilter
{
public:
    using Path = boost::filesystem::path;
    
    RandomForestFilter() = delete;
    
    RandomForestFilter(FacetFactory facet_factory,
                       std::vector<MeasureWrapper> measures,
                       OutputOptions output_config,
                       ConcurrencyPolicy threading,
                       Path ranger_forest,
                       Path temp_directory = "/tmp",
                       boost::optional<ProgressMeter&> progress = boost::none);
    
    RandomForestFilter(const RandomForestFilter&)            = delete;
    RandomForestFilter& operator=(const RandomForestFilter&) = delete;
    RandomForestFilter(RandomForestFilter&&)                 = default;
    RandomForestFilter& operator=(RandomForestFilter&&)      = default;
    
    virtual ~RandomForestFilter() override = default;

private:
    struct File
    {
        std::ofstream handle;
        Path path;
        template <typename F, typename P>
        File(F&& handle, P&& path) : handle {std::forward<F>(handle)}, path {std::forward<P>(path)} {};
    };
    
    std::unique_ptr<Forest> forest_;
    Path ranger_forest_, temp_dir_;
    
    mutable std::vector<File> data_;
    mutable std::size_t num_records_;
    mutable std::vector<std::vector<double>> data_buffer_;
    
    void annotate(VcfHeader::Builder& header) const override;
    void prepare_for_registration(const SampleList& samples) const override;
    void record(std::size_t call_idx, std::size_t sample_idx, MeasureVector measures) const override;
    void prepare_for_classification(boost::optional<Log>& log) const override;
    Classification classify(std::size_t call_idx, std::size_t sample_idx) const override;
};

} // namespace csr
} // namespace octopus

#endif
