// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef conditional_random_forest_filter_hpp
#define conditional_random_forest_filter_hpp

#include <vector>
#include <cstddef>
#include <memory>
#include <fstream>
#include <functional>
#include <cstddef>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "ranger/Forest.h"

#include "double_pass_variant_call_filter.hpp"

namespace octopus { namespace csr {

class ConditionalRandomForestFilter : public DoublePassVariantCallFilter
{
public:
    using Path = boost::filesystem::path;
    
    ConditionalRandomForestFilter() = delete;
    
    ConditionalRandomForestFilter(FacetFactory facet_factory,
                                  std::vector<MeasureWrapper> measures,
                                  std::vector<MeasureWrapper> chooser_measures,
                                  std::function<std::int8_t(std::vector<Measure::ResultType>)> chooser,
                                  std::vector<Path> ranger_forests,
                                  OutputOptions output_config,
                                  ConcurrencyPolicy threading,
                                  Path temp_directory = "/tmp",
                                  boost::optional<ProgressMeter&> progress = boost::none);
    
    ConditionalRandomForestFilter(const ConditionalRandomForestFilter&)            = delete;
    ConditionalRandomForestFilter& operator=(const ConditionalRandomForestFilter&) = delete;
    ConditionalRandomForestFilter(ConditionalRandomForestFilter&&)                 = default;
    ConditionalRandomForestFilter& operator=(ConditionalRandomForestFilter&&)      = default;
    
    virtual ~ConditionalRandomForestFilter() override = default;

private:
    struct File
    {
        std::ofstream handle;
        Path path;
        template <typename F, typename P>
        File(F&& handle, P&& path) : handle {std::forward<F>(handle)}, path {std::forward<P>(path)} {};
    };
    
    std::vector<Path> forest_paths_;
    Path temp_dir_;
    std::vector<std::unique_ptr<ranger::Forest>> forests_;
    std::function<std::int8_t(std::vector<Measure::ResultType>)> chooser_;
    std::size_t num_chooser_measures_;
    
    mutable std::vector<std::vector<File>> data_;
    mutable std::size_t num_records_;
    mutable std::vector<std::vector<std::vector<double>>> data_buffer_;
    mutable std::vector<std::deque<std::int8_t>> choices_;
    mutable std::deque<std::size_t> hard_filtered_record_indices_;
    mutable std::vector<bool> hard_filtered_;
    
    const static std::string call_qual_name_;
    
    boost::optional<std::string> call_quality_name() const override;
    void annotate(VcfHeader::Builder& header) const override;
    std::int8_t choose_forest(const MeasureVector& measures) const;
    void prepare_for_registration(const SampleList& samples) const override;
    void record(std::size_t call_idx, std::size_t sample_idx, MeasureVector measures) const override;
    void close_data_files() const;
    void prepare_for_classification(boost::optional<Log>& log) const override;
    std::size_t get_forest_choice(std::size_t call_idx, std::size_t sample_idx) const;
    Classification classify(std::size_t call_idx, std::size_t sample_idx) const override;
};

} // namespace csr
} // namespace octopus

#endif
