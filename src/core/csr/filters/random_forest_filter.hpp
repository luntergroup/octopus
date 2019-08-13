// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef random_forest_filter_hpp
#define random_forest_filter_hpp

#include <vector>
#include <cstddef>
#include <memory>
#include <fstream>
#include <functional>
#include <cstddef>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "ranger/Forest.h"

#include "basics/phred.hpp"
#include "double_pass_variant_call_filter.hpp"

namespace octopus { namespace csr {

class RandomForestFilter : public DoublePassVariantCallFilter
{
public:
    using Path = boost::filesystem::path;
    
    struct Options
    {
        Phred<double> min_forest_quality = probability_false_to_phred(0.5);
    };
    
    RandomForestFilter() = delete;
    
    RandomForestFilter(FacetFactory facet_factory,
                       std::vector<MeasureWrapper> measures,
                       Path ranger_forest,
                       OutputOptions output_config,
                       ConcurrencyPolicy threading,
                       Path temp_directory,
                       Options options,
                       boost::optional<ProgressMeter&> progress = boost::none);
    
    RandomForestFilter(FacetFactory facet_factory,
                       std::vector<MeasureWrapper> measures,
                       std::vector<MeasureWrapper> chooser_measures,
                       std::function<std::int8_t(std::vector<Measure::ResultType>)> chooser,
                       std::vector<Path> ranger_forests,
                       OutputOptions output_config,
                       ConcurrencyPolicy threading,
                       Path temp_directory,
                       Options options,
                       boost::optional<ProgressMeter&> progress = boost::none);
    
    RandomForestFilter(const RandomForestFilter&)            = delete;
    RandomForestFilter& operator=(const RandomForestFilter&) = delete;
    RandomForestFilter(RandomForestFilter&&)                 = delete;
    RandomForestFilter& operator=(RandomForestFilter&&)      = delete;
    
    virtual ~RandomForestFilter() override = default;

protected:
    virtual void annotate(VcfHeader::Builder& header) const override;
    
    Phred<double> min_soft_genotype_quality() const noexcept;
    Phred<double> min_soft_call_quality() const noexcept;

private:
    struct File
    {
        std::ofstream handle;
        Path path;
        template <typename F, typename P>
        File(F&& handle, P&& path) : handle {std::forward<F>(handle)}, path {std::forward<P>(path)} {};
    };
    
    std::vector<Path> forest_paths_;
    std::vector<std::unique_ptr<ranger::Forest>> forests_;
    std::function<std::int8_t(std::vector<Measure::ResultType>)> chooser_;
    std::size_t num_chooser_measures_;
    Options options_;
    
    mutable std::vector<std::vector<File>> data_;
    mutable std::size_t num_records_;
    mutable std::vector<std::vector<std::vector<double>>> data_buffer_;
    mutable std::vector<std::deque<std::int8_t>> choices_;
    mutable std::deque<std::size_t> hard_filtered_record_indices_;
    mutable std::vector<bool> hard_filtered_;
    
    const static std::string genotype_quality_name_;
    const static std::string call_quality_name_;
    
    std::string do_name() const override;
    virtual boost::optional<std::string> call_quality_name() const override;
    virtual bool is_soft_filtered(const ClassificationList& sample_classifications, boost::optional<Phred<double>> joint_quality,
                                  const MeasureVector& measures, std::vector<std::string>& reasons) const override;
    
    std::unique_ptr<ranger::Forest> make_forest() const;
    boost::optional<std::string> genotype_quality_name() const override;
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
