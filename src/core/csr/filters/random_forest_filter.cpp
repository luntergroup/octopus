// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "random_forest_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cassert>

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>

#include "ranger/ForestProbability.h"

#include "basics/phred.hpp"

namespace octopus { namespace csr {

RandomForestFilter::RandomForestFilter(FacetFactory facet_factory,
                                       std::vector<MeasureWrapper> measures,
                                       OutputOptions output_config,
                                       ConcurrencyPolicy threading,
                                       Path ranger_forest, Path temp_directory,
                                       boost::optional<ProgressMeter&> progress)
: DoublePassVariantCallFilter {std::move(facet_factory), std::move(measures), std::move(output_config), threading, progress}
, forest_ {std::make_unique<ForestProbability>()}
, ranger_forest_ {std::move(ranger_forest)}
, temp_dir_ {std::move(temp_directory)}
, num_records_ {0}
, data_buffer_ {}
{}

void RandomForestFilter::annotate(VcfHeader::Builder& header) const
{
    header.add_filter("RF", "Random Forest filtered");
}

void RandomForestFilter::prepare_for_registration(const SampleList& samples) const
{
    std::vector<std::string> data_header {};
    data_header.reserve(measures_.size());
    for (const auto& measure : measures_) {
        data_header.push_back(measure.name());
    }
    data_header.push_back("TP");
    data_.reserve(samples.size());
    for (const auto& sample : samples) {
        auto data_path = temp_dir_;
        Path fname {"octopus_ranger_temp_forest_data_" + sample + ".dat"};
        data_path /= fname;
        data_.emplace_back(data_path.string(), data_path);
        std::copy(std::cbegin(data_header), std::prev(std::cend(data_header)),
                  std::ostream_iterator<std::string> {data_.back().handle, " "});
        data_.back().handle << data_header.back() << '\n';
    }
    data_buffer_.resize(samples.size());
}

namespace {

struct MeasureDoubleVisitor : boost::static_visitor<>
{
    double result;
    template <typename T> void operator()(const T& value)
    {
        result = boost::lexical_cast<double>(value);
    }
    template <typename T> void operator()(const boost::optional<T>& value)
    {
        if (value) {
            (*this)(*value);
        } else {
            result = 0;
        }
    }
    template <typename T> void operator()(const std::vector<T>& values)
    {
        throw std::runtime_error {"Vector cast not supported"};
    }
    void operator()(boost::any value)
    {
        throw std::runtime_error {"Any cast not supported"};
    }
};

auto cast_to_double(const Measure::ResultType& value)
{
    MeasureDoubleVisitor vis {};
    boost::apply_visitor(vis, value);
    return vis.result;
}

void skip_lines(std::istream& in, int n = 1)
{
    for (; n > 0; --n) in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

} // namespace

void RandomForestFilter::record(const std::size_t call_idx, std::size_t sample_idx, MeasureVector measures) const
{
    assert(!measures.empty());
    std::transform(std::cbegin(measures), std::cend(measures), std::back_inserter(data_buffer_[sample_idx]), cast_to_double);
    std::copy(std::cbegin(data_buffer_[sample_idx]), std::cend(data_buffer_[sample_idx]),
              std::ostream_iterator<double> {data_[sample_idx].handle, " "});
    data_[sample_idx].handle << 0 << '\n'; // dummy TP value
    data_buffer_[sample_idx].clear();
    if (call_idx >= num_records_) ++num_records_;
}

void RandomForestFilter::prepare_for_classification(boost::optional<Log>& log) const
{
    const Path ranger_prefix {temp_dir_ / "octopus_ranger_temp"};
    const Path ranger_prediction_fname {ranger_prefix.string() + ".prediction"};
    data_buffer_.resize(num_records_);
    for (auto& file : data_) {
        file.handle.close();
        std::vector<std::string> tmp {}, cat_vars {};
        forest_->initCpp("TP", MemoryMode::MEM_DOUBLE, file.path.string(), 0, ranger_prefix.string(),
                         1000, nullptr, 12, 1, ranger_forest_.string(), ImportanceMode::IMP_GINI, 1, "", tmp, "", true,
                         cat_vars, false, SplitRule::LOGRANK, "", false, 1.0, DEFAULT_ALPHA, DEFAULT_MINPROP, false,
                         PredictionType::RESPONSE, DEFAULT_NUM_RANDOM_SPLITS);
        forest_->run(false);
        forest_->writePredictionFile();
        std::ifstream prediction_file {ranger_prediction_fname.string()};
        skip_lines(prediction_file, 3); // header
        std::string line;
        std::size_t i {0};
        while (std::getline(prediction_file, line)) {
            line.erase(std::find(std::cbegin(line), std::cend(line), ' '), std::cend(line));
            data_buffer_[i++].push_back(boost::lexical_cast<double>(line));
        }
        boost::filesystem::remove(file.path);
    }
    boost::filesystem::remove(ranger_prediction_fname);
    data_.clear();
    data_.shrink_to_fit();
}

VariantCallFilter::Classification RandomForestFilter::classify(const std::size_t call_idx, std::size_t sample_idx) const
{
    assert(call_idx < data_buffer_.size() && sample_idx < data_buffer_[call_idx].size());
    const auto prob_true = data_buffer_[call_idx][sample_idx];
    Classification result {};
    if (prob_true >= 0.5) {
        result.category = Classification::Category::unfiltered;
    } else {
        result.category = Classification::Category::soft_filtered;
        result.reasons.assign({"RF"});
    }
    result.quality = probability_to_phred(1 - prob_true);
    return result;
}

} // namespace csr
} // namespace octopus
