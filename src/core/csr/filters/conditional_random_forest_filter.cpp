// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "conditional_random_forest_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>

#include "ranger/ForestProbability.h"

#include "utils/concat.hpp"
#include "utils/maths.hpp"
#include "exceptions/missing_file_error.hpp"
#include "exceptions/program_error.hpp"
#include "exceptions/malformed_file_error.hpp"

namespace octopus { namespace csr {

namespace {

class MissingForestFile : public MissingFileError
{
    std::string do_where() const override { return "ConditionalRandomForestFilter"; }
public:
    MissingForestFile(boost::filesystem::path p) : MissingFileError {std::move(p), ".forest"} {};
};

void check_all_exists(const std::vector<ConditionalRandomForestFilter::Path>& forests)
{
    for (const auto& forest : forests) {
        if (!boost::filesystem::exists(forest)) {
            throw MissingForestFile {forest};
        }
    }
}

} // namespace

ConditionalRandomForestFilter::ConditionalRandomForestFilter(FacetFactory facet_factory,
                                                             std::vector<MeasureWrapper> measures,
                                                             std::vector<MeasureWrapper> chooser_measures,
                                                             std::function<std::int8_t(std::vector<Measure::ResultType>)> chooser,
                                                             std::vector<Path> ranger_forests,
                                                             OutputOptions output_config,
                                                             ConcurrencyPolicy threading,
                                                             Path temp_directory,
                                                             boost::optional<ProgressMeter&> progress)
: ConditionalRandomForestFilter {std::move(facet_factory), std::move(measures), std::move(chooser_measures),
                                 std::move(chooser), std::move(ranger_forests), probability_to_phred(0.5),
                                 std::move(output_config), threading, std::move(temp_directory), progress} {}

ConditionalRandomForestFilter::ConditionalRandomForestFilter(FacetFactory facet_factory,
                                                             std::vector<MeasureWrapper> measures,
                                                             std::vector<MeasureWrapper> chooser_measures,
                                                             std::function<std::int8_t(std::vector<Measure::ResultType>)> chooser,
                                                             std::vector<Path> ranger_forests,
                                                             Phred<double> min_forest_quality,
                                                             OutputOptions output_config,
                                                             ConcurrencyPolicy threading,
                                                             Path temp_directory,
                                                             boost::optional<ProgressMeter&> progress)
: DoublePassVariantCallFilter {std::move(facet_factory), concat(std::move(measures), chooser_measures),
                               std::move(output_config), threading, std::move(temp_directory), progress}
, forest_paths_ {std::move(ranger_forests)}
, chooser_ {std::move(chooser)}
, num_chooser_measures_ {chooser_measures.size()}
, min_forest_quality_ {min_forest_quality}
, num_records_ {0}
, data_buffer_ {}
{
    check_all_exists(forest_paths_);
    forests_.reserve(forest_paths_.size());
    std::generate_n(std::back_inserter(forests_), forest_paths_.size(),
                    [] () { return std::make_unique<ranger::ForestProbability>(); });
}

const std::string ConditionalRandomForestFilter::genotype_quality_name_ = "RFQUAL";

boost::optional<std::string> ConditionalRandomForestFilter::genotype_quality_name() const
{
    return genotype_quality_name_;
}

void ConditionalRandomForestFilter::annotate(VcfHeader::Builder& header) const
{
    header.add_format(genotype_quality_name_, "1", "Float", "Empirical quality score from random forest classifier");
    header.add_filter("RF", "Random Forest filtered");
}

std::int8_t ConditionalRandomForestFilter::choose_forest(const MeasureVector& measures) const
{
    const MeasureVector chooser_measures(std::prev(std::cend(measures), num_chooser_measures_), std::cend(measures));
    return chooser_(chooser_measures);
}

template <typename T>
static void write_line(const std::vector<T>& data, std::ostream& out)
{
    std::copy(std::cbegin(data), std::prev(std::cend(data)), std::ostream_iterator<T> {out, " "});
    out << data.back() << '\n';
}

void ConditionalRandomForestFilter::prepare_for_registration(const SampleList& samples) const
{
    std::vector<std::string> data_header {};
    data_header.reserve(measures_.size() - num_chooser_measures_);
    std::transform(std::cbegin(measures_), std::prev(std::cend(measures_), num_chooser_measures_), std::back_inserter(data_header),
                   [] (const auto& measure) { return measure.name(); });
    data_header.push_back("TP");
    const auto num_forests = forest_paths_.size();
    data_.resize(num_forests);
    for (std::size_t forest_idx {0}; forest_idx < num_forests; ++forest_idx) {
        data_[forest_idx].reserve(samples.size());
        for (const auto& sample : samples) {
            auto data_path = temp_directory();
            Path fname {"octopus_ranger_temp_forest_data_" + std::to_string(forest_idx) + "_" + sample + ".dat"};
            data_path /= fname;
            data_[forest_idx].emplace_back(data_path.string(), data_path);
            write_line(data_header, data_[forest_idx].back().handle);
        }
    }
    data_buffer_.resize(num_forests);
    for (auto& buffer : data_buffer_) buffer.resize(samples.size());
    choices_.resize(samples.size());
}

namespace {

template <typename T>
bool is_subnormal(const T x) noexcept
{
    return std::fpclassify(x) == FP_SUBNORMAL;
}

template <typename T>
double lexical_cast_to_double(const T& value)
{
    auto result = boost::lexical_cast<double>(value);
    if (is_subnormal(result)) {
        result = 0;
    }
    return result;
}

struct MeasureDoubleVisitor : boost::static_visitor<>
{
    double result;
    template <typename T> void operator()(const T& value)
    {
        result = lexical_cast_to_double(value);
    }
    template <typename T> void operator()(const boost::optional<T>& value)
    {
        if (value) {
            (*this)(*value);
        } else {
            result = -1;
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

class NanMeasure : public ProgramError
{
    std::string do_where() const override { return "ConditionalRandomForestFilter"; }
    std::string do_why() const override { return "detected a nan measure"; }
    std::string do_help() const override { return "submit an error report"; }
};

void check_nan(const std::vector<double>& values)
{
    if (std::any_of(std::cbegin(values), std::cend(values), [] (auto v) { return std::isnan(v); })) {
        throw NanMeasure {};
    }
}

void skip_lines(std::istream& in, int n = 1)
{
    for (; n > 0; --n) in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

} // namespace

void ConditionalRandomForestFilter::record(const std::size_t call_idx, std::size_t sample_idx, MeasureVector measures) const
{
    assert(!measures.empty());
    const auto forest_idx = choose_forest(measures);
    const auto num_forests = static_cast<std::remove_const_t<decltype(forest_idx)>>(data_buffer_.size());
    if (forest_idx >= 0 && forest_idx < num_forests) {
        auto& buffer = data_buffer_[forest_idx][sample_idx];
        std::transform(std::cbegin(measures), std::prev(std::cend(measures), num_chooser_measures_),
                       std::back_inserter(buffer), cast_to_double);
        buffer.push_back(0); // dummy TP value
        check_nan(buffer);
        write_line(buffer, data_[forest_idx][sample_idx].handle);
        buffer.clear();
    } else {
        hard_filtered_record_indices_.push_back(call_idx);
    }
    if (call_idx >= num_records_) ++num_records_;
    choices_[sample_idx].push_back(forest_idx);
}

void ConditionalRandomForestFilter::close_data_files() const
{
    for (auto& forest : data_) {
        for (auto& sample : forest) {
            sample.handle.close();
        }
    }
}

namespace {

bool read_header(std::ifstream& prediction_file)
{
    skip_lines(prediction_file);
    std::string order;
    std::getline(prediction_file, order);
    skip_lines(prediction_file);
    return order.front() == '1';
}

static double get_prob_false(std::string& prediction_line, const bool tp_first)
{
    using std::cbegin; using std::cend;
    if (tp_first) {
        prediction_line.erase(cbegin(prediction_line), std::next(std::find(cbegin(prediction_line), cend(prediction_line), ' ')));
        prediction_line.erase(std::find(cbegin(prediction_line), cend(prediction_line), ' '), cend(prediction_line));
    } else {
        prediction_line.erase(std::find(cbegin(prediction_line), cend(prediction_line), ' '), cend(prediction_line));
    }
    return boost::lexical_cast<double>(prediction_line);
}

} // namespace

class MalformedForestFile : public MalformedFileError
{
    std::string do_where() const override { return "ConditionalRandomForestFilter"; }
    std::string do_help() const override
    {
        return "make sure the forest was trained with the same measures and in the same order as the prediction measures";
    }
public:
    MalformedForestFile(boost::filesystem::path file) : MalformedFileError {std::move(file)} {}
};

void ConditionalRandomForestFilter::prepare_for_classification(boost::optional<Log>& log) const
{
    close_data_files();
    if (num_records_ == 0) return;
    const Path ranger_prefix {temp_directory() / "octopus_ranger_temp"};
    const Path ranger_prediction_fname {ranger_prefix.string() + ".prediction"};
    data_buffer_.resize(1);
    auto& predictions = data_buffer_[0];
    predictions.resize(num_records_);
    const auto num_samples = choices_.size();
    for (std::size_t forest_idx {0}; forest_idx < forest_paths_.size(); ++forest_idx) {
        for (std::size_t sample_idx {0}; sample_idx < num_samples; ++sample_idx) {
            auto forest_choice_itr = std::find(std::cbegin(choices_[sample_idx]), std::cend(choices_[sample_idx]), forest_idx);
            if (forest_choice_itr != std::cend(choices_[sample_idx])) {
                const auto& file = data_[forest_idx][sample_idx];
                std::vector<std::string> tmp {}, cat_vars {};
                auto& forest = forests_[forest_idx];
                try {
                    forest->initCpp("TP", ranger::MemoryMode::MEM_DOUBLE, file.path.string(), 0, ranger_prefix.string(),
                                    1000, nullptr, 12, 1, forest_paths_[forest_idx].string(), ranger::ImportanceMode::IMP_GINI, 1, "",
                                    tmp, "", true, cat_vars, false, ranger::SplitRule::LOGRANK, "", false, 1.0,
                                    ranger::DEFAULT_ALPHA, ranger::DEFAULT_MINPROP, false,
                                    ranger::PredictionType::RESPONSE, ranger::DEFAULT_NUM_RANDOM_SPLITS);
                } catch (const std::runtime_error& e) {
                    throw MalformedForestFile {forest_paths_[forest_idx]};
                }
                forest->run(false);
                forest->writePredictionFile();
                std::ifstream prediction_file {ranger_prediction_fname.string()};
                const auto tp_first = read_header(prediction_file);
                std::string line;
                while (std::getline(prediction_file, line)) {
                    if (!line.empty()) {
                        const auto record_idx = std::distance(std::cbegin(choices_[sample_idx]), forest_choice_itr);
                        predictions[record_idx].push_back(get_prob_false(line, tp_first));
                        assert(forest_choice_itr != std::cend(choices_[sample_idx]));
                        forest_choice_itr = std::find(std::next(forest_choice_itr), std::cend(choices_[sample_idx]), forest_idx);
                    }
                }
                boost::filesystem::remove(file.path);
            }
        }
    }
    boost::filesystem::remove(ranger_prediction_fname);
    data_.clear();
    data_.shrink_to_fit();
    choices_.clear();
    choices_.shrink_to_fit();
    if (!hard_filtered_record_indices_.empty()) {
        hard_filtered_.resize(num_records_, false);
        for (auto idx : hard_filtered_record_indices_) {
            hard_filtered_[idx] = true;
        }
        hard_filtered_record_indices_.clear();
        hard_filtered_record_indices_.shrink_to_fit();
    }
}

std::size_t ConditionalRandomForestFilter::get_forest_choice(std::size_t call_idx, std::size_t sample_idx) const
{
    return choices_.empty() ? 0 : choices_[sample_idx][call_idx];
}

VariantCallFilter::Classification ConditionalRandomForestFilter::classify(const std::size_t call_idx, std::size_t sample_idx) const
{
    Classification result {};
    if (hard_filtered_.empty() || !hard_filtered_[call_idx]) {
        const auto& predictions = data_buffer_[0];
        assert(call_idx < predictions.size() && sample_idx < predictions[call_idx].size());
        const auto prob_false = predictions[call_idx][sample_idx];
        result.quality = probability_to_phred(std::max(prob_false, 1e-10));
        if (result.quality >= min_forest_quality_) {
            result.category = Classification::Category::unfiltered;
        } else {
            result.category = Classification::Category::soft_filtered;
            result.reasons.assign({"RF"});
        }
    } else {
        result.category = Classification::Category::hard_filtered;
    }
    return result;
}

} // namespace csr
} // namespace octopus
