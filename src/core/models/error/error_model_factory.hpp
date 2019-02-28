// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef error_model_factory_hpp
#define error_model_factory_hpp

#include <memory>
#include <string>

#include <boost/filesystem/path.hpp>

#include "snv_error_model.hpp"
#include "indel_error_model.hpp"

namespace octopus {

enum class LibraryPreparation { pcr, pcr_free, tenx };
enum class Sequencer { hiseq_2000, hiseq_2500, hiseq_4000, xten, novaseq, bgiseq_500 };

struct ModelConfig
{
    LibraryPreparation library;
    Sequencer sequencer;
};

static constexpr LibraryPreparation default_library_preparation {LibraryPreparation::pcr_free};
static constexpr Sequencer default_sequencer {Sequencer::hiseq_2500};
static constexpr ModelConfig default_model_config {default_library_preparation, default_sequencer};

struct ErrorModel
{
    std::unique_ptr<IndelErrorModel> indel;
    std::unique_ptr<SnvErrorModel> snv;
};

std::unique_ptr<SnvErrorModel> make_snv_error_model(ModelConfig config = default_model_config);
std::unique_ptr<IndelErrorModel> make_indel_error_model(ModelConfig config = default_model_config);

ErrorModel make_error_model(const std::string& label);
ErrorModel make_error_model(const boost::filesystem::path& model_file_name);

} // namespace octopus

#endif
