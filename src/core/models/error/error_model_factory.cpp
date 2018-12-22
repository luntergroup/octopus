// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_model_factory.hpp"

#include <fstream>
#include <sstream>

#include "utils/string_utils.hpp"
#include "hiseq_snv_error_model.hpp"
#include "x10_snv_error_model.hpp"
#include "hiseq_indel_error_model.hpp"
#include "x10_indel_error_model.hpp"
#include "novaseq_indel_error_model.hpp"
#include "bgiseq_indel_error_model.hpp"
#include "custom_repeat_based_indel_error_model.hpp"

namespace octopus {

std::unique_ptr<SnvErrorModel> make_snv_error_model()
{
    return std::make_unique<HiSeqSnvErrorModel>();
}

std::unique_ptr<IndelErrorModel> make_indel_error_model()
{
    return std::make_unique<HiSeqIndelErrorModel>();
}

bool is_xten(const std::string& model_name) noexcept
{
    return model_name == "XTEN" || model_name == "X10";
}

bool is_novaseq(const std::string& model_name) noexcept
{
    return model_name == "NOVASEQ";
}

bool is_bgiseq(const std::string& model_name) noexcept
{
    return model_name == "BGISEQ" || model_name == "BGISEQ-500" || model_name == "BGISEQ500";
}

bool is_hiseq(const std::string& model_name) noexcept
{
    return model_name == "HISEQ";
}

std::unique_ptr<SnvErrorModel> make_snv_error_model(const std::string& model_name)
{
    auto model_code_name = utils::capitalise(model_name);
    if (is_xten(model_code_name)) {
        return std::make_unique<X10SnvErrorModel>();
    } else {
        return std::make_unique<HiSeqSnvErrorModel>();
    }
}

std::unique_ptr<IndelErrorModel> make_indel_error_model(const std::string& model_name)
{
    auto model_code_name = utils::capitalise(model_name);
    if (is_xten(model_code_name)) {
        return std::make_unique<X10IndelErrorModel>();
    } else if (is_novaseq(model_code_name)) {
        return std::make_unique<NovaSeqIndelErrorModel>();
    } else if (is_bgiseq(model_code_name)) {
        return std::make_unique<BGISeqIndelErrorModel>();
    } else if (is_hiseq(model_code_name)) {
        return std::make_unique<HiSeqIndelErrorModel>();
    } else {
        std::ifstream model_file {model_name};
        std::string model_str {static_cast<std::stringstream const&>(std::stringstream() << model_file.rdbuf()).str()};
        auto open_model = make_penalty_map(std::move(model_str));
        return std::make_unique<CustomRepeatBasedIndelErrorModel>(std::move(open_model), 3);
    }
}

} // namespace octopus
