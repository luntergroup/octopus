// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_model_factory.hpp"

#include "utils/string_utils.hpp"
#include "hiseq_snv_error_model.hpp"
#include "x10_snv_error_model.hpp"
#include "hiseq_indel_error_model.hpp"
#include "x10_indel_error_model.hpp"
#include "novaseq_indel_error_model.hpp"
#include "bgiseq_indel_error_model.hpp"

namespace octopus {

std::unique_ptr<SnvErrorModel> make_snv_error_model()
{
    return std::make_unique<HiSeqSnvErrorModel>();
}

std::unique_ptr<IndelErrorModel> make_indel_error_model()
{
    return std::make_unique<HiSeqIndelErrorModel>();
}

bool is_xten(const std::string& sequencer) noexcept
{
    return sequencer == "XTEN" || sequencer == "X10";
}

bool is_novaseq(const std::string& sequencer) noexcept
{
    return sequencer == "NOVASEQ";
}

bool is_bgiseq(const std::string& sequencer) noexcept
{
    return sequencer == "BGISEQ" || sequencer == "BGISEQ-500" || sequencer == "BGISEQ500";
}

std::unique_ptr<SnvErrorModel> make_snv_error_model(std::string sequencer)
{
    sequencer = utils::capitalise(sequencer);
    if (is_xten(sequencer)) {
        return std::make_unique<X10SnvErrorModel>();
    } else {
        return std::make_unique<HiSeqSnvErrorModel>();
    }
}

std::unique_ptr<IndelErrorModel> make_indel_error_model(std::string sequencer)
{
    sequencer = utils::capitalise(sequencer);
    if (is_xten(sequencer)) {
        return std::make_unique<X10IndelErrorModel>();
    } else if (is_novaseq(sequencer)) {
        return std::make_unique<NovaSeqIndelErrorModel>();
    } else if (is_bgiseq(sequencer)) {
        return std::make_unique<BGISeqIndelErrorModel>();
    } else {
        return std::make_unique<HiSeqIndelErrorModel>();
    }
}

} // namespace octopus
