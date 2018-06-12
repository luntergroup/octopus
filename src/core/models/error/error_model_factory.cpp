// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_model_factory.hpp"

#include "hiseq_snv_error_model.hpp"
#include "x10_snv_error_model.hpp"
#include "hiseq_indel_error_model.hpp"
#include "x10_indel_error_model.hpp"

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
    return sequencer == "xTen" || sequencer == "x10";
}

std::unique_ptr<SnvErrorModel> make_snv_error_model(const std::string& sequencer)
{
    if (is_xten(sequencer)) {
        return std::make_unique<X10SnvErrorModel>();
    }
    return std::make_unique<HiSeqSnvErrorModel>();
}

std::unique_ptr<IndelErrorModel> make_indel_error_model(const std::string& sequencer)
{
    if (is_xten(sequencer)) {
        return std::make_unique<X10IndelErrorModel>();
    }
    return std::make_unique<HiSeqIndelErrorModel>();
}
    
} // namespace octopus
