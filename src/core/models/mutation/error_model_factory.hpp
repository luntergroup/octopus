// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef error_model_factory_hpp
#define error_model_factory_hpp

#include <memory>
#include <string>

#include "snv_error_model.hpp"
#include "indel_error_model.hpp"

namespace octopus {

std::unique_ptr<SnvErrorModel> make_snv_error_model(const std::string& sequencer);
std::unique_ptr<IndelErrorModel> make_indel_error_model(const std::string& sequencer);

} // namespace octopus

#endif
