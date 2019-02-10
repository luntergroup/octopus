// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef measure_factory_hpp
#define measure_factory_hpp

#include <string>

#include "measure.hpp"

namespace octopus { namespace csr {

MeasureWrapper make_measure(std::string name);

std::vector<std::string> get_all_measure_names();

} // namespace csr
} // namespace octopus

#endif
