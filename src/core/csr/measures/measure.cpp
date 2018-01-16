// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "measure.hpp"

#include <boost/lexical_cast.hpp>

namespace octopus { namespace csr {

struct MeasureSerialiseVisitor : boost::static_visitor<>
{
    std::string str;
    
    template <typename T>
    void operator()(const T& value)
    {
        str = boost::lexical_cast<std::string>(value);
    }
    template <typename T>
    void operator()(const boost::optional<T>& value)
    {
        if (value) {
            str = boost::lexical_cast<std::string>(*value);
        } else {
            str = ".";
        }
    }
};

std::string Measure::do_serialise(const ResultType& value) const
{
    MeasureSerialiseVisitor vis {};
    boost::apply_visitor(vis, value);
    return vis.str;
}

} // namespace csr
} // namespace octopus
