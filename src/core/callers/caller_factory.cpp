// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caller_factory.hpp"

#include <utility>

#include "io/reference/reference_genome.hpp"
#include "readpipe/read_pipe.hpp"

namespace octopus {

CallerFactory::CallerFactory(CallerBuilder template_builder)
: template_builder_ {std::move(template_builder)}
{}

CallerFactory& CallerFactory::set_reference(const ReferenceGenome& reference) noexcept
{
    template_builder_.set_reference(reference);
    return *this;
}

CallerFactory& CallerFactory::set_read_pipe(ReadPipe& read_pipe) noexcept
{
    template_builder_.set_read_pipe(read_pipe);
    return *this;
}

std::unique_ptr<Caller> CallerFactory::make(const ContigName& contig) const
{
    return template_builder_.build(contig);
}

} // namespace octopus
