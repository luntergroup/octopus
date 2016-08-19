// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caller_factory.hpp"

#include <utility>

#include <io/reference/reference_genome.hpp>
#include <readpipe/read_pipe.hpp>

namespace octopus {

CallerFactory::CallerFactory(CallerBuilder template_builder, const unsigned default_ploidy)
: template_builder_ {std::move(template_builder)}
, contig_ploidies_ {}
, default_ploidy_ {default_ploidy}
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

CallerFactory& CallerFactory::set_contig_ploidy(const ContigName& contig, const unsigned ploidy)
{
    contig_ploidies_[contig] = ploidy;
    return *this;
}

std::unique_ptr<Caller> CallerFactory::make(const ContigName& contig) const
{
    if (contig_ploidies_.count(contig) == 1) {
        template_builder_.set_ploidy(contig_ploidies_.at(contig));
    } else {
        template_builder_.set_ploidy(default_ploidy_);
    }
    return template_builder_.build();
}

} // namespace octopus
