// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "call_wrapper.hpp"

#include <algorithm>

#include "basics/genomic_region.hpp"
#include "call.hpp"

namespace octopus {

const GenomicRegion& CallWrapper::mapped_region() const noexcept
{
    return call->mapped_region();
}

std::vector<std::unique_ptr<Call>> unwrap(std::deque<CallWrapper>&& calls)
{
    std::vector<std::unique_ptr<Call>> result {};
    result.reserve(calls.size());
    std::transform(std::make_move_iterator(std::begin(calls)),
                   std::make_move_iterator(std::end(calls)),
                   std::back_inserter(result),
                   [] (CallWrapper&& wrapped_call) -> std::unique_ptr<Call> {
                       return std::move(wrapped_call.call);
                   });
    calls.clear();
    calls.shrink_to_fit();
    return result;
}

bool operator==(const CallWrapper& lhs, const CallWrapper& rhs)
{
    return lhs.call == rhs.call;
}

CallWrapper clone(const CallWrapper& call)
{
    return CallWrapper {call->clone()};
}

} // namespace octopus
