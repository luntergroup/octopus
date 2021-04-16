// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef call_wrapper_hpp
#define call_wrapper_hpp

#include <memory>
#include <vector>
#include <deque>
#include <iterator>

#include "concepts/equitable.hpp"
#include "concepts/mappable.hpp"

namespace octopus {

class GenomicRegion;
class Call;

// Wrap the pointer so can use mappable algorithms
struct CallWrapper : public Equitable<CallWrapper>, public Mappable<CallWrapper>
{
    CallWrapper(std::unique_ptr<Call> call) : call {std::move(call) } {}
    operator const std::unique_ptr<Call>&()  { return call; }
    operator std::unique_ptr<Call>&() noexcept  { return call; }
    std::unique_ptr<Call>::pointer operator->() const noexcept { return call.get(); }
    const GenomicRegion& mapped_region() const noexcept;
    std::unique_ptr<Call> call;
};

template <typename T>
auto wrap(std::vector<std::unique_ptr<T>>&& calls)
{
    return std::vector<CallWrapper> {std::make_move_iterator(std::begin(calls)), std::make_move_iterator(std::end(calls))};
}

std::vector<std::unique_ptr<Call>> unwrap(std::deque<CallWrapper>&& calls);

bool operator==(const CallWrapper& lhs, const CallWrapper& rhs);

CallWrapper clone(const CallWrapper& call);

} // namespace octopus

#endif
