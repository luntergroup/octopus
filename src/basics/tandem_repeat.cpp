// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "tandem_repeat.hpp"

namespace octopus {

unsigned count_periods(const TandemRepeat& repeat) noexcept
{
    return region_size(repeat) / repeat.period();
}

bool operator==(const TandemRepeat& lhs, const TandemRepeat& rhs) noexcept
{
    return is_same_region(lhs, rhs) && lhs.period() == rhs.period();
}

bool operator<(const TandemRepeat& lhs, const TandemRepeat& rhs) noexcept
{
    return is_same_region(lhs, rhs) ? lhs.period() < rhs.period() : mapped_region(lhs) < mapped_region(rhs);
}

} // namespace octopus
