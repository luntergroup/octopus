// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef timing_hpp
#define timing_hpp

#include <chrono>
#include <ctime>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>

namespace octopus { namespace utils
{
inline std::ostream& operator<<(std::ostream& os, const std::chrono::system_clock::time_point& t)
{
    const std::time_t t_c {std::chrono::system_clock::to_time_t(t)};
    os << std::put_time(std::localtime(&t_c), "%F %T");
    return os;
}

struct TimeInterval
{
    using TimePoint = std::chrono::system_clock::time_point;
    TimePoint start, end;
};

template <typename T>
auto duration(const TimeInterval& interval)
{
    return std::chrono::duration_cast<T>(interval.end - interval.start);
}

inline std::ostream& operator<<(std::ostream& os, const TimeInterval& interval)
{
    const auto duration_ms = duration<std::chrono::milliseconds>(interval);
    if (duration_ms.count() < 1000) {
        os << duration_ms.count() << "ms";
        return os;
    }
    const auto duration_s = duration<std::chrono::seconds>(interval);
    const auto secs = duration_s.count();
    if (secs < 60) {
        os << secs << 's';
        return os;
    }
    const auto duration_m = duration<std::chrono::minutes>(interval);
    const auto mins = duration_m.count();
    if (mins < 60) {
        os << mins << 'm';
        const auto remainder_secs = secs % 60;
        if (remainder_secs > 0) os << ' ' << remainder_secs << 's';
    } else {
        const auto duration_h = duration<std::chrono::hours>(interval);
        if (duration_h.count() <= 24) {
            const auto hours = duration_h.count();
            os << hours << 'h';
            const auto remainder_mins = mins % 60;
            if (remainder_mins > 0) os << ' ' << remainder_mins << 'm';
        } else {
            using H = std::chrono::hours::rep;
            constexpr H num_hours_in_day {24};
            const auto days = std::div(duration_h.count(), num_hours_in_day);
            if (days.quot < 7) {
                os << days.quot << 'd';
                if (days.rem > 0) os << ' ' << days.rem << 'h';
            } else {
                constexpr H num_hours_in_week {7 * num_hours_in_day};
                const auto weeks = std::div(duration_h.count(), num_hours_in_week);
                os << weeks.quot << 'w';
                const auto rem_days = weeks.rem / num_hours_in_day;
                if (rem_days > 0) os << ' ' << rem_days << 'd';
            }
        }
    }
    return os;
}
} // namespace utils
} // namespace octopus

#endif
