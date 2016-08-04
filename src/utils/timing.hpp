// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef timing_hpp
#define timing_hpp

#include <chrono>
#include <ctime>
#include <iomanip>

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
    
    if (duration_s.count() < 60) {
        os << duration_s.count() << 's';
        return os;
    }
    
    const auto duration_m = duration<std::chrono::minutes>(interval);
    
    if (duration_m.count() < 60) {
        os << duration_m.count() << '.' << std::setw(2) << std::setfill('0') << ((100 * (duration_s.count() % 60)) / 60) << 'm';
    } else {
        const auto duration_h = duration<std::chrono::hours>(interval);
        os << duration_h.count() << '.' << std::setw(2) << std::setfill('0') << ((100 * (duration_m.count() % 60)) / 60) << 'h';
    }
    
    return os;
}
} // namespace utils
} // namespace octopus

#endif /* timing_hpp */
