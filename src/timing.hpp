//
//  timing.hpp
//  Octopus
//
//  Created by Daniel Cooke on 20/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef timing_hpp
#define timing_hpp

#include <chrono>
#include <ctime>
#include <iomanip>

inline std::ostream& operator<<(std::ostream& os, const std::chrono::system_clock::time_point& t)
{
    const std::time_t t_c {std::chrono::system_clock::to_time_t(t)};
    os << std::put_time(std::localtime(&t_c), "%F %T");
    return os;
}

struct TimeInterval
{
    using TimePoint = std::chrono::system_clock::time_point;
    
    TimeInterval() = default;
    TimeInterval(TimePoint start, TimePoint end) : start {start}, end {end} {}
    
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
        os << duration_s.count() << "s";
        return os;
    }
    
    const auto duration_m = duration<std::chrono::minutes>(interval);
    
    if (duration_m.count() < 60) {
        os << duration_m.count() << "." << ((100 * (duration_s.count() % 60)) / 60) << "m";
    } else {
        const auto duration_h = duration<std::chrono::hours>(interval);
        os << duration_h.count() << "." << ((100 * (duration_m.count() % 60)) / 60) << "h";
    }
    
    return os;
}

#endif /* timing_hpp */
