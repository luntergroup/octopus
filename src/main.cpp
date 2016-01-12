//
//  main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MODULE Main
//#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <utility>

#include "program_options.hpp"
#include "octopus.hpp"

#include "mock_options.hpp"

std::ostream& operator<<(std::ostream& os, const std::chrono::system_clock::time_point& t)
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

std::ostream& operator<<(std::ostream& os, TimeInterval interval)
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
        os << duration_m.count() << "." << (duration_s.count() % 60) << "m";
    } else {
        const auto duration_h = duration<std::chrono::hours>(interval);
        os << duration_h.count() << "." << (duration_m.count() % 60) << "h";
    }
    
    return os;
}

int main(int argc, const char **argv)
{
    using std::cout; using std::cerr; using std::endl;
    
    try {
        //const auto options = Octopus::parse_options(argc, argv);
        const auto options = get_basic_mock_options();
        
        if (options) {
            auto start = std::chrono::system_clock::now();
            
            cout << "Started running Octopus at " << start << endl;
            
            Octopus::run_octopus(options.get());
            
            auto end = std::chrono::system_clock::now();
            
            cout << "Finished running Octopus at " << end << endl;
            
            cout << "Elapsed time: " << TimeInterval {start, end} << endl;
        } else {
            cout << "Could not parse input options. Did not run Octopus." << endl;
        }
    } catch (std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
//    catch (...) {
//        cerr << "Error: encountered unknown error. Quiting now" << endl;
//        return EXIT_FAILURE;
//    }
    
    return EXIT_SUCCESS;
}
