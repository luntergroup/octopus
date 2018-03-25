// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef timers_hpp
#define timers_hpp

//#define BENCHMARK

#include <array>

#include <boost/timer/timer.hpp>

// variant caller timers
extern boost::timer::cpu_timer init_timer;
extern boost::timer::cpu_timer haplotype_likelihood_timer;
extern boost::timer::cpu_timer latent_timer;
extern boost::timer::cpu_timer calling_timer;
extern boost::timer::cpu_timer phasing_timer;
extern boost::timer::cpu_timer output_timer;

using TimerArray = std::array<boost::timer::cpu_timer, 12>;

extern TimerArray misc_timer;

inline void resume(boost::timer::cpu_timer& timer)
{
    #ifdef BENCHMARK
    timer.resume();
    #endif
}

inline void pause(boost::timer::cpu_timer& timer)
{
    #ifdef BENCHMARK
    timer.stop();
    #endif
}

void init_timers();

void print_all_timers();

#endif
