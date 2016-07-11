//
//  timers.hpp
//  Octopus
//
//  Created by Daniel Cooke on 24/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef timers_h
#define timers_h

#define BENCHMARK

#include <array>

#include <boost/timer/timer.hpp>

// variant caller timers
extern boost::timer::cpu_timer init_timer;
extern boost::timer::cpu_timer haplotype_generation_timer;
extern boost::timer::cpu_timer haplotype_likelihood_timer;
extern boost::timer::cpu_timer haplotype_fitler_timer;
extern boost::timer::cpu_timer latent_timer;
extern boost::timer::cpu_timer calling_timer;
extern boost::timer::cpu_timer phasing_timer;
extern boost::timer::cpu_timer output_timer;

using TimerArray = std::array<boost::timer::cpu_timer, 12>;

extern TimerArray misc_timer;

inline void resume_timer(boost::timer::cpu_timer& timer)
{
    #ifdef BENCHMARK
    timer.resume();
    #endif
}

inline void pause_timer(boost::timer::cpu_timer& timer)
{
    #ifdef BENCHMARK
    timer.stop();
    #endif
}

void init_timers();

void print_all_timers();

#endif /* timers_h */
