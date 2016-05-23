//
//  timers.cpp
//  Octopus
//
//  Created by Daniel Cooke on 24/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "timers.hpp"

#include <iostream>

boost::timer::cpu_timer init_timer {};
boost::timer::cpu_timer haplotype_generation_timer {};
boost::timer::cpu_timer haplotype_likelihood_timer {};
boost::timer::cpu_timer haplotype_fitler_timer {};
boost::timer::cpu_timer latent_timer {};
boost::timer::cpu_timer calling_timer {};
boost::timer::cpu_timer phasing_timer {};
boost::timer::cpu_timer output_timer {};

TimerArray misc_timer {};

void init_timers()
{
    init_timer.start(); init_timer.stop();
    haplotype_generation_timer.start(); haplotype_generation_timer.stop();
    haplotype_likelihood_timer.start(); haplotype_likelihood_timer.stop();
    haplotype_fitler_timer.start(); haplotype_fitler_timer.stop();
    latent_timer.start(); latent_timer.stop();
    phasing_timer.start(); phasing_timer.stop();
    calling_timer.start(); calling_timer.stop();
    
    for (auto& timer : misc_timer) {
        timer.stop();
    }
}

void print_all_timers()
{
    std::cout << "init timer" << '\n';
    std::cout << init_timer.format() << std::endl;
    
    std::cout << "haplotype generation timer" << '\n';
    std::cout << haplotype_generation_timer.format() << std::endl;
    
    std::cout << "likelihood timer" << '\n';
    std::cout << haplotype_likelihood_timer.format() << std::endl;
    
    std::cout << "haplotype fitler timer" << '\n';
    std::cout << haplotype_fitler_timer.format() << std::endl;
    
    std::cout << "latent timer" << '\n';
    std::cout << latent_timer.format() << std::endl;
    
    std::cout << "calling timer" << '\n';
    std::cout << calling_timer.format() << std::endl;
    
    std::cout << "phasing timer" << '\n';
    std::cout << phasing_timer.format() << std::endl;
    
    for (std::size_t i {0}; i < misc_timer.size(); ++i) {
        std::cout << "misc timer " << i << '\n';
        std::cout << misc_timer[i].format() << '\n';
    }
}
