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

boost::timer::cpu_timer misc_timer1 {};
boost::timer::cpu_timer misc_timer2 {};

void init_timers()
{
    init_timer.start(); init_timer.stop();
    haplotype_generation_timer.start(); haplotype_generation_timer.stop();
    haplotype_likelihood_timer.start(); haplotype_likelihood_timer.stop();
    haplotype_fitler_timer.start(); haplotype_fitler_timer.stop();
    latent_timer.start(); latent_timer.stop();
    phasing_timer.start(); phasing_timer.stop();
    calling_timer.start(); calling_timer.stop();
    
    misc_timer1.start(); misc_timer1.stop();
    misc_timer2.start(); misc_timer2.stop();
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
    
    std::cout << "misc timer 1" << '\n';
    std::cout << misc_timer1.format() << std::endl;
    
    std::cout << "misc timer 2" << '\n';
    std::cout << misc_timer2.format() << std::endl;
}
