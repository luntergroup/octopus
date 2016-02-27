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

#include <iostream>

#include <boost/timer/timer.hpp> // BENCHMARK

// variant caller timers
extern boost::timer::cpu_timer init_timer;
extern boost::timer::cpu_timer haplotype_generation_timer;
extern boost::timer::cpu_timer likelihood_timer;
extern boost::timer::cpu_timer haplotype_fitler_timer;
extern boost::timer::cpu_timer prior_model_timer;
extern boost::timer::cpu_timer latent_timer;
extern boost::timer::cpu_timer phasing_timer;
extern boost::timer::cpu_timer allele_generator_timer;
extern boost::timer::cpu_timer calling_timer;

// population model timers
extern boost::timer::cpu_timer genotype_generation_timer;
extern boost::timer::cpu_timer genotype_likelihood_timer;
extern boost::timer::cpu_timer prior_count_timer;
extern boost::timer::cpu_timer frequency_init_timer;
extern boost::timer::cpu_timer marginal_init_timer;
extern boost::timer::cpu_timer posterior_init_timer;
extern boost::timer::cpu_timer frequency_update_timer;
extern boost::timer::cpu_timer marginal_update_timer;
extern boost::timer::cpu_timer posterior_update_timer;
extern boost::timer::cpu_timer em_timer;

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

inline void init_timers()
{
    init_timer.stop();
    haplotype_generation_timer.stop();
    likelihood_timer.stop();
    haplotype_fitler_timer.stop();
    prior_model_timer.stop();
    latent_timer.stop();
    phasing_timer.stop();
    allele_generator_timer.stop();
    calling_timer.stop();
    
    genotype_generation_timer.stop();
    genotype_likelihood_timer.stop();
    prior_count_timer.stop();
    frequency_init_timer.stop();
    marginal_init_timer.stop();
    posterior_init_timer.stop();
    frequency_update_timer.stop();
    marginal_update_timer.stop();
    posterior_update_timer.stop();
    em_timer.stop();
}

inline void print_all_timers()
{
    std::cout << "init timer" << '\n';
    std::cout << init_timer.format() << std::endl;
    
    std::cout << "haplotype generation timer" << '\n';
    std::cout << haplotype_generation_timer.format() << std::endl;
    
    std::cout << "likelihood timer" << '\n';
    std::cout << likelihood_timer.format() << std::endl;
    
    std::cout << "haplotype fitler timer" << '\n';
    std::cout << haplotype_fitler_timer.format() << std::endl;
    
    std::cout << "prior model timer timer" << '\n';
    std::cout << prior_model_timer.format() << std::endl;
    
    std::cout << "latent timer" << '\n';
    std::cout << latent_timer.format() << std::endl;
    
    std::cout << "phasing timer" << '\n';
    std::cout << phasing_timer.format() << std::endl;
    
    std::cout << "allele generator timer" << '\n';
    std::cout << allele_generator_timer.format() << std::endl;
    
    std::cout << "calling timer" << '\n';
    std::cout << calling_timer.format() << std::endl;
    
    std::cout << "genotype generation timer" << '\n';
    std::cout << calling_timer.format() << std::endl;
    
    std::cout << "genotype likelihood timer" << '\n';
    std::cout << genotype_likelihood_timer.format() << std::endl;
    
    std::cout << "prior count timer" << '\n';
    std::cout << prior_count_timer.format() << std::endl;
    
    std::cout << "freqency init timer" << '\n';
    std::cout << frequency_init_timer.format() << std::endl;
    
    std::cout << "marginal init timer" << '\n';
    std::cout << marginal_init_timer.format() << std::endl;
    
    std::cout << "posterior init timer" << '\n';
    std::cout << posterior_init_timer.format() << std::endl;
    
    std::cout << "frequency update timer" << '\n';
    std::cout << frequency_update_timer.format() << std::endl;
    
    std::cout << "marginal update timer" << '\n';
    std::cout << marginal_update_timer.format() << std::endl;
    
    std::cout << "posterior update timer" << '\n';
    std::cout << posterior_update_timer.format() << std::endl;
    
    std::cout << "em timer" << '\n';
    std::cout << em_timer.format() << std::endl;
}

inline void print_caller_timers()
{
    std::cout << "init timer" << '\n';
    std::cout << init_timer.format() << std::endl;
    
    std::cout << "haplotype generation timer" << '\n';
    std::cout << haplotype_generation_timer.format() << std::endl;
    
    std::cout << "likelihood timer" << '\n';
    std::cout << likelihood_timer.format() << std::endl;
    
    std::cout << "haplotype fitler timer" << '\n';
    std::cout << haplotype_fitler_timer.format() << std::endl;
    
    std::cout << "prior model timer timer" << '\n';
    std::cout << prior_model_timer.format() << std::endl;
    
    std::cout << "latent timer" << '\n';
    std::cout << latent_timer.format() << std::endl;
    
    std::cout << "phasing timer" << '\n';
    std::cout << phasing_timer.format() << std::endl;
    
    std::cout << "allele generator timer" << '\n';
    std::cout << allele_generator_timer.format() << std::endl;
    
    std::cout << "calling timer" << '\n';
    std::cout << calling_timer.format() << std::endl;
}

inline void print_model_timers()
{
    std::cout << "genotype generation timer" << '\n';
    std::cout << calling_timer.format() << std::endl;
    
    std::cout << "genotype likelihood timer" << '\n';
    std::cout << genotype_likelihood_timer.format() << std::endl;
    
    std::cout << "prior count timer" << '\n';
    std::cout << prior_count_timer.format() << std::endl;
    
    std::cout << "freqency init timer" << '\n';
    std::cout << frequency_init_timer.format() << std::endl;
    
    std::cout << "marginal init timer" << '\n';
    std::cout << marginal_init_timer.format() << std::endl;
    
    std::cout << "posterior init timer" << '\n';
    std::cout << posterior_init_timer.format() << std::endl;
    
    std::cout << "frequency update timer" << '\n';
    std::cout << frequency_update_timer.format() << std::endl;
    
    std::cout << "marginal update timer" << '\n';
    std::cout << marginal_update_timer.format() << std::endl;
    
    std::cout << "posterior update timer" << '\n';
    std::cout << posterior_update_timer.format() << std::endl;
    
    std::cout << "em timer" << '\n';
    std::cout << em_timer.format() << std::endl;
}

#endif /* timers_h */
