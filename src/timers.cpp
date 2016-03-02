//
//  timers.cpp
//  Octopus
//
//  Created by Daniel Cooke on 24/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "timers.hpp"

boost::timer::cpu_timer init_timer {};
boost::timer::cpu_timer haplotype_generation_timer {};
boost::timer::cpu_timer likelihood_timer {};
boost::timer::cpu_timer haplotype_fitler_timer {};
boost::timer::cpu_timer prior_model_timer {};
boost::timer::cpu_timer latent_timer {};
boost::timer::cpu_timer phasing_timer {};
boost::timer::cpu_timer allele_generator_timer {};
boost::timer::cpu_timer calling_timer {};

// population model timers
boost::timer::cpu_timer genotype_generation_timer {};
boost::timer::cpu_timer genotype_likelihood_timer {};
boost::timer::cpu_timer prior_count_timer {};
boost::timer::cpu_timer frequency_init_timer {};
boost::timer::cpu_timer marginal_init_timer {};
boost::timer::cpu_timer posterior_init_timer {};
boost::timer::cpu_timer frequency_update_timer {};
boost::timer::cpu_timer marginal_update_timer {};
boost::timer::cpu_timer posterior_update_timer {};
boost::timer::cpu_timer em_timer {};

boost::timer::cpu_timer kmer_mapping_timer {};
boost::timer::cpu_timer likelihood_cache_timer {};

boost::timer::cpu_timer misc_timer1 {};
boost::timer::cpu_timer misc_timer2 {};
