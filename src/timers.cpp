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
boost::timer::cpu_timer em_timer {};
