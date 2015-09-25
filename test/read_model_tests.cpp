//
//  read_model_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <string>
//#include <cstddef>
//#include <unordered_map>
//
//#include "test_common.hpp"
//#include "reference_genome.hpp"
//#include "read_manager.hpp"
//#include "allele.hpp"
//#include "variant.hpp"
//#include "variant_utils.hpp"
//#include "candidate_variant_generator.hpp"
//#include "alignment_candidate_variant_generator.hpp"
//#include "haplotype.hpp"
//#include "genotype.hpp"
//#include "read_model.hpp"
//
//using std::cout;
//using std::endl;
//
//using Octopus::ReadModel;

//BOOST_AUTO_TEST_CASE(partially_overlapped_reads_evaluate_correctly)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    auto a_region = parse_region("11:27282186-27282290", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    auto reads = a_read_manager.fetch_reads(samples[0], a_region);
//    
//    //cout << reads.size() << endl;
//    
////    std::vector<AlignedRead> read_3_copies(25, reads[3]);
////    reads.insert(reads.end(), read_3_copies.cbegin(), read_3_copies.cend());
//    
//    Haplotype reference_haplotype {human, a_region};
//    
//    Haplotype true_haplotype {human, a_region};
//    true_haplotype.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    Haplotype false_haplotype1 {human, a_region};
//    false_haplotype1.push_back(Allele {parse_region("11:27282258-27282259", human), "C"});
//    false_haplotype1.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    Haplotype false_haplotype2 {human, a_region};
//    false_haplotype2.push_back(Allele {parse_region("11:27282199-27282200", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282207-27282208", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282218-27282219", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282239-27282240", human), "T"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    unsigned ploidy {2};
//    
//    ReadModel a_read_model {ploidy};
//    
//    //cout << reads[4] << endl;
//    
////    cout << a_read_model.log_probability(reads[3], reference_haplotype, 0) << endl;
////    cout << a_read_model.log_probability(reads[3], true_haplotype, 0) << endl;
////    cout << a_read_model.log_probability(reads[3], false_haplotype1, 0) << endl;
////    cout << a_read_model.log_probability(reads[3], false_haplotype2, 0) << endl;
////    cout << endl;
//    
////    cout << a_read_model.log_probability(reads[7], reference_haplotype, 0) << endl;
////    cout << a_read_model.log_probability(reads[7], true_haplotype, 0) << endl;
////    cout << a_read_model.log_probability(reads[7], false_haplotype1, 0) << endl;
////    cout << a_read_model.log_probability(reads[7], false_haplotype2, 0) << endl;
////    cout << endl;
//    
////    for (const auto& read : reads) {
////        cout << a_read_model.log_probability(read, reference_haplotype, 0) << " "
////            << a_read_model.log_probability(read, true_haplotype, 0) << endl;
////    }
//    
//    Genotype<Haplotype> ref_true {};
//    ref_true.emplace(reference_haplotype);
//    ref_true.emplace(true_haplotype);
//    
//    Genotype<Haplotype> ref_false1 {};
//    ref_false1.emplace(reference_haplotype);
//    ref_false1.emplace(false_haplotype1);
//    
//    Genotype<Haplotype> ref_false2 {};
//    ref_false2.emplace(reference_haplotype);
//    ref_false2.emplace(false_haplotype2);
//    
//    Genotype<Haplotype> true_false1 {};
//    true_false1.emplace(true_haplotype);
//    true_false1.emplace(false_haplotype2);
//    
//    Genotype<Haplotype> true_false2 {};
//    true_false2.emplace(true_haplotype);
//    true_false2.emplace(false_haplotype2);
//    
//    Genotype<Haplotype> false1_false2 {};
//    false1_false2.emplace(false_haplotype1);
//    false1_false2.emplace(false_haplotype2);
//    
////    cout << "ref_true      = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), ref_true, 0) << endl;
////    cout << "ref_false1    = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), ref_false1, 0) << endl;
////    cout << "ref_false2    = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), ref_false2, 0) << endl;
////    cout << "true_false1   = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), true_false1, 0) << endl;
////    cout << "true_false2   = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), true_false2, 0) << endl;
////    cout << "false1_false2 = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), false1_false2, 0) << endl;
//}
//
//BOOST_AUTO_TEST_CASE(obviously_homozygous_sites_evaluate_correctly)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam2});
//    
//    auto a_region = parse_region("11:67503118-67503253", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    auto reads = a_read_manager.fetch_reads(samples[0], a_region);
//    
//    Haplotype reference_haplotype {human, a_region};
//    
//    Haplotype true_haplotype {human, a_region};
//    true_haplotype.push_back(Allele {parse_region("11:67503147-67503148", human), "A"});
//    true_haplotype.push_back(Allele {parse_region("11:67503214-67503215", human), "A"});
//    
//    unsigned ploidy {2};
//    
//    ReadModel a_read_model {ploidy};
//    
//    Genotype<Haplotype> hom_ref {};
//    hom_ref.emplace(reference_haplotype);
//    hom_ref.emplace(reference_haplotype);
//    
//    Genotype<Haplotype> het_alt {};
//    het_alt.emplace(reference_haplotype);
//    het_alt.emplace(true_haplotype);
//    
//    Genotype<Haplotype> hom_alt {};
//    hom_alt.emplace(true_haplotype);
//    hom_alt.emplace(true_haplotype);
//    
////    cout << "genotype likelihoods:" << endl;
////    cout << "hom_ref = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), hom_ref, samples[0]) << endl;
////    cout << "het_alt = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), het_alt, samples[0]) << endl;
////    cout << "hom_alt = " << a_read_model.log_probability(reads.cbegin(), reads.cend(), hom_alt, samples[0]) << endl;
//}
//
//BOOST_AUTO_TEST_CASE(ReadModel_works_on_haploid_genomes)
//{
//    unsigned ploidy {1};
//    
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome ecoli {a_factory.make(ecoli_reference_fasta)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {ecoli_bam});
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(ecoli, 0));
//    
//    auto a_region = parse_region("R00000042:99640-99745", ecoli);
//    
//    auto reference_sequence = ecoli.get_sequence(a_region);
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//    
//    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
//    
//    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
//    
//    auto variants = candidate_generator.get_candidates(a_region);
//    
//    Haplotype reference_haplotype {ecoli, a_region}; // no fully supporting read, just a read with all N's
//    
//    Haplotype best_haplotype {ecoli, a_region}; // all reads fully support this
//    for (const auto& variant : variants) {
//        if (is_snp(variant)) {
//            add_to_back(variant, best_haplotype);
//        }
//    }
//    
//    // Edge case - first snp missing
//    Haplotype missing_snp_haplotype1 {ecoli, a_region};
//    std::for_each(variants.cbegin() + 1, variants.cend(), [&missing_snp_haplotype1] (const auto& v) {
//        if (is_snp(v)) {
//            add_to_back(v, missing_snp_haplotype1);
//        }
//    });
//    
//    // Edge case - last snp missing
//    Haplotype missing_snp_haplotype2 {ecoli, a_region};
//    std::for_each(variants.cbegin(), variants.cend() - 1, [&missing_snp_haplotype2] (const auto& v) {
//        if (is_snp(v)) {
//            add_to_back(v, missing_snp_haplotype2);
//        }
//    });
//    
//    Haplotype okay_haplotype {ecoli, a_region}; // Bad insertion and 3 missing snps
//    add_to_back(variants[0], okay_haplotype);
//    add_to_back(variants[1], okay_haplotype);
//    add_to_back(variants[3], okay_haplotype);
//    add_to_back(variants[4], okay_haplotype);
//    add_to_back(variants[5], okay_haplotype);
//    add_to_back(variants[6], okay_haplotype);
//    add_to_back(variants[11], okay_haplotype);
//    
//    unsigned num_haplotypes {3};
//    std::vector<Haplotype> haplotypes {reference_haplotype, best_haplotype, missing_snp_haplotype1,
//                                        missing_snp_haplotype2, okay_haplotype};
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    ReadModel the_model {ploidy};
//    
//    std::unordered_map<Genotype<Haplotype>, double> genotype_log_probabilities {};
//    
//    for (const auto& genotype : genotypes) {
//        genotype_log_probabilities[genotype] = the_model.log_probability(some_reads.cbegin(), some_reads.cend(), genotype, the_sample_id);
//    }
//    
//    std::sort(genotypes.begin(), genotypes.end(), [&genotype_log_probabilities] (const auto& g1, const auto& g2) {
//        return genotype_log_probabilities[g1] > genotype_log_probabilities[g2];
//    });
//    
////    for (const auto genotype : genotypes) {
////        cout << genotype_log_probabilities.at(genotype) << endl;
////    }
//    
//    BOOST_CHECK(genotypes.at(0).at(0) == best_haplotype);
//    
//    // There is one additional read supporting the last snp so we should expect the haplotype
//    // with this snp missing to be less probable than the one with the first snp missing
//    BOOST_CHECK(genotypes.at(1).at(0) == missing_snp_haplotype1);
//    BOOST_CHECK(genotypes.at(2).at(0) == missing_snp_haplotype2);
//    
//    BOOST_CHECK(genotypes.at(3).at(0) == okay_haplotype);
//    BOOST_CHECK(genotypes.at(4).at(0) == reference_haplotype);
//}
//
//BOOST_AUTO_TEST_CASE(diploid_read_model_test)
//{
//    unsigned ploidy {2};
//    
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    auto a_region = parse_region("2:104142870-104142884", human);
//    
//    auto reference_sequence = human.get_sequence(a_region);
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//    
//    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
//    
//    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
//    
//    auto variants = candidate_generator.get_candidates(a_region);
//    
//    BOOST_CHECK(variants.size() == 3);
//    
//    Haplotype reference_haplotype {human, a_region}; // there are no reads completely supporting the reference
//    
//    Haplotype hap1 {human, a_region};
//    add_to_back(variants[0], hap1); // well supported insert
//    add_to_back(variants[2], hap1); // well supported snp
//    
//    Haplotype hap2 {human, a_region};
//    add_to_back(variants[1], hap2); // this is a low quality snp
//    
//    Haplotype hap3 {human, a_region};
//    add_to_back(variants[0], hap3);
//    add_to_back(variants[1], hap3);
//    add_to_back(variants[2], hap3);
//    
//    unsigned num_haplotypes {4};
//    std::vector<Haplotype> haplotypes {reference_haplotype, hap1, hap2, hap3};
//    
//    ReadModel the_model {ploidy};
//    
////    const auto& hap2_supporting_read =  some_reads.at(1);
////    auto ref_log_prob  = the_model.log_probability(hap2_supporting_read, reference_haplotype, sample);
////    auto hap2_log_prob = the_model.log_probability(hap2_supporting_read, hap2, sample);
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    std::unordered_map<Genotype<Haplotype>, double> genotype_log_probabilities {};
//    
//    for (const auto& genotype : genotypes) {
//        genotype_log_probabilities[genotype] = the_model.log_probability(some_reads.cbegin(), some_reads.cend(), genotype, the_sample_id);
//    }
//    
//    std::sort(genotypes.begin(), genotypes.end(), [&genotype_log_probabilities] (const auto& g1, const auto& g2) {
//        return genotype_log_probabilities[g1] > genotype_log_probabilities[g2];
//    });
//    
////    for (const auto& genotype : genotypes) {
////        std::cout << genotype << " " << genotype_log_probabilities[genotype] << std::endl;
////    }
//    
//    BOOST_CHECK(genotypes.at(0).num_occurences(hap1) == 1);
//    BOOST_CHECK(genotypes.at(0).num_occurences(hap2) == 1);
//    
//    BOOST_CHECK(genotypes.at(1).num_occurences(hap1) == 1);
//    BOOST_CHECK(genotypes.at(1).num_occurences(reference_haplotype) == 1);
//}

//BOOST_AUTO_TEST_CASE(two_sample_diploid_read_model_test)
//{
//    unsigned ploidy {2};
//    
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1, human_1000g_bam2});
//    
//    VariantCandidateGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    auto a_region = parse_region("2:104142870-104142884", human);
//    
//    auto reference_sequence = human.get_sequence(a_region);
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    
//    auto some_reads = a_read_manager.fetch_reads(sample_ids, a_region);
//    
//    candidate_generator.add_reads(some_reads[sample_ids[0]].cbegin(), some_reads[sample_ids[0]].cend());
//    candidate_generator.add_reads(some_reads[sample_ids[1]].cbegin(), some_reads[sample_ids[1]].cend());
//    
//    auto variants = candidate_generator.get_candidates(a_region);
//    
//    BOOST_CHECK(variants.size() == 3);
//    
//    Haplotype reference_haplotype {human, a_region}; // there are no reads completely supporting the reference
//    
//    Haplotype hap1 {human, a_region};
//    add_to_back(variants[0], hap1); // well supported insert
//    add_to_back(variants[2], hap1); // well supported snp
//    
//    Haplotype hap2 {human, a_region};
//    add_to_back(variants[1], hap2); // this is a low quality snp
//    
//    Haplotype hap3 {human, a_region};
//    add_to_back(variants[0], hap3);
//    add_to_back(variants[1], hap3);
//    add_to_back(variants[2], hap3);
//    
//    unsigned num_haplotypes {4};
//    std::vector<Haplotype> haplotypes {reference_haplotype, hap1, hap2, hap3};
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    BOOST_CHECK(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
//    
//    ReadModel a_read_model {ploidy};
//    
//    std::unordered_map<Genotype, double> genotype_log_probabilities0 {};
//    
//    for (const auto& genotype : genotypes) {
//        genotype_log_probabilities0[genotype] = a_read_model.log_probability(some_reads[sample_ids[0]].cbegin(),
//                                                                             some_reads[sample_ids[0]].cend(), genotype, sample_ids[0]);
//    }
//    
//    std::unordered_map<Genotype, double> genotype_log_probabilities1 {};
//    
//    for (const auto& genotype : genotypes) {
//        genotype_log_probabilities1[genotype] = a_read_model.log_probability(some_reads[sample_ids[0]].cbegin(),
//                                                                             some_reads[sample_ids[0]].cend(), genotype, 1);
//    }
//    
////    for (const auto& genotype : genotypes) {
////        std::cout << genotype << " " << genotype_log_probabilities0[genotype] << " " << genotype_log_probabilities1[genotype] << std::endl;
////    }
//}
