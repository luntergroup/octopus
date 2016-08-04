//// Copyright (c) 2016 Daniel Cooke
//// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
//
//#define BOOST_TEST_DYN_LINK
//
//#include <boost/test/unit_test.hpp>
//
//#include <iostream>
//#include <string>
//#include <cstddef>
//#include <set>
//#include <vector>
//
//#include "test_common.hpp"
//#include "test_utils.hpp"
//
//#include <io/reference/reference_genome.hpp>
//#include <io/read/read_manager.hpp>
//#include <core/types/variant.hpp>
//#include "composer.hpp"
//#include "cigar_scanner.hpp"
//#include <core/types/haplotype.hpp>
//#include "haplotype_tree.hpp"
//#include <core/types/genotype.hpp>
//
//using std::cout;
//using std::endl;
//
//using test::make_haplotype;
//
//BOOST_AUTO_TEST_SUITE(Components)
//
//BOOST_AUTO_TEST_CASE(can_iterate_Genotype_Haplotypes_with_range_based_for)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele {region, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele {region, "C"}});
//    
//    Genotype<Haplotype> genotype {hap1, hap2};
//    
//    std::vector<Haplotype> r {};
//    
//    for (const auto& haplotype : genotype) {
//        r.push_back(haplotype);
//    }
//    
//    BOOST_CHECK(r.front() == hap1);
//    BOOST_CHECK(r.back() == hap2);
//}
//
//BOOST_AUTO_TEST_CASE(Genotype_can_be_tested_for_haplotype_occurence)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele {region, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele {region, "C"}});
//    const auto hap3 = make_haplotype(human, region, {Allele {region, "G"}});
//    const auto hap4 = make_haplotype(human, region, {Allele {region, "T"}});
//    
//    Genotype<Haplotype> g1 {hap1, hap2, hap3};
//    
//    BOOST_CHECK(g1.contains(hap1));
//    BOOST_CHECK(g1.contains(hap2));
//    BOOST_CHECK(g1.contains(hap3));
//    BOOST_CHECK(!g1.contains(hap4));
//    
//    BOOST_CHECK(g1.count(hap1) == 1);
//    BOOST_CHECK(g1.count(hap2) == 1);
//    BOOST_CHECK(g1.count(hap3) == 1);
//    BOOST_CHECK(g1.count(hap4) == 0);
//    
//    Genotype<Haplotype> g2 {hap1, hap1, hap2};
//    
//    BOOST_CHECK(g2.contains(hap1));
//    BOOST_CHECK(g2.contains(hap2));
//    BOOST_CHECK(!g2.contains(hap3));
//    BOOST_CHECK(!g2.contains(hap4));
//    
//    BOOST_CHECK(g2.count(hap1) == 2);
//    BOOST_CHECK(g2.count(hap2) == 1);
//    BOOST_CHECK(g2.count(hap3) == 0);
//    BOOST_CHECK(g2.count(hap4) == 0);
//    
//    Genotype<Haplotype> g3 {hap1, hap3, hap4};
//    
//    BOOST_CHECK(g3.contains(hap1));
//    BOOST_CHECK(!g3.contains(hap2));
//    BOOST_CHECK(g3.contains(hap3));
//    BOOST_CHECK(g3.contains(hap4));
//    
//    BOOST_CHECK(g3.count(hap1) == 1);
//    BOOST_CHECK(g3.count(hap2) == 0);
//    BOOST_CHECK(g3.count(hap3) == 1);
//    BOOST_CHECK(g3.count(hap4) == 1);
//    
//    Genotype<Haplotype> g4 {hap4, hap4, hap4};
//    
//    BOOST_CHECK(!g4.contains(hap1));
//    BOOST_CHECK(!g4.contains(hap2));
//    BOOST_CHECK(!g4.contains(hap3));
//    BOOST_CHECK(g4.contains(hap4));
//    
//    BOOST_CHECK(g4.count(hap1) == 0);
//    BOOST_CHECK(g4.count(hap2) == 0);
//    BOOST_CHECK(g4.count(hap3) == 0);
//    BOOST_CHECK(g4.count(hap4) == 3);
//}
//
//BOOST_AUTO_TEST_CASE(Genotypes_are_equal_if_they_contain_the_same_haplotypes)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele {region, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele {region, "C"}});
//    const auto hap3 = make_haplotype(human, region, {Allele {region, "G"}});
//    
//    Genotype<Haplotype> g1 {};
//    g1.emplace(hap1);
//    g1.emplace(hap2);
//    g1.emplace(hap3);
//    
//    Genotype<Haplotype> g2 {};
//    g2.emplace(hap1);
//    g2.emplace(hap2);
//    g2.emplace(hap2);
//    
//    Genotype<Haplotype> g3 {};
//    g3.emplace(hap1);
//    g3.emplace(hap2);
//    g3.emplace(hap3);
//    
//    Genotype<Haplotype> g4 {};
//    g4.emplace(hap1);
//    g4.emplace(hap3);
//    g4.emplace(hap3);
//    
//    Genotype<Haplotype> g5 {};
//    g5.emplace(hap1);
//    g5.emplace(hap2);
//    g5.emplace(hap2);
//    
//    BOOST_CHECK(g1 == g1);
//    BOOST_CHECK(g1 != g2);
//    BOOST_CHECK(g1 == g3);
//    BOOST_CHECK(g1 != g4);
//    BOOST_CHECK(g1 != g5);
//    
//    BOOST_CHECK(g2 == g2);
//    BOOST_CHECK(g2 != g3);
//    BOOST_CHECK(g2 != g4);
//    BOOST_CHECK(g2 == g5);
//    
//    BOOST_CHECK(g3 == g3);
//    BOOST_CHECK(g3 != g4);
//    BOOST_CHECK(g3 != g5);
//    
//    BOOST_CHECK(g4 == g4);
//    BOOST_CHECK(g4 != g5);
//    
//    BOOST_CHECK(g5 == g5);
//}
//
//BOOST_AUTO_TEST_CASE(Genotypes_are_not_influenced_by_haplotype_entry_order)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele {region, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele {region, "T"}});
//    
//    Genotype<Haplotype> g1 {};
//    g1.emplace(hap1);
//    g1.emplace(hap2);
//    g1.emplace(hap2);
//    
//    Genotype<Haplotype> g2 {};
//    g2.emplace(hap2);
//    g2.emplace(hap1);
//    g2.emplace(hap2);
//    
//    BOOST_CHECK(g1.count(hap1) == g2.count(hap1));
//    BOOST_CHECK(g1.count(hap2) == g2.count(hap2));
//    
//    BOOST_CHECK(g1 == g2);
//    BOOST_CHECK(std::hash<Genotype<Haplotype>>()(g1) == std::hash<Genotype<Haplotype>>()(g2));
//}
//
//BOOST_AUTO_TEST_CASE(generate_all_genotypes_works_when_num_elements_is_less_than_ploidy)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele {region, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele {region, "T"}});
//    
//    std::vector<Haplotype> haplotypes {hap1};
//    
//    auto genotypes = generate_all_genotypes(haplotypes, 2);
//    
//    BOOST_CHECK(genotypes.size() == 1);
//    
//    genotypes = generate_all_genotypes(haplotypes, 3);
//    
//    BOOST_CHECK(genotypes.size() == 1);
//    
//    genotypes = generate_all_genotypes(haplotypes, 4);
//    
//    BOOST_CHECK(genotypes.size() == 1);
//    
//    genotypes = generate_all_genotypes(haplotypes, 5);
//    
//    BOOST_CHECK(genotypes.size() == 1);
//    
//    haplotypes.push_back(hap2);
//    
//    genotypes = generate_all_genotypes(haplotypes, 3);
//    
//    BOOST_CHECK(genotypes.size() == 4);
//    
//    genotypes = generate_all_genotypes(haplotypes, 4);
//    
//    BOOST_CHECK(genotypes.size() == 5);
//}
//
//BOOST_AUTO_TEST_CASE(generate_all_genotypes_returns_all_possible_unique_genotypes)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele {region, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele {region, "C"}});
//    const auto hap3 = make_haplotype(human, region, {Allele {region, "G"}});
//    const auto hap4 = make_haplotype(human, region, {Allele {region, "T"}});
//    
//    std::vector<Haplotype> haplotypes {hap1, hap2, hap3, hap4};
//    
//    unsigned num_haplotypes {4};
//    
//    auto genotypes_1 = generate_all_genotypes(haplotypes, 1);
//    
//    BOOST_CHECK(genotypes_1.size() == num_genotypes(num_haplotypes, 1));
//    
//    std::unordered_set<Genotype<Haplotype>> unique_1 {genotypes_1.cbegin(), genotypes_1.cend()};
//    
//    BOOST_CHECK(genotypes_1.size() == unique_1.size());
//    
//    auto genotypes_2 = generate_all_genotypes(haplotypes, 2);
//    
//    BOOST_CHECK(genotypes_2.size() == num_genotypes(num_haplotypes, 2));
//    
//    std::unordered_set<Genotype<Haplotype>> unique_2 {genotypes_2.cbegin(), genotypes_2.cend()};
//    
//    BOOST_CHECK(genotypes_2.size() == unique_2.size());
//    
//    auto genotypes_3 = generate_all_genotypes(haplotypes, 3);
//    
//    BOOST_CHECK(genotypes_3.size() == num_genotypes(num_haplotypes, 3));
//    
//    std::unordered_set<Genotype<Haplotype>> unique_3 {genotypes_3.cbegin(), genotypes_3.cend()};
//    
//    BOOST_CHECK(genotypes_3.size() == unique_3.size());
//    
//    auto genotypes_4 = generate_all_genotypes(haplotypes, 4);
//    
//    BOOST_CHECK(genotypes_4.size() == num_genotypes(num_haplotypes, 4));
//    
//    std::unordered_set<Genotype<Haplotype>> unique_4 {genotypes_4.cbegin(), genotypes_4.cend()};
//    
//    BOOST_CHECK(genotypes_4.size() == unique_4.size());
//}
//
//BOOST_AUTO_TEST_CASE(copy_unique_returns_all_the_unique_Haplotypes_in_a_Genotype)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele {region, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele {region, "C"}});
//    const auto hap3 = make_haplotype(human, region, {Allele {region, "G"}});
//    
//    Genotype<Haplotype> g1 {};
//    g1.emplace(hap1);
//    g1.emplace(hap2);
//    g1.emplace(hap3);
//    
//    auto g1_unique = g1.copy_unique();
//    
//    BOOST_CHECK(std::count(g1_unique.cbegin(), g1_unique.cend(), hap1) == 1);
//    BOOST_CHECK(std::count(g1_unique.cbegin(), g1_unique.cend(), hap2) == 1);
//    BOOST_CHECK(std::count(g1_unique.cbegin(), g1_unique.cend(), hap3) == 1);
//    
//    Genotype<Haplotype> g2 {};
//    g2.emplace(hap1);
//    g2.emplace(hap3);
//    g2.emplace(hap3);
//    
//    auto g2_unique = g2.copy_unique();
//    
//    BOOST_CHECK(std::count(g2_unique.cbegin(), g2_unique.cend(), hap1) == 1);
//    BOOST_CHECK(std::count(g2_unique.cbegin(), g2_unique.cend(), hap2) == 0);
//    BOOST_CHECK(std::count(g2_unique.cbegin(), g2_unique.cend(), hap3) == 1);
//    
//    Genotype<Haplotype> g3 {};
//    g3.emplace(hap3);
//    g3.emplace(hap3);
//    g3.emplace(hap3);
//    
//    auto g3_unique = g3.copy_unique();
//    
//    BOOST_CHECK(std::count(g3_unique.cbegin(), g3_unique.cend(), hap1) == 0);
//    BOOST_CHECK(std::count(g3_unique.cbegin(), g3_unique.cend(), hap2) == 0);
//    BOOST_CHECK(std::count(g3_unique.cbegin(), g3_unique.cend(), hap3) == 1);
//}
//
//BOOST_AUTO_TEST_CASE(generate_all_genotypes_results_in_correct_ploidy)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("3:1000000-1000011", human);
//    
//    const auto region1 = parse_region("3:1000000-1000001", human);
//    const auto region2 = parse_region("3:1000010-1000011", human);
//    
//    const auto hap1 = make_haplotype(human, region, {Allele{region1, "A"}, Allele{region2, "A"}});
//    const auto hap2 = make_haplotype(human, region, {Allele{region1, "C"}, Allele{region2, "C"}});
//    const auto hap3 = make_haplotype(human, region, {Allele{region1, "G"}, Allele{region2, "G"}});
//    const auto hap4 = make_haplotype(human, region, {Allele{region1, "A"}, Allele{region2, "C"}});
//    const auto hap5 = make_haplotype(human, region, {Allele{region1, "C"}, Allele{region2, "G"}});
//    const auto hap6 = make_haplotype(human, region, {Allele{region1, "G"}, Allele{region2, "C"}});
//    
//    std::vector<Haplotype> haplotypes {hap1, hap2, hap3, hap4, hap5, hap6};
//    
//    auto genotypes1 = generate_all_genotypes(haplotypes, 1);
//    
//    for (const auto& genotype : genotypes1) {
//        BOOST_CHECK(genotype.ploidy() == 1);
//    }
//    
//    auto genotypes2 = generate_all_genotypes(haplotypes, 2);
//    
//    for (const auto& genotype : genotypes2) {
//        BOOST_CHECK(genotype.ploidy() == 2);
//    }
//    
//    auto genotypes3 = generate_all_genotypes(haplotypes, 3);
//    
//    for (const auto& genotype : genotypes3) {
//        BOOST_CHECK(genotype.ploidy() == 3);
//    }
//}
//
//BOOST_AUTO_TEST_CASE(extract_all_elements_works_correctly)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto variant1 = make_variant("6:31235411-31235412", "A", human);
//    const auto variant2 = make_variant("6:31235412-31235413", "A", human);
//    const auto variant3 = make_variant("6:31235413-31235414", "A", human);
//    const auto variant4 = make_variant("6:31235414-31235415", "A", human);
//    
//    const std::vector<Variant> variants {variant1, variant2, variant3, variant4};
//    
//    auto haplotypes = octopus::generate_all_haplotypes(variants, human);
//    
//    std::sort(std::begin(haplotypes), std::end(haplotypes));
//    
//    const auto haploid_genotypes = generate_all_genotypes(haplotypes, 1);
//    
//    auto haploid_haplotypes = extract_all_elements(haploid_genotypes);
//    
//    std::sort(std::begin(haploid_haplotypes), std::end(haploid_haplotypes));
//    
//    BOOST_CHECK(haplotypes == haploid_haplotypes);
//    
//    const auto diploid_genotypes = generate_all_genotypes(haplotypes, 2);
//    
//    auto diploid_haplotypes = extract_all_elements(diploid_genotypes);
//    
//    std::sort(std::begin(diploid_haplotypes), std::end(diploid_haplotypes));
//    
//    BOOST_CHECK(haplotypes == diploid_haplotypes);
//    
//    const auto triploid_genotypes = generate_all_genotypes(haplotypes, 3);
//    
//    auto triploid_haplotypes = extract_all_elements(triploid_genotypes);
//    
//    std::sort(std::begin(triploid_haplotypes), std::end(triploid_haplotypes));
//    
//    BOOST_CHECK(haplotypes == triploid_haplotypes);
//    
//    const auto tetraploid_genotypes = generate_all_genotypes(haplotypes, 4);
//    
//    auto tetraploid_haplotypes = extract_all_elements(tetraploid_genotypes);
//    
//    std::sort(std::begin(tetraploid_haplotypes), std::end(tetraploid_haplotypes));
//    
//    BOOST_CHECK(haplotypes == tetraploid_haplotypes);
//}
//
//BOOST_AUTO_TEST_CASE(splice_works_correctly)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto variant1 = make_variant("6:31235411-31235412", "A", human);
//    const auto variant2 = make_variant("6:31235412-31235413", "A", human);
//    const auto variant3 = make_variant("6:31235413-31235414", "A", human);
//    const auto variant4 = make_variant("6:31235414-31235414", "A", human);
//    const auto variant5 = make_variant("6:31235414-31235415", "A", human);
//    
//    const std::vector<Variant> variants {variant1, variant2, variant3, variant4, variant5};
//    
//    const auto region = encompassing_region(variants);
//    
//    Haplotype ref_haplotype {region, human};
//    const auto alt_haplotype = make_haplotype(human, region, {
//        variant1.alt_allele(), variant2.alt_allele(), variant3.alt_allele(),
//        variant4.alt_allele(), variant5.alt_allele()
//    });
//    
//    Genotype<Haplotype> genotype {ref_haplotype, alt_haplotype};
//    
//    const auto allele_splice_snp1 = splice<Allele>(genotype, mapped_region(variant1));
//    
//    BOOST_CHECK(!allele_splice_snp1.is_homozygous());
//    BOOST_CHECK(allele_splice_snp1.contains(variant1.ref_allele()));
//    BOOST_CHECK(allele_splice_snp1.contains(variant1.alt_allele()));
//    
//    const auto allele_splice_snp2 = splice<Allele>(genotype, mapped_region(variant2));
//    
//    BOOST_CHECK(!allele_splice_snp2.is_homozygous());
//    BOOST_CHECK(allele_splice_snp2.contains(variant2.ref_allele()));
//    BOOST_CHECK(allele_splice_snp2.contains(variant2.alt_allele()));
//    
//    const auto allele_splice_snp3 = splice<Allele>(genotype, mapped_region(variant3));
//    
//    BOOST_CHECK(!allele_splice_snp3.is_homozygous());
//    BOOST_CHECK(allele_splice_snp3.contains(variant3.ref_allele()));
//    BOOST_CHECK(allele_splice_snp3.contains(variant3.alt_allele()));
//    
//    const auto allele_splice_snp4 = splice<Allele>(genotype, mapped_region(variant5));
//    
//    BOOST_CHECK(!allele_splice_snp4.is_homozygous());
//    BOOST_CHECK(allele_splice_snp4.contains(variant5.ref_allele()));
//    BOOST_CHECK(allele_splice_snp4.contains(variant5.alt_allele()));
//    
//    const auto allele_splice_insertion = splice<Allele>(genotype, mapped_region(variant4));
//    
//    BOOST_CHECK(!allele_splice_insertion.is_homozygous());
//    BOOST_CHECK(allele_splice_insertion.contains(variant4.ref_allele()));
//    BOOST_CHECK(allele_splice_insertion.contains(variant4.alt_allele()));
//}
//
//BOOST_AUTO_TEST_CASE(splice_all_works_correctly)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto variant1 = make_variant("6:31235411-31235412", "A", human);
//    const auto variant2 = make_variant("6:31235412-31235413", "A", human);
//    const auto variant3 = make_variant("6:31235413-31235414", "A", human);
//    const auto variant4 = make_variant("6:31235414-31235415", "A", human);
//    
//    const std::vector<Variant> variants {variant1, variant2, variant3, variant4};
//    
//    const auto haplotypes = octopus::generate_all_haplotypes(variants, human);
//    const auto genotypes  = generate_all_genotypes(haplotypes, 2);
//    
//    const auto allele_splices = splice_all<Allele>(genotypes, mapped_region(variant2));
//    
//    BOOST_REQUIRE(allele_splices.size() == 4);
//    
//    const auto haplotype_splices1 = splice_all<Haplotype>(genotypes, mapped_region(variant2));
//    
//    BOOST_REQUIRE(haplotype_splices1.size() == num_genotypes(2, 2));
//    
//    const auto haplotype_splices2 = splice_all<Haplotype>(genotypes, encompassing_region(variant2, variant3));
//    
//    BOOST_REQUIRE(haplotype_splices2.size() == num_genotypes(4, 2));
//    
//    // TODO
//}
//
//BOOST_AUTO_TEST_CASE(non_member_contains_works_correctly)
//{
//    // TODO
//}
//
//BOOST_AUTO_TEST_CASE(non_member_includes_works_correctly)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto variant1 = make_variant("6:31235411-31235412", "A", human);
//    const auto variant2 = make_variant("6:31235412-31235413", "A", human);
//    const auto variant3 = make_variant("6:31235413-31235414", "A", human);
//    const auto variant4 = make_variant("6:31235414-31235414", "A", human);
//    const auto variant5 = make_variant("6:31235414-31235415", "A", human);
//    
//    const std::vector<Variant> variants {variant1, variant2, variant3, variant4, variant5};
//    
//    const auto region = encompassing_region(variants);
//    
//    Haplotype ref_haplotype {region, human};
//    const auto alt_haplotype = make_haplotype(human, region, {
//        variant1.alt_allele(), variant2.alt_allele(), variant3.alt_allele(),
//        variant4.alt_allele(), variant5.alt_allele()
//    });
//    
//    Genotype<Haplotype> genotype {ref_haplotype, alt_haplotype};
//    
//    BOOST_CHECK(includes(genotype, variant1.alt_allele()));
//    BOOST_CHECK(includes(genotype, variant2.alt_allele()));
//    BOOST_CHECK(includes(genotype, variant3.alt_allele()));
//    BOOST_CHECK(includes(genotype, variant4.alt_allele()));
//    BOOST_CHECK(includes(genotype, variant5.alt_allele()));
//}
//
//BOOST_AUTO_TEST_CASE(are_equal_in_region_works_correctly)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto variant1 = make_variant("6:31235411-31235412", "A", human);
//    const auto variant2 = make_variant("6:31235412-31235413", "",  human);
//    const auto variant3 = make_variant("6:31235413-31235414", "A", human);
//    const auto variant4 = make_variant("6:31235414-31235414", "A", human);
//    const auto variant5 = make_variant("6:31235414-31235415", "A", human);
//    const auto variant6 = make_variant("6:31235415-31235415", "A", human);
//    
//    const std::vector<Variant> variants {variant1, variant2, variant3, variant4, variant5, variant6};
//    
//    const auto region = expand(encompassing_region(variants), 10);
//    
//    const auto ref_haplotype = make_haplotype(human, region, {
//        variant1.ref_allele(), variant2.ref_allele(), variant3.ref_allele(),
//        variant4.ref_allele(), variant5.ref_allele(), variant6.ref_allele()
//    });
//    
//    const auto alt_haplotype = make_haplotype(human, region, {
//        variant1.alt_allele(), variant2.alt_allele(), variant3.alt_allele(),
//        variant4.alt_allele(), variant5.alt_allele(), variant6.alt_allele()
//    });
//    
//    const auto mix_haplotype1 = make_haplotype(human, region, {
//        variant1.alt_allele(), variant2.alt_allele(), variant3.alt_allele(),
//        variant4.ref_allele(), variant5.ref_allele(), variant6.ref_allele()
//    });
//    
//    const auto mix_haplotype2 = make_haplotype(human, region, {
//        variant1.ref_allele(), variant2.ref_allele(), variant3.ref_allele(),
//        variant4.alt_allele(), variant5.alt_allele(), variant6.alt_allele()
//    });
//    
//    Genotype<Haplotype> genotype1 {ref_haplotype, alt_haplotype};
//    Genotype<Haplotype> genotype2 {mix_haplotype1, mix_haplotype2};
//    
////    std::cout << "\n\n";
////    print_alleles(genotype1);
////    std::cout << "\n\n";
////    print_alleles(genotype2);
////    std::cout << "\n\n";
////    print_alleles(splice<Haplotype>(genotype1, mapped_region(variant3)));
////    std::cout << "\n\n";
////    print_alleles(splice<Haplotype>(genotype2, mapped_region(variant3)));
////    std::cout << "\n\n";
////    print_alleles(splice<Haplotype>(genotype1, mapped_region(variant4)));
////    std::cout << "\n\n";
////    print_alleles(splice<Haplotype>(genotype2, mapped_region(variant4)));
////    std::cout << "\n\n";
////    print_alleles(splice<Haplotype>(genotype1, mapped_region(variant5)));
////    std::cout << "\n\n";
////    print_alleles(splice<Haplotype>(genotype2, mapped_region(variant5)));
////    std::cout << "\n\n";
////    
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype1, genotype2, mapped_region(variant1)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype2, genotype1, mapped_region(variant1)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype1, genotype2, mapped_region(variant2)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype2, genotype1, mapped_region(variant2)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype1, genotype2, mapped_region(variant3)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype2, genotype1, mapped_region(variant3)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype1, genotype2, mapped_region(variant4)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype2, genotype1, mapped_region(variant4)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype1, genotype2, mapped_region(variant5)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype2, genotype1, mapped_region(variant5)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype1, genotype2, mapped_region(variant6)));
////    BOOST_CHECK(are_equal_in_region<Haplotype>(genotype2, genotype1, mapped_region(variant6)));
//}
//
//BOOST_AUTO_TEST_SUITE_END()
