//
//  indel_prior_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

//void biggest_repeats_in_huamn()
//{
//    auto reference = make_reference(human_reference_fasta);
//    
//    auto contigs = reference.get_contig_names();
//    
//    std::vector<GenomicRegion> biggest_regions {};
//    
//    auto start = std::chrono::system_clock::now();
//    
//    for (const auto& contig : contigs) {
//        cout << "finding repeats in chromosome " << contig << endl;
//        
//        auto region = reference.get_contig_region(contig);
//        auto repeats = find_exact_tandem_repeats(reference.get_sequence(region), region, 100);
//        
//        if (!repeats.empty()) {
//            auto it = std::max_element(repeats.cbegin(), repeats.cend(), [] (const auto& lhs, const auto& rhs) { return size(lhs.region) < size(rhs.region); });
//            biggest_regions.push_back(it->region);
//        }
//    }
//    
//    std::sort(biggest_regions.begin(), biggest_regions.end(), [] (const auto& lhs, const auto& rhs) { return size(lhs) > size(rhs); });
//    
//    cout << "biggest exact repeats in each chromosome" << endl;
//    for (auto& region : biggest_regions) {
//        cout << region << " " << size(region) << endl;
//    }
//    
//    auto end = std::chrono::system_clock::now();
//    
//    auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();
//    
//    cout << "took " << duration << " minutes" << endl;
//}
