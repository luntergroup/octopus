//
//  haplotype_tree.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_tree.h"

#include "reference_genome.h"
#include "genotype_model.h"
#include "variant.h"
#include "genomic_region.h"

HaplotypeTree::HaplotypeTree(ReferenceGenome& the_reference, ReadManager& the_reads,
                             GenotypeModel& the_genotype_model,
                             const std::vector<ReadManager::SampleIdType>& the_sample_ids,
                             size_t max_num_haplotypes, double min_posterior)
:
the_reference_ {the_reference},
the_reads_ {the_reads},
the_genotype_model_ {the_genotype_model},
the_sample_ids_ {the_sample_ids},
max_num_haplotypes_ {max_num_haplotypes},
min_posterior_ {min_posterior},
num_extensions_since_posterior_update_ {}
{}

HaplotypeTree::Haplotypes HaplotypeTree::get_haplotypes(const std::vector<Variant>& ordered_variants)
{
    for (const auto& variant : ordered_variants) {
        extend_tree(variant);
        
        if (num_haplotypes() > max_num_haplotypes_) {
            update_posteriors();
            prune_low_probability_haplotypes(num_haplotypes() - max_num_haplotypes_);
        }
    }
    
    return get_highest_probability_haplotypes();
}

size_t HaplotypeTree::num_haplotypes() const
{
    return haplotype_branch_ends_.size();
}

void HaplotypeTree::extend_tree(const Variant& a_variant)
{
    for (auto haplotype_branch : haplotype_branch_ends_) {
        extend_haplotype(haplotype_branch, a_variant);
    }
    ++num_extensions_since_posterior_update_;
}

void HaplotypeTree::extend_haplotype(Vertex haplotype_branch_end, const Variant& a_variant)
{
    auto new_branch_end = boost::add_vertex(the_tree_);
    the_tree_[new_branch_end].the_variant = a_variant;
    the_tree_[new_branch_end].haplotype_probability = the_tree_[new_branch_end].haplotype_probability;
    boost::add_edge(haplotype_branch_end, new_branch_end, the_tree_);
}

HaplotypeTree::Haplotypes HaplotypeTree::get_unupdated_haplotypes()
{
    Haplotypes result {};
    
    for (auto haplotype_end : haplotype_branch_ends_) {
        Haplotype a_haplotype {};
        
        for (size_t i {0}; i < num_extensions_since_posterior_update_; ++i) {
            a_haplotype.emplace_front(the_tree_[haplotype_end].the_variant);
            haplotype_end = *boost::inv_adjacent_vertices(haplotype_end, the_tree_).first;
        }
        
        result.emplace_back(std::move(a_haplotype));
    }
    
    return result;
}

void HaplotypeTree::update_posteriors()
{
    auto unupdated_haplotypes = get_unupdated_haplotypes();
    
    GenomicRegion the_haplotype_region {
        unupdated_haplotypes.front().front().get_reference_allele_region().get_contig_name(),
        unupdated_haplotypes.front().front().get_reference_allele_region().get_begin(),
        unupdated_haplotypes.front().back().get_reference_allele_region().get_end()
    };
    
    auto sample_reads_map = the_reads_.fetch_reads(the_sample_ids_, the_haplotype_region);
    
    GenotypeModel::SampleReads sample_reads {};
    for (const auto& sample_id : the_sample_ids_) {
        sample_reads.emplace_back(std::move(sample_reads_map.at(sample_id)));
    }
    
    auto haplotype_probabilities = the_genotype_model_.get_haplotype_probabilities(unupdated_haplotypes,
                                                                                   sample_reads);
    
    auto it = haplotype_probabilities.begin();
    
    for (auto haplotype : haplotype_branch_ends_) {
        the_tree_[haplotype].haplotype_probability = it->population_probability;
        ++it;
    }
    
    num_extensions_since_posterior_update_ = 0;
}

void HaplotypeTree::prune_low_probability_haplotypes(size_t n)
{
    auto it = haplotype_branch_ends_.begin();
    auto end = haplotype_branch_ends_.end();
    while (n > 0 && it != end) {
        if (the_tree_[*it].haplotype_probability < min_posterior_) {
            prune_haplotype(*it);
            it = haplotype_branch_ends_.erase(it);
        } else {
            ++it;
        }
        --n;
    }
}

void HaplotypeTree::prune_haplotype(Vertex haplotype_end)
{
    while (boost::out_degree(haplotype_end, the_tree_) == 1) {
        auto new_haplotype_end = *boost::inv_adjacent_vertices(haplotype_end, the_tree_).first;
        boost::remove_edge(new_haplotype_end, haplotype_end, the_tree_);
        boost::remove_vertex(haplotype_end, the_tree_);
        haplotype_end = new_haplotype_end;
    }
}

HaplotypeTree::Haplotypes HaplotypeTree::get_highest_probability_haplotypes()
{
    return Haplotypes {};
}
