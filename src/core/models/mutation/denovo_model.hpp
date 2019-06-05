// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_model_hpp
#define denovo_model_hpp

#include <cstddef>
#include <unordered_map>
#include <string>
#include <utility>

#include <boost/optional.hpp>
#include <boost/functional/hash.hpp>

#include "core/types/haplotype.hpp"
#include "../pairhmm/pair_hmm.hpp"
#include "indel_mutation_model.hpp"

namespace octopus {

class DeNovoModel
{
public:
    using LogProbability = double;
    
    struct Parameters
    {
        double snv_mutation_rate, indel_mutation_rate;
    };
    
    enum class CachingStrategy { none, value, address };
    
    DeNovoModel() = delete;
    
    DeNovoModel(Parameters parameters,
                std::size_t num_haplotypes_hint = 1000,
                CachingStrategy caching = CachingStrategy::value);
    
    DeNovoModel(const DeNovoModel&)            = default;
    DeNovoModel& operator=(const DeNovoModel&) = default;
    DeNovoModel(DeNovoModel&&)                 = default;
    DeNovoModel& operator=(DeNovoModel&&)      = default;
    
    ~DeNovoModel() = default;
    
    Parameters parameters() const;
    
    void prime(std::vector<Haplotype> haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    // ln p(target | given)
    LogProbability evaluate(const Haplotype& target, const Haplotype& given) const;
    LogProbability evaluate(unsigned target, unsigned given) const;
    
private:
    struct AddressPairHash
    {
        std::size_t operator()(const std::pair<const Haplotype*, const Haplotype*>& p) const noexcept
        {
            auto result = boost::hash_value(p.first);
            boost::hash_combine(result, p.second);
            return result;
        }
    };
    
    using PenaltyVector = hmm::PenaltyVector;
    struct LocalIndelModel
    {
        IndelMutationModel::ContextIndelModel indel;
        PenaltyVector open, extend;
    };
    
    using HMM = hmm::PairHMM<hmm::VariableGapExtendMutationModel, 32, int>;
    
    Parameters params_;
    std::int8_t snv_penalty_;
    IndelMutationModel indel_model_;
    boost::optional<LogProbability> min_ln_probability_;
    std::size_t num_haplotypes_hint_;
    std::vector<Haplotype> haplotypes_;
    CachingStrategy caching_;
    
    mutable hmm::Alignment alignment_;
    mutable LocalIndelModel tmp_indel_model_;
    mutable LocalIndelModel* local_indel_model_;
    mutable std::vector<boost::optional<LocalIndelModel>> gap_model_index_cache_;
    mutable std::unordered_map<Haplotype, std::unordered_map<Haplotype, LogProbability>> value_cache_;
    mutable std::unordered_map<std::pair<const Haplotype*, const Haplotype*>, LogProbability, AddressPairHash> address_cache_;
    mutable std::vector<std::vector<boost::optional<double>>> guarded_index_cache_;
    mutable std::vector<std::vector<LogProbability>> unguarded_index_cache_;
    mutable std::string padded_given_;
    mutable bool use_unguarded_;
    mutable HMM hmm_;
    
    LocalIndelModel generate_local_indel_model(const Haplotype& given) const;
    void set_local_indel_model(unsigned given) const;
    void update_hmm_from_cache() const;
    bool can_try_align_with_hmm(const Haplotype& target, const Haplotype& given) const noexcept;
    void align_with_hmm(const Haplotype& target, const Haplotype& given) const;
    LogProbability evaluate_uncached(const Haplotype& target, const Haplotype& given, bool gap_penalties_cached = false) const;
    LogProbability evaluate_uncached(unsigned target, unsigned given) const;
    LogProbability evaluate_basic_cache(const Haplotype& target, const Haplotype& given) const;
    LogProbability evaluate_address_cache(const Haplotype& target, const Haplotype& given) const;
};

} // namespace octopus
 
#endif
