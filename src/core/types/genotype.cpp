// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "genotype.hpp"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "utils/maths.hpp"

namespace octopus {

// non-member methods

Genotype<Haplotype> remap(const Genotype<Haplotype>& genotype, const GenomicRegion& region)
{
    Genotype<Haplotype> result {genotype.ploidy()};
    for (const auto& haplotype : genotype) result.emplace(remap(haplotype, region));
    return result;
}

std::size_t num_genotypes(const unsigned num_elements, const unsigned ploidy)
{
    return boost::math::binomial_coefficient<double>(num_elements + ploidy - 1, ploidy);
}

boost::optional<std::size_t> num_genotypes_noexcept(const unsigned num_elements, const unsigned ploidy) noexcept
{
    boost::optional<std::size_t> result {};
    try {
        result = num_genotypes(num_elements, ploidy);
    } catch (...) {}
    return result;
}

std::size_t max_num_elements(const std::size_t num_genotypes, const unsigned ploidy)
{
    if (num_genotypes == 0 || ploidy == 0) return 0;
    auto y = maths::factorial<std::size_t>(ploidy);
    if (y >= num_genotypes) return 1;
    const auto t = num_genotypes * y;
    unsigned j {1};
    for (; j < num_genotypes; ++j) {
        y /= j;
        y *= j + ploidy;
        if (y >= t) break;
    }
    return j + 1;
}

std::size_t element_cardinality_in_genotypes(const unsigned num_elements, const unsigned ploidy)
{
    return ploidy * (num_genotypes(num_elements, ploidy) / num_elements);
}

std::size_t num_max_zygosity_genotypes(const unsigned num_elements, const unsigned ploidy)
{
    namespace bmp = boost::math::policies;
    using policy = bmp::policy<bmp::overflow_error<bmp::throw_on_error>>;
    try {
        return boost::numeric_cast<std::size_t>(boost::math::binomial_coefficient<double>(num_elements, ploidy, policy {}));
    } catch (const boost::numeric::positive_overflow& e) {
        throw std::overflow_error {e.what()};
    }
}

boost::optional<std::size_t> num_max_zygosity_genotypes_noexcept(unsigned num_elements, unsigned ploidy) noexcept
{
    assert(num_elements >= ploidy);
    boost::optional<std::size_t> result {};
    try {
        result = num_max_zygosity_genotypes(num_elements, ploidy);
    } catch (...) {}
    return result;
}

} // namespace octopus
