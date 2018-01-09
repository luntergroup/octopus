// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "global_aligner.hpp"

#include <vector>
#include <deque>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <functional>
#include <cassert>

namespace octopus { namespace coretools {

namespace {

struct Cell
{
    int score      = 0;
    char traceback = '$';
};

bool operator<(const Cell& lhs, const Cell& rhs) noexcept
{
    return lhs.score < rhs.score;
}

using DPMatrix = std::vector<std::vector<Cell>>;

auto ncols(const DPMatrix& matrix) noexcept
{
    return matrix.size();
}

auto nrows(const DPMatrix& matrix) noexcept
{
    return matrix.front().size();
}

auto init_dp_matrix(const std::string& target, const std::string& query,
                    const Model& model)
{
    assert(!(target.empty() || query.empty()));
    DPMatrix result(target.size() + 1, DPMatrix::value_type(query.size() + 1));
    result[0][0] = Cell {};
    using S = decltype(Cell::score);
    for (std::size_t i {1}; i < ncols(result); ++i) {
        result[i][0].score = model.gap_open + static_cast<S>(i - 1) * model.gap_extend;
        result[i][0].traceback = 'D';
    }
    for (std::size_t j {1}; j < nrows(result); ++j) {
        result[0][j].score = model.gap_open + static_cast<S>(j - 1) * model.gap_extend;
        result[0][j].traceback = 'I';
    }
    return result;
}

auto match(const std::string& target, const std::string& query,
           const DPMatrix& matrix, const std::size_t i, const std::size_t j,
           const Model& model) noexcept
{
    if (target[i - 1] == query[j - 1]) {
        return Cell {matrix[i - 1][j - 1].score + model.match, '='};
    } else {
        return Cell {matrix[i - 1][j - 1].score + model.mismatch, 'X'};
    }
}

auto insertion(const DPMatrix& matrix, const std::size_t i, const std::size_t j,
               const Model& model) noexcept
{
    const auto score = matrix[i][j - 1].traceback == 'I' ? model.gap_extend : model.gap_open;
    return Cell {matrix[i][j - 1].score + score, 'I'};
}

auto deletion(const DPMatrix& matrix, const std::size_t i, const std::size_t j,
              const Model& model) noexcept
{
    const auto score = matrix[i - 1][j].traceback == 'D' ? model.gap_extend : model.gap_open;
    return Cell {matrix[i - 1][j].score + score, 'D'};
}

Cell extract_max(const std::string& target, const std::string& query,
                 const DPMatrix& matrix, const std::size_t i, const std::size_t j,
                 const Model& model) noexcept
{
    return std::max({
        insertion(matrix, i, j, model),
        deletion(matrix, i, j, model),
        match(target, query, matrix, i, j, model)
    });
}

void fill(DPMatrix& matrix, const std::string& target, const std::string& query,
          const Model& model) noexcept
{
    for (std::size_t i {1}; i < ncols(matrix); ++i) {
        for (std::size_t j {1}; j < nrows(matrix); ++j) {
            matrix[i][j] = extract_max(target, query, matrix, i, j, model);
        }
    }
}

using AlignmentString = std::deque<CigarOperation::Flag>;

auto make_cigar(const AlignmentString& alignment)
{
    CigarString result {};
    result.reserve(alignment.size());
    auto itr = std::cbegin(alignment);
    const auto last = std::cend(alignment);
    while (itr != last) {
        auto next_unique = std::adjacent_find(itr, last, std::not_equal_to<> {});
        if (next_unique != last) ++next_unique;
        assert(std::distance(itr, next_unique) > 0);
        result.emplace_back(static_cast<CigarOperation::Size>(std::distance(itr, next_unique)), *itr);
        itr = next_unique;
    }
    result.shrink_to_fit();
    return result;
}

auto extract_alignment(const DPMatrix& matrix)
{
    AlignmentString alignment {};
    auto i = ncols(matrix) - 1;
    auto j = nrows(matrix) - 1;
    while(i > 0 || j > 0) {
        using Flag = CigarOperation::Flag;
        switch(matrix[i][j].traceback) {
            case '=':
            {
                assert(i > 0 && j > 0);
                alignment.push_front(Flag::sequenceMatch);
                --i;
                --j;
                break;
            }
            case 'X':
            {
                assert(i > 0 && j > 0);
                alignment.push_front(Flag::substitution);
                --i;
                --j;
                break;
            }
            case 'I':
            {
                assert(j > 0);
                alignment.push_front(Flag::insertion);
                --j;
                break;
            }
            case 'D':
            {
                assert(i > 0);
                alignment.push_front(Flag::deletion);
                --i;
                break;
            }
        }
    }
    return make_cigar(alignment);
}

} // namespace

Alignment align(const std::string& target, const std::string& query, Model model)
{
    using Flag = CigarOperation::Flag;
    using Size = CigarOperation::Size;
    if (target.empty()) {
        if (query.empty()) return {CigarString {}, 0};
        return {CigarString {CigarOperation {static_cast<Size>(query.size()), Flag::insertion}},
                model.gap_open + static_cast<int>(query.size() - 1) * model.gap_extend};
    }
    if (query.empty()) {
        return {CigarString {CigarOperation {static_cast<Size>(target.size()), Flag::deletion}},
                model.gap_open + static_cast<int>(target.size() - 1) * model.gap_extend};
    }
    auto matrix = init_dp_matrix(target, query, model);
    fill(matrix, target, query, model);
    return {extract_alignment(matrix), matrix.back().back().score};
}

} // namespace coretools
} // namespace octopus
