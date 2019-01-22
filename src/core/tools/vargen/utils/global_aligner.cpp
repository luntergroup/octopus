// Copyright (c) 2015-2018 Daniel Cooke
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
    int match = 0, insertion = 0, deletion = 0;
};

using DPMatrix = std::vector<std::vector<Cell>>;

auto ncols(const DPMatrix& matrix) noexcept
{
    return matrix.size();
}

auto nrows(const DPMatrix& matrix) noexcept
{
    return matrix.front().size();
}

auto init_dp_matrix(const std::string& target, const std::string& query, const Model& model)
{
    assert(!(target.empty() || query.empty()));
    DPMatrix result(target.size() + 1, DPMatrix::value_type(query.size() + 1));
    result[0][0] = Cell {};
    using S = decltype(Cell::match);
    const S min_score {std::min({model.mismatch, model.gap_open})};
    const S inf {min_score * static_cast<S>(std::max({ncols(result), nrows(result)}))};
    for (std::size_t i {1}; i < ncols(result); ++i) {
        result[i][0].match     = inf;
        result[i][0].insertion = inf;
        result[i][0].deletion  = model.gap_open + static_cast<S>(i - 1) * model.gap_extend;
    }
    for (std::size_t j {1}; j < nrows(result); ++j) {
        result[0][j].match     = inf;
        result[0][j].insertion = model.gap_open + static_cast<S>(j - 1) * model.gap_extend;
        result[0][j].deletion  = inf;
    }
    return result;
}

auto match(const std::string& target, const std::string& query,
           const DPMatrix& matrix, const std::size_t i, const std::size_t j,
           const Model& model) noexcept
{
    const auto penalty = target[i - 1] == query[j - 1] ? model.match : model.mismatch;
    const auto& prev = matrix[i - 1][j - 1];
    return std::max({prev.match, prev.insertion, prev.deletion}) + penalty;
}

auto insertion(const DPMatrix& matrix, const std::size_t i, const std::size_t j, const Model& model) noexcept
{
    const auto& prev = matrix[i][j - 1];
    return std::max(prev.insertion + model.gap_extend, prev.match + model.gap_open);
}

auto deletion(const DPMatrix& matrix, const std::size_t i, const std::size_t j, const Model& model) noexcept
{
    const auto& prev = matrix[i - 1][j];
    return std::max(prev.deletion + model.gap_extend, prev.match + model.gap_open);
}

void fill(const std::string& target, const std::string& query,
          DPMatrix& matrix, const std::size_t i, const std::size_t j,
          const Model& model) noexcept
{
    auto& curr = matrix[i][j];
    curr.match     = match(target, query, matrix, i, j, model);
    curr.insertion = insertion(matrix, i, j, model);
    curr.deletion  = deletion(matrix, i, j, model);
}

void fill(DPMatrix& matrix, const std::string& target, const std::string& query, const Model& model) noexcept
{
    for (std::size_t i {1}; i < ncols(matrix); ++i) {
        for (std::size_t j {1}; j < nrows(matrix); ++j) {
            fill(target, query, matrix, i, j, model);
        }
    }
}

auto build_dp_matrix(const std::string& target, const std::string& query, const Model& model)
{
    auto result = init_dp_matrix(target, query, model);
    fill(result, target, query, model);
    return result;
}

char traceback(const std::string& target, const std::string& query, const DPMatrix& matrix,
               const std::size_t i, const std::size_t j, const Model& model, const char prev_state) noexcept
{
    const auto& curr = matrix[i][j];
    if (prev_state == '$') {
        if (curr.match >= curr.deletion) {
            if (curr.match >= curr.insertion) {
                return target[i - 1] == query[j - 1] ? '=' : 'X';
            } else {
                return 'I';
            }
        } else {
            return curr.deletion >= curr.insertion ? 'D' : 'I';
        }
    } else {
        if (prev_state == '=' || prev_state == 'X') {
            if (curr.match >= curr.deletion) {
                if (curr.match >= curr.insertion) {
                    return target[i - 1] == query[j - 1] ? '=' : 'X';
                } else {
                    return 'I';
                }
            } else {
                return curr.deletion >= curr.insertion ? 'D' : 'I';
            }
        } else if (prev_state == 'D') {
            const auto& prev = matrix[i + 1][j];
            if (prev.deletion == curr.match + model.gap_open) {
                return target[i - 1] == query[j - 1] ? '=' : 'X';
            } else {
                return 'D';
            }
        } else {
            assert(prev_state == 'I');
            const auto& prev = matrix[i][j + 1];
            if (prev.insertion == curr.match + model.gap_open) {
                return target[i - 1] == query[j - 1] ? '=' : 'X';
            } else {
                return 'I';
            }
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

auto extract_alignment(const std::string& target, const std::string& query, const DPMatrix& matrix, const Model& model)
{
    AlignmentString alignment {};
    auto i = ncols(matrix) - 1;
    auto j = nrows(matrix) - 1;
    char state {'$'};
    while(i > 0 || j > 0) {
        using Flag = CigarOperation::Flag;
        state = traceback(target, query, matrix, i, j, model, state);
        switch(state) {
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

auto score(const DPMatrix& matrix) noexcept
{
    const auto& last = matrix.back().back();
    return std::max({last.match, last.insertion, last.deletion});
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
    const auto matrix = build_dp_matrix(target, query, model);
    return {extract_alignment(target, query, matrix, model), score(matrix)};
}

} // namespace coretools
} // namespace octopus
