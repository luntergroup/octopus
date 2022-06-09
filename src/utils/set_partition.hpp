// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <vector>
#include <cstddef>

// Based on algorithm presented in https://codereview.stackexchange.com/a/1944/52596

namespace octopus {

namespace detail {

template <typename UnaryOperation>
void visit(const int m, UnaryOperation op,
           int n, std::vector<int>& a)
{
    std::vector<std::vector<std::size_t>> partitions(m);
    for (int j {0}; j < n; ++j) {
        partitions[a[j + 1]].push_back(j);
    }
    op(partitions);
}

template <typename UnaryOperation>
void b(const int m, UnaryOperation op,
       int mu, int nu, int sigma, int n,
       std::vector<int>& a);

template <typename UnaryOperation>
void f(const int m, UnaryOperation op,
       int mu, int nu, int sigma, int n,
       std::vector<int>& a)
{
    if (mu == 2) {
        visit(m, op, n, a);
    } else {
        f(m, op, mu - 1, nu - 1, (mu + sigma) % 2, n, a);
    }
    if (nu == mu + 1) {
        a[mu] = mu - 1;
        visit(m, op, n, a);
        while (a[nu] > 0) {
            --a[nu];
            visit(m, op, n, a);
        }
    } else if (nu > mu + 1) {
        if ((mu + sigma) % 2 == 1) {
            a[nu - 1] = mu - 1;
        } else {
            a[mu] = mu - 1;
        }
        if ((a[nu] + sigma) % 2 == 1) {
            b(m, op, mu, nu - 1, 0, n, a);
        } else {
            f(m, op, mu, nu - 1, 0, n, a);
        }
        while (a[nu] > 0) {
            --a[nu];
            if ((a[nu] + sigma) % 2 == 1) {
                b(m, op, mu, nu - 1, 0, n, a);
            } else {
                f(m, op, mu, nu - 1, 0, n, a);
            }
        }
    }
}

template <typename UnaryOperation>
void b(const int m, UnaryOperation op,
       int mu, int nu, int sigma, int n,
       std::vector<int>& a)
{
    if (nu == mu + 1) {
        while (a[nu] < mu - 1) {
            visit(m, op, n, a);
            ++a[nu];
        }
        visit(m, op, n, a);
        a[mu] = 0;
    } else if (nu > mu + 1) {
        if ((a[nu] + sigma) % 2 == 1) {
            f(m, op, mu, nu - 1, 0, n, a);
        } else {
            b(m, op, mu, nu - 1, 0, n, a);
        }
        while (a[nu] < mu - 1) {
            ++a[nu];
            if ((a[nu] + sigma) % 2 == 1) {
                f(m, op, mu, nu - 1, 0, n, a);
            } else {
                b(m, op, mu, nu - 1, 0, n, a);
            }
        }
        if ((mu + sigma) % 2 == 1) {
            a[nu - 1] = 0;
        } else {
            a[mu] = 0;
        }
    }
    if (mu == 2) {
        visit(m, op, n, a);
    } else {
        b(m, op, mu - 1, nu - 1, (mu + sigma) % 2, n, a);
    }
}

} // namespace detail

template <typename UnaryOperation>
void set_partitions(int n, int m, UnaryOperation op)
{
    std::vector<int> a(n + 1);
    detail::visit(m, op, n, a);
    for (int j {1}; j < m + 1; ++j) {
        a[n - m + j] = j - 1;
    }
    detail::f(m, op, m, n, 0, n, a);
}

} // namespace octopus
