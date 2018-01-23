// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "thread_pool.hpp"

namespace octopus {

ThreadPool::ThreadPool() : ThreadPool {0} {}

ThreadPool::ThreadPool(const std::size_t n_threads)
: stop_ {false}
, n_idle_ {n_threads}
{
    workers_.reserve(n_threads);
    for (std::size_t i {0}; i < n_threads; ++i) {
        workers_.emplace_back([this] {
            std::function<void()> task;
            while (true) {
                std::unique_lock<std::mutex> lk {mutex_};
                cv_.wait(lk, [this] () { return stop_ || !tasks_.empty(); });
                if (stop_ && tasks_.empty()) return;
                task = std::move(tasks_.front());
                tasks_.pop();
                --n_idle_;
                lk.unlock();
                task();
                lk.lock();
                ++n_idle_;
            }
        });
    }
}

ThreadPool::~ThreadPool() noexcept
{
    {
        std::lock_guard<std::mutex> lk {mutex_};
        stop_ = true;
    }
    cv_.notify_all();
    for (auto& worker : workers_) {
        if (worker.joinable()) worker.join();
    }
}

std::size_t ThreadPool::size() const noexcept
{
    return workers_.size();
}

bool ThreadPool::empty() const noexcept
{
    return workers_.empty();
}

std::size_t ThreadPool::n_idle() const noexcept
{
    return n_idle_;
}

void ThreadPool::clear() noexcept
{
    std::lock_guard<std::mutex> lk {mutex_};
    while (!tasks_.empty()) tasks_.pop();
}

} // namespace octopus
