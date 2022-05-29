// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

// This thread pool implementation is mostly derived from https://github.com/progschj/ThreadPool

#ifndef thread_pool_hpp
#define thread_pool_hpp

#include <cstddef>
#include <vector>
#include <queue>
#include <functional>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <atomic>
#include <type_traits>
#include <utility>
#include <exception>

namespace octopus {

class ThreadPool
{
public:
    ThreadPool();
    explicit ThreadPool(std::size_t n_threads);
    
    ThreadPool(const ThreadPool&)             = delete;
    ThreadPool& operator=(const ThreadPool&)  = delete;
    ThreadPool(ThreadPool&& other) noexcept   = delete;
    ThreadPool& operator=(ThreadPool&& other) = delete;
    
    ~ThreadPool() noexcept;
    
    std::size_t size() const noexcept;
    bool empty() const noexcept;
    std::size_t n_idle() const noexcept;
    
    void clear() noexcept;
    
    template <typename F, typename... Args>
    auto push(F&& f, Args&&... args) -> std::future<std::result_of_t<F(Args...)>>;
    template <typename F, typename... Args>
    auto try_push(F&& f, Args&&... args) -> std::future<std::result_of_t<F(Args...)>>;
    
private:
    std::mutex mutex_;
    std::condition_variable cv_;
    std::atomic<bool> stop_;
    std::atomic<std::size_t> n_idle_;
    
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
};

template <typename F, typename... Args>
auto ThreadPool::push(F&& f, Args&&... args) -> std::future<std::result_of_t<F(Args...)>>
{
    using f_result_type = std::result_of_t<F(Args...)>;
    auto task = std::make_shared<std::packaged_task<f_result_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    auto result = task->get_future();
    {
        std::lock_guard<std::mutex> lk {mutex_};
        if (stop_) throw std::runtime_error {"ThreadPool: calling push on stopped pool"};
        tasks_.emplace([task] () { (*task)(); });
    }
    cv_.notify_one();
    return result;
}

template <typename F, typename... Args>
auto ThreadPool::try_push(F&& f, Args&&... args) -> std::future<std::result_of_t<F(Args...)>>
{
    using f_result_type = std::result_of_t<F(Args...)>;
    auto task = std::make_shared<std::packaged_task<f_result_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    auto result = task->get_future();
    bool pushed {false};
    {
        std::lock_guard<std::mutex> lk {mutex_};
        if (stop_) throw std::runtime_error {"ThreadPool: calling push on stopped pool"};
        if (n_idle_ > 0) {
            tasks_.emplace([task] () { (*task)(); });
            pushed = true;
        }
    }
    if (pushed) {
        cv_.notify_one();
    } else {
        (*task)(); // run task in calling thread
    }
    return result;
}

} // namespace octopus

#endif
