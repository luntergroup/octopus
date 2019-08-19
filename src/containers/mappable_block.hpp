// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mappable_block_hpp
#define mappable_block_hpp

#include <vector>
#include <initializer_list>
#include <iterator>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"

namespace octopus {

// A MappableBlock is a container of Mappable objects that all have the
// same mapped_region

template <typename MappableTp,
          typename Container = std::vector<MappableTp>>
class MappableBlock : public Mappable<MappableBlock<MappableTp, Container>>,
                      public Comparable<MappableBlock<MappableTp, Container>>
{
public:
    using MappableType = MappableTp;
    using RegionType   = octopus::RegionType<MappableType>;
    
    using allocator_type  = typename Container::allocator_type;
    using value_type      = typename Container::value_type;
    using reference       = typename Container::reference ;
    using const_reference = typename Container::const_reference;
    using difference_type = typename Container::difference_type ;
    using size_type       = typename Container::size_type ;

    using iterator               = typename Container::iterator;
    using const_iterator         = typename Container::const_iterator;
    using reverse_iterator       = typename Container::reverse_iterator;
    using const_reverse_iterator = typename Container::const_reverse_iterator;
    
    struct BadMappableBlock : public std::runtime_error
    {
    public:
        BadMappableBlock(MappableType mappable, RegionType region)
        : runtime_error {"BadMappableBlock"}
        , mappable_ {std::move(mappable)}
        , region_ {std::move(region)}
        {}
        
        virtual ~BadMappableBlock() noexcept = default;
        
        virtual const char* what() const noexcept override { return "bad MappableBlock"; }
        const MappableType& mappable() const noexcept { return mappable_; }
        const RegionType& block_region() const noexcept { return region_; }
    
    private:
        MappableType mappable_;
        RegionType region_;
    };
    
    MappableBlock() = default;
    
    // Empty block (still requires a region)
    MappableBlock(RegionType region)
    : base_ {}
    , region_ {std::move(region)}
    {}
    
    template <typename InputIterator>
    MappableBlock(InputIterator first, InputIterator last)
    : base_ {std::move(first), std::move(last)}
    {
        region_ = octopus::mapped_region(this->front());
        this->check();
    }
    template <typename InputIterator>
    MappableBlock(InputIterator first, InputIterator last, RegionType region)
    : base_ {std::move(first), std::move(last)}
    , region_ {std::move(region)}
    {
        this->check();
    }
    MappableBlock(std::initializer_list<MappableType> values)
    : base_ {values}
    {
        region_ = octopus::mapped_region(this->front());
        this->check();
    }
    MappableBlock(std::initializer_list<MappableType> values, RegionType region)
    : base_ {values}
    , region_ {std::move(region)}
    {
        this->check();
    }
    MappableBlock(const Container& values)
    : base_ {values}
    {
        region_ = octopus::mapped_region(this->front());
        this->check();
    }
    MappableBlock(const Container& values, RegionType region)
    : base_ {values}
    , region_ {std::move(region)}
    {
        this->check();
    }
    MappableBlock(Container&& values)
    : base_ {std::move(values)}
    {
        region_ = octopus::mapped_region(this->front());
        this->check();
    }
    MappableBlock(Container&& values, RegionType region)
    : base_ {std::move(values)}
    , region_ {std::move(region)}
    {
        this->check();
    }
    
    MappableBlock(const MappableBlock&)            = default;
    MappableBlock& operator=(const MappableBlock&) = default;
    MappableBlock(MappableBlock&&)                 = default;
    MappableBlock& operator=(MappableBlock&&)      = default;
    
    MappableBlock& operator=(Container values)
    {
        if (!values.empty()) {
            base_ = std::move(values);
            this->check();
            region_ = octopus::mapped_region(this->front());
        } else {
            this->clear();
        }
        return *this;
    }
    
    ~MappableBlock() = default;
    
    const RegionType& mapped_region() const noexcept { return region_; }
    
    operator Container&() noexcept { return base_; }
    operator const Container&() const noexcept { return base_; }
    
    reference at(size_type pos) { return base_.at(pos); }
    const_reference at(size_type pos) const { return base_.at(pos); }
    reference operator[](size_type pos) noexcept { return base_[pos]; }
    const_reference operator[](size_type pos) const noexcept { return base_[pos]; }
    reference front() noexcept { return base_.front(); }
    const_reference front() const noexcept { return base_.front(); }
    reference back() noexcept { return base_.back(); }
    const_reference back() const noexcept { return base_.back(); }
    
    iterator begin() noexcept { return base_.begin(); }
    const_iterator begin() const noexcept { return base_.begin(); }
    const_iterator cbegin() const noexcept { return base_.cbegin(); }
    iterator end() noexcept { return base_.end(); }
    const_iterator end() const noexcept { return base_.end(); }
    const_iterator cend() const noexcept { return base_.cend(); }
    reverse_iterator rbegin() noexcept { return base_.rbegin(); }
    const_reverse_iterator rbegin() const noexcept { return base_.rbegin(); }
    const_reverse_iterator crbegin() const noexcept { return base_.crbegin(); }
    reverse_iterator rend() noexcept { return base_.rend(); }
    const_reverse_iterator rend() const noexcept { return base_.rend(); }
    const_reverse_iterator crend() const noexcept { return base_.crend(); }
    
    bool empty() const noexcept { return base_.empty(); }
    size_type size() const noexcept { return base_.size(); }
    size_type max_size() const noexcept { return base_.max_size(); }
    void reserve(size_type n) { base_.reserve(n); }
    size_type capacity() const noexcept { return base_.capacity(); }
    void shrink_to_fit() { base_.shrink_to_fit(); }
    
    void clear() noexcept { base_.clear(); }
    iterator insert(const_iterator pos, const MappableType& value) { return base_.insert(std::move(pos), value); }
    iterator insert(const_iterator pos, MappableType&& value) { return base_.insert(std::move(pos), std::move(value)); }
    template <typename InputIterator>
    iterator insert(const_iterator pos, InputIterator first, InputIterator last) { return base_.insert(std::move(pos), std::move(first), std::move(last)); }
    iterator insert(const_iterator pos, std::initializer_list<MappableType> values) { return base_.insert(std::move(pos), values); }
    iterator erase(const_iterator pos) { return base_.erase(pos); }
    iterator erase(const_iterator first, const_iterator last) { return base_.erase(std::move(first), std::move(last)); }
    void push_back(const MappableType& value) { base_.push_back(value); }
    void push_back(MappableType&& value) { base_.push_back(std::move(value)); }
    template <typename ...Args>
    void emplace_back(Args&&... args) { base_.emplace_back(std::forward<Args>(args)...); }
    void pop_back() { base_.pop_back(); }
    void resize(size_type n) { base_.resize(n); }
    void resize(size_type n, const MappableType& value) { base_.resize(n, value); }
    
    void swap(MappableBlock& other) noexcept
    {
        std::swap(this->base_, other.base_);
        std::swap(this->region_, other.region_);
    }
    
    friend bool operator==(const MappableBlock<MappableTp, Container>& lhs, const MappableBlock<MappableTp, Container>& rhs)
    {
        return lhs.region_ == rhs.region_ && lhs.base_ == rhs.base_;
    }
    friend bool operator<(const MappableBlock<MappableTp, Container>& lhs, const MappableBlock<MappableTp, Container>& rhs)
    {
        return lhs.region_ < rhs.region_ || (lhs.region_ == rhs.region_ && lhs.base_ < rhs.base_);
    }
    
private:
    Container base_;
    RegionType region_;
    
    void check() const
    {
        const auto bad_itr = std::find_if(this->cbegin(), this->cend(), [this] (const MappableType& value) { return octopus::mapped_region(value) != region_; });
        if (bad_itr != this->cend()) {
            throw BadMappableBlock {*bad_itr, region_};
        }
    }
};

template <typename MappableTp, typename Container>
void swap(MappableBlock<MappableTp, Container>& lhs, MappableBlock<MappableTp, Container>& rhs) noexcept
{
    lhs.swap(rhs);
}

} // namespace octopus

#endif
