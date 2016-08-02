//
//  matrix_map.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef matrix_map_hpp
#define matrix_map_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <cstddef>
#include <utility>
#include <iterator>
#include <type_traits>
#include <algorithm>
#include <memory>

template <
typename Key1,
typename Key2,
typename T,
typename Hash1 = std::hash<Key1>,
typename Hash2 = std::hash<Key2>,
typename KeyEqual1 = std::equal_to<Key1>,
typename KeyEqual2 = std::equal_to<Key2>
> class MatrixMap
{
    using Key2ContainerType = std::vector<Key2>;
    using ValueContainerType = std::vector<T>;
    
    struct Key2RefHash
    {
        std::size_t operator()(std::reference_wrapper<const Key2> k) const
        {
            return Hash2()(k.get());
        }
    };
    
    struct Key2RefEqual
    {
        bool operator()(const Key2& lhs, std::reference_wrapper<const Key2> rhs) const
        {
            return KeyEqual2()(lhs, rhs.get());
        }
    };
    
    using ValueMap = std::unordered_map<Key1, ValueContainerType, Hash1, KeyEqual1>;
    using IndexSizeType = typename ValueMap::mapped_type::size_type;
    using IndiceMap = std::unordered_map<std::reference_wrapper<const Key2>, IndexSizeType, Key2RefHash, Key2RefEqual>;
    
    using Key2Iterator  = typename Key2ContainerType::const_iterator;
    using ValueIterator = typename ValueContainerType::const_iterator;
    
public:
    using size_type = std::size_t;
    
    class InnerMap;
    class InnerSlice;
    
    MatrixMap()  = default;
    ~MatrixMap() = default;
    
    template <typename InputIt>
    MatrixMap(InputIt first, InputIt last) : key2s_ {first, last}
    {
        this->generate_indice_map();
    }
    
    MatrixMap(const MatrixMap& other)
    :
    key2s_ {other.key2s_},
    values_ {other.values_}
    {
        this->generate_indice_map();
    }
    
    MatrixMap& operator=(const MatrixMap& other)
    {
        if (&other == this) {
            return *this;
        }
        
        key2s_  = other.key2s_;
        values_ = other.values_;
        
        this->regenerate_indice_map();
        
        return *this;
    }
    
    MatrixMap(MatrixMap&& other)
    :
    key2s_ {std::move(other.key2s_)},
    values_ {std::move(other.values_)}
    {
        this->generate_indice_map();
    }
    
    MatrixMap& operator=(MatrixMap&& other)
    {
        if (&other == this) {
            return *this;
        }
        
        key2s_  = std::move(other.key2s_);
        values_ = std::move(other.values_);
        
        this->regenerate_indice_map();
        
        return *this;
    }
    
    T& operator()(const Key1& key1, const Key2& key2)
    {
        return values_.at(key1)[key2_indices_.at(key2)];
    }
    
    const T& operator()(const Key1& key1, const Key2& key2) const
    {
        return values_.at(key1)[key2_indices_.at(key2)];
    }
    
    InnerSlice operator()(const Key1& key) const
    {
        const auto& key_values = values_.at(key);
        return InnerSlice {std::cbegin(key_values), std::cend(key_values)};
    }
    
    InnerMap operator[](const Key1& key) const
    {
        return InnerMap {this->begin(key), this->end(key), key2_indices_};
    }
    
    bool empty1() const noexcept
    {
        return values_.empty();
    }
    
    bool empty2() const noexcept
    {
        return key2s_.empty();
    }
    
    size_type size1() const noexcept
    {
        return values_.size();
    }
    
    size_type size2() const noexcept
    {
        return key2s_.size();
    }
    
    void reserve1(size_type n)
    {
        values_.reserve(n);
    }
    
    void reserve2(size_type n)
    {
        key2s_.reserve(n);
        key2_indices_.reserve(n);
    }
    
    void reserve(size_type n1, size_type n2)
    {
        reserve1(n1);
        reserve2(n2);
    }
    
    void clear() noexcept
    {
        key2s_.clear();
        values_.clear();
        key2_indices_.clear();
    }
    
    template <typename InputIt>
    bool assign_keys(InputIt first, InputIt last)
    {
        key2s_.assign(first, last);
        
        this->regenerate_indice_map();
        
        if (!values_.empty()) {
            values_.clear();
            return true;
        }
        
        return false;;
    }
    
    template <typename K>
    bool push_back(K&& key)
    {
        this->push_back_reallocate(std::forward<K>(key));
        
        if (!values_.empty()) {
            values_.clear();
            return true;
        }
        
        return false;
    }
    
    template <typename... Args>
    bool emplace_back(Args&&... args)
    {
        this->emplace_back_reallocate(std::forward<Args>(args)...);
        
        if (!values_.empty()) {
            values_.clear();
            return true;
        }
        
        return false;
    }
    
    template <typename K, typename InputIt>
    bool insert_at(K&& key, InputIt first, InputIt last)
    {
        if (std::distance(first, last) != this->size2()) {
            throw std::out_of_range {"MatrixMap::insert_at called with value range of different"
                " length to Key2 range in this MatrixMap"};
        }
        
        return values_.emplace(std::piecewise_construct,
                               std::forward_as_tuple<K>(std::forward<K>(key)),
                               std::make_tuple(first, last)).second;
    }
    
    template <typename K, typename InputIt>
    bool insert_or_assign_at(K&& key, InputIt first, InputIt last)
    {
        if (values_.count(key) == 0) {
            return insert_at(std::forward<K>(key), first, last);
        }
        
        if (std::distance(first, last) != this->size2()) {
            throw std::out_of_range {"MatrixMap::insert_at called with value range of different"
                " length to Key2 range in this MatrixMap"};
        }
        
        values_[std::forward<K>(key)].assign(first, last);
        
        return false;
    }
    
    template <typename K, typename InputIt>
    bool insert_each(K&& key, InputIt first, InputIt last)
    {
        if (std::distance(first, last) != this->size1()) {
            throw std::out_of_range {"MatrixMap::insert_each called with value range of different"
                " length to Key1 range in this MatrixMap"};
        }
        
        // TODO: not working for all forwards
        //        if (key2_indices_.count(key) != 0) {
        //            return false;
        //        }
        
        for (auto& p : values_) {
            p.second.push_back(*first++);
        }
        
        this->push_back_reallocate(std::forward<K>(key));
        
        return true;
    }
    
    template <typename K, typename InputIt>
    bool insert_or_assign_each(K&& key, InputIt first, InputIt last)
    {
        if (std::distance(first, last) != this->size1()) {
            throw std::out_of_range {"MatrixMap::insert_each called with value range of different"
                " length to Key1 range in this MatrixMap"};
        }
        
        if (key2_indices_.count(key) == 0) {
            for (auto& p : values_) {
                p.second.push_back(*first++);
            }
            
            this->push_back_reallocate(std::forward<K>(key));
            
            return true;
        }
        
        const auto index = key2_indices_[key];
        
        for (auto& p : values_) {
            p.second[index] = *first++;
        }
        
        this->push_back_reallocate(std::forward<K>(key));
        
        return false;
    }
    
    void pop_back()
    {
        if (key2s_.empty()) {
            return;
        }
        
        key2_indices_.erase(key2s_.back());
        
        key2s_.pop_back();
        
        for (auto& p : values_) {
            p.second.pop_back();
        }
    }
    
    bool erase1(const Key1& key)
    {
        return values_.erase(key) == 1;
    }
    
    bool erase2(const Key2& key)
    {
        if (key2_indices_.count(key) == 0) {
            return false;
        }
        
        const auto key_index = key2_indices_[key];
        
        for (auto& p : values_) {
            p.second.erase(std::next(std::begin(p.second), key_index));
        }
        
        const auto it = key2s_.erase(std::next(std::begin(key2s_), key_index));
        
        key2_indices_.erase(key);
        
        std::for_each(it, std::end(key2s_), [this] (const auto& key) { --key2_indices_[key]; });
        
        return true;
    }
    
    // TODO
    template <typename InputIt>
    size_type remove_keys(InputIt first, InputIt last)
    {
        const auto key2_begin = std::begin(key2s_);
        const auto key2_end   = std::end(key2s_);
        
        const auto it = std::find_first_of(key2_begin, key2_end, first, last);
        
        if (it == key2_end) return 0;
        
        std::vector<IndexSizeType> indices_removed {};
        indices_removed.reserve(std::distance(it, key2_end));
        
        auto i = std::distance(key2_begin, it);
        
        indices_removed.push_back(i++);
        
        for (auto it2 = it; ++it != key2_end; ++i) {
            if (std::find(first, last, *it) == last) {
                *it++ = std::move(*it2);
            } else {
                indices_removed.push_back(i);
            }
        }
        
        indices_removed.shrink_to_fit();
        
        return 1;
    }
    
private:
    Key2ContainerType key2s_;
    ValueMap values_;
    IndiceMap key2_indices_;
    
    void generate_indice_map()
    {
        key2_indices_.reserve(key2s_.size());
        
        IndexSizeType i {0};
        
        for (const auto& key : key2s_) {
            key2_indices_.emplace(key, i++);
        }
    }
    
    void regenerate_indice_map()
    {
        key2_indices_.clear();
        generate_indice_map();
    }
    
    template <typename K>
    void push_back_reallocate(K&& key)
    {
        if (key2s_.capacity() == key2s_.size()) {
            key2s_.push_back(std::forward<K>(key));
            this->regenerate_indice_map();
        } else {
            key2s_.push_back(std::forward<K>(key));
            key2_indices_.emplace(key2s_.back(), key2s_.size() - 1);
        }
    }
    
    template <typename... Args>
    void emplace_back_reallocate(Args&&... args)
    {
        if (key2s_.capacity() == key2s_.size()) {
            key2s_.emplace_back(std::forward<Args>(args)...);
            this->regenerate_indice_map();
        } else {
            key2s_.emplace_back(std::forward<Args>(args)...);
            key2_indices_.emplace(key2s_.back(), key2s_.size() - 1);
        }
    }
    
public:
    class InnerSlice
    {
    public:
        InnerSlice() = delete;
        ~InnerSlice() = default;
        
        InnerSlice(ValueIterator begin, ValueIterator end) : begin_ {begin}, end_ {end} {}
        
        ValueIterator begin() const { return begin_; }
        ValueIterator end() const { return end_; }
        ValueIterator cbegin() const { return begin_; }
        ValueIterator cend() const { return end_; }
    private:
        ValueIterator begin_, end_;
    };
    
    class ZipIterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = std::pair<const Key2&, const T&>;
        using difference_type   = typename std::iterator_traits<Key2Iterator>::difference_type;
        using reference         = value_type;
        using pointer           = const value_type*;
        
        ZipIterator() = delete;
        
        explicit ZipIterator(Key2Iterator key2_itr, ValueIterator value_itr)
        : key2_itr_ {key2_itr}, value_itr_ {value_itr} {}
        
        ~ZipIterator() = default;
        
        ZipIterator& operator++()
        {
            ++key2_itr_;
            ++value_itr_;
            return *this;
        }
        
        value_type operator*() const
        {
            return std::make_pair(std::ref(*key2_itr_), std::ref(*value_itr_));
        }
        
        auto operator->() const
        {
            return std::make_unique<value_type>(*key2_itr_, *value_itr_);
        }
        
        friend bool operator==(const ZipIterator& lhs, const ZipIterator& rhs)
        {
            return lhs.key2_itr_ == rhs.key2_itr_;
        }
        
        friend bool operator!=(const ZipIterator& lhs, const ZipIterator& rhs)
        {
            return !operator==(lhs, rhs);
        }
        
        friend InnerMap;
        
    private:
        Key2Iterator key2_itr_;
        ValueIterator value_itr_;
    };
    
    ZipIterator begin(const Key1& key) const
    {
        return ZipIterator {std::begin(key2s_), std::begin(values_.at(key))};
    }
    
    ZipIterator end(const Key1& key) const
    {
        return ZipIterator {std::end(key2s_), std::end(values_.at(key))};
    }
    
    ZipIterator cbegin(const Key1& key) const
    {
        return begin(key);
    }
    
    ZipIterator cend(const Key1& key) const
    {
        return end(key);
    }
    
private:
    using ValueMapIterator = typename ValueMap::const_iterator;
    
public:
    class InnerMap
    {
    public:
        using key_type = Key2;
        using mapped_type = T;
        using value_type = std::pair<const key_type, T>;
        
        InnerMap() = delete;
        ~InnerMap() = default;
        
        InnerMap(ZipIterator begin, ZipIterator end, const IndiceMap& key2_indices)
        : begin_ {begin}, end_ {end}, key2_indices_ {key2_indices} {}
        
        ZipIterator begin() const { return begin_; }
        ZipIterator end() const { return end_; }
        ZipIterator cbegin() const { return begin_; }
        ZipIterator cend() const { return end_; }
        
        IndexSizeType size() const noexcept { return key2_indices_.get().size(); }
        
        //        T& operator[](const Key2& key)
        //        {
        //            return *std::next(begin_.value_itr_, key2_indices_.get().at(key));
        //        }
        
        const T& operator[](const Key2& key) const
        {
            return *std::next(begin_.value_itr_, key2_indices_.get().at(key));
        }
        
        const T& at(const Key2& key) const
        {
            return *std::next(begin_.value_itr_, key2_indices_.get().at(key));
        }
        
    private:
        ZipIterator begin_, end_;
        std::reference_wrapper<const IndiceMap> key2_indices_;
    };
    
    class Iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = std::pair<const Key1&, InnerMap>;
        using difference_type   = typename std::iterator_traits<Key2Iterator>::difference_type;
        using reference         = value_type;
        using pointer           = value_type*;
        
        Iterator() = delete;
        
        explicit Iterator(ValueMapIterator map_itr, Key2Iterator key2_begin_itr,
                          const IndiceMap& key2_indices)
        : map_itr_ {map_itr}, key2_begin_itr_ {key2_begin_itr}, key2_indices_ {key2_indices} {}
        
        ~Iterator() = default;
        
        Iterator& operator++()
        {
            map_itr_++;
            return *this;
        }
        
        value_type operator*() const
        {
            const auto zip1 = ZipIterator {key2_begin_itr_, std::begin(map_itr_->second)};
            const auto zip2 = ZipIterator {std::next(key2_begin_itr_, key2_indices_.get().size()),
                std::end(map_itr_->second)};
            return std::make_pair(std::ref(map_itr_->first), InnerMap {zip1, zip2, key2_indices_});
        }
        
        auto operator->() const
        {
            const auto zip1 = ZipIterator {key2_begin_itr_, std::begin(map_itr_->second)};
            const auto zip2 = ZipIterator {std::next(key2_begin_itr_, key2_indices_.get().size()),
                std::end(map_itr_->second)};
            return std::make_unique<value_type>(map_itr_->first, InnerMap {zip1, zip2, key2_indices_});
        }
        
        friend bool operator==(const Iterator& lhs, const Iterator& rhs)
        {
            return lhs.map_itr_ == rhs.map_itr_;
        }
        
        friend bool operator!=(const Iterator& lhs, const Iterator& rhs)
        {
            return !operator==(lhs, rhs);
        }
        
    private:
        ValueMapIterator map_itr_;
        Key2Iterator key2_begin_itr_;
        std::reference_wrapper<const IndiceMap> key2_indices_;
    };
    
    Iterator begin() const { return Iterator {std::cbegin(values_), std::begin(key2s_), key2_indices_}; }
    Iterator end() const { return Iterator {std::cend(values_), std::begin(key2s_), key2_indices_}; }
    Iterator cbegin() const { return begin(); }
    Iterator cend() const { return end(); }
};

#endif /* matrix_map_hpp */
