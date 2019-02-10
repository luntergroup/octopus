// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef emplace_iterator_hpp
#define emplace_iterator_hpp

#include <iterator>
#include <utility>

namespace octopus {
    
template <typename Container>
class front_emplace_iterator : public std::iterator<std::output_iterator_tag, void, void, void, void>
{
protected:
    Container* container;
public:
    using container_type = Container;
    
    front_emplace_iterator() = default;
    
    explicit front_emplace_iterator(Container& x) noexcept : container(&x) {}
    
    front_emplace_iterator(const front_emplace_iterator&)            = default;
    front_emplace_iterator& operator=(const front_emplace_iterator&) = default;
    front_emplace_iterator(front_emplace_iterator&&)                 = default;
    front_emplace_iterator& operator=(front_emplace_iterator&&)      = default;
    
    ~front_emplace_iterator() = default;

    template <typename T>
    front_emplace_iterator<Container>& operator=(T&& t)
    {
        container->emplace_front(std::forward<T>(t));
        return *this;
    }
    
    front_emplace_iterator& operator*()     { return *this; }
    front_emplace_iterator& operator++()    { return *this; }
    front_emplace_iterator& operator++(int) { return *this; }
};

template <typename Container>
inline front_emplace_iterator<Container> front_emplacer(Container& c) noexcept
{
    return front_emplace_iterator<Container>(c);
}

template <typename Container>
class back_emplace_iterator : public std::iterator<std::output_iterator_tag, void, void, void, void>
{
protected:
    Container* container;
public:
    using container_type = Container;
    
    back_emplace_iterator() = default;
    
    explicit back_emplace_iterator(Container& x) noexcept : container(&x) {}
    
    back_emplace_iterator(const back_emplace_iterator&)            = default;
    back_emplace_iterator& operator=(const back_emplace_iterator&) = default;
    back_emplace_iterator(back_emplace_iterator&&)                 = default;
    back_emplace_iterator& operator=(back_emplace_iterator&&)      = default;
    
    ~back_emplace_iterator() = default;
    
    template <typename T>
    back_emplace_iterator<Container>& operator=(T&& t)
    {
        container->emplace_back(std::forward<T>(t));
        return *this;
    }
    
    back_emplace_iterator& operator*()     { return *this; }
    back_emplace_iterator& operator++()    { return *this; }
    back_emplace_iterator& operator++(int) { return *this; }
};

template <typename Container>
inline back_emplace_iterator<Container> back_emplacer(Container& c) noexcept
{
    return back_emplace_iterator<Container>(c);
}

} // namespace octopus

#endif
