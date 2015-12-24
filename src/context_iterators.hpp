//
//  context_iterators.hpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_context_iterators_hpp
#define Octopus_context_iterators_hpp

#include <iterator>
#include <memory>

/*
 ContextIterator's are iterator adapters that allow iteration over a container, but also allow
 access to the underlying containers iterators, to provide 'context'.
 */

namespace Octopus
{
    template <typename Container>
    class ContextBackInsertIterator : public std::iterator<std::output_iterator_tag, void, void, void, void>
    {
    protected:
        Container* container;
        
    public:
        using container_type = Container;
        
        using iterator       = typename container_type::iterator;
        using const_iterator = typename container_type::const_iterator;
        
        explicit ContextBackInsertIterator (Container& x) : container(std::addressof(x)) {}
        ContextBackInsertIterator<Container>& operator= (const typename Container::value_type& value)
        { container->push_back(value); return *this; }
        ContextBackInsertIterator<Container>& operator= (typename Container::value_type&& value)
        { container->push_back(std::move(value)); return *this; }
        ContextBackInsertIterator<Container>& operator* ()
        { return *this; }
        ContextBackInsertIterator<Container>& operator++ ()
        { return *this; }
        ContextBackInsertIterator<Container> operator++ (int)
        { return *this; }
        
        iterator begin() { return container->begin(); }
        const_iterator begin() const { return container->begin(); }
        iterator end() { return container->end(); }
        const_iterator end() const { return container->end(); }
        const_iterator cbegin() { return container->cbegin(); }
        const_iterator cbegin() const { return container->cbegin(); }
        const_iterator cend() { return container->cend(); }
        const_iterator cend() const { return container->cend(); }
    };
    
    template <class Container>
    ContextBackInsertIterator<Container> ContextBackInserter(Container& x)
    {
        return ContextBackInsertIterator<Container>(x);
    }
    
    template <typename Container>
    class ContextInsertIterator : public std::iterator<std::output_iterator_tag, void, void, void, void>
    {
    protected:
        Container* container;
        
    public:
        using container_type = Container;
        
        using iterator       = typename container_type::iterator;
        using const_iterator = typename container_type::const_iterator;
        
        explicit ContextInsertIterator (Container& x) : container(std::addressof(x)) {}
        ContextInsertIterator<Container>& operator= (const typename Container::value_type& value)
        { container->insert(value); return *this; }
        ContextInsertIterator<Container>& operator= (typename Container::value_type&& value)
        { container->insert(std::move(value)); return *this; }
        ContextInsertIterator<Container>& operator* ()
        { return *this; }
        ContextInsertIterator<Container>& operator++ ()
        { return *this; }
        ContextInsertIterator<Container> operator++ (int)
        { return *this; }
        
        iterator begin() { return container->begin(); }
        const_iterator begin() const { return container->begin(); }
        iterator end() { return container->end(); }
        const_iterator end() const { return container->end(); }
        const_iterator cbegin() { return container->cbegin(); }
        const_iterator cbegin() const { return container->cbegin(); }
        const_iterator cend() { return container->cend(); }
        const_iterator cend() const { return container->cend(); }
    };
    
    template <class Container>
    ContextInsertIterator<Container> ContextInserter(Container& x)
    {
        return ContextInsertIterator<Container>(x);
    }
} // namespace Octopus

#endif
