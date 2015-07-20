//
//  context_back_insert_iterator.h
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_context_back_insert_iterator_h
#define Octopus_context_back_insert_iterator_h

#include <iterator>
#include <memory> // std::addressof

template <typename Container>
class ContextBackInsertIterator : public std::iterator<std::output_iterator_tag, void, void, void, void>
{
protected:
    Container* container;
    
public:
    using container_type = Container;
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
    
    typename container_type::const_iterator begin() const { return container->begin(); }
    typename container_type::const_iterator end() const { return container->end(); }
    typename container_type::const_iterator cbegin() const { return container->cbegin(); }
    typename container_type::const_iterator cend() const { return container->cend(); }
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
    
    typename container_type::const_iterator begin() const { return container->begin(); }
    typename container_type::const_iterator end() const { return container->end(); }
    typename container_type::const_iterator cbegin() const { return container->cbegin(); }
    typename container_type::const_iterator cend() const { return container->cend(); }
};

template <class Container>
ContextInsertIterator<Container> ContextInserter(Container& x)
{
    return ContextInsertIterator<Container>(x);
}

#endif
