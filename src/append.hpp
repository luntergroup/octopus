#ifndef append_hpp
#define append_hpp

#include <vector>
#include <iterator>
#include <utility>

template <typename T>
typename std::vector<T>::iterator append(const std::vector<T>& src, std::vector<T>& dest)
{
    typename std::vector<T>::iterator result;

    if (dest.empty()) {
        dest = src;
        result = std::begin(dest);
    } else {
        result = dest.insert(std::end(dest), std::cbegin(src), std::cend(src));
    }

    return result;
}

template <typename T>
typename std::vector<T>::iterator append(std::vector<T>&& src, std::vector<T>& dest)
{
    typename std::vector<T>::iterator result;

    if (dest.empty()) {
        dest = std::move(src);
        result = std::begin(dest);
    } else {
        result = dest.insert(std::end(dest), std::make_move_iterator(std::begin(src)), std::make_move_iterator(std::end(src)));
    }

    src.clear();
    src.shrink_to_fit();

    return result;
}

#endif
