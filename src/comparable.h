//
//  comparable.h
//  Octopus
//
//  Created by Daniel Cooke on 13/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_comparable_h
#define Octopus_comparable_h

/**
    A little bit of voodoo template magic. See the curiously recurring template pattern!
    The crux of this is if the class T provides the '==' and '<' operators, then it is
    fully comparable, and by deriving from this class, will automatically have the
    '!=', '>', '<=', & '>=' provided.
 */

template <typename T>
class Comparable {};

template <typename T>
inline bool operator!=(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return !operator==(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

template <typename T>
inline bool operator>(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return operator<(static_cast<const T&>(rhs), static_cast<const T&>(lhs));
}

template <typename T>
inline bool operator<=(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return !operator>(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

template <typename T>
inline bool operator>=(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return !operator<(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

#endif
