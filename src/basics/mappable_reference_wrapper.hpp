//
//  mappable_reference_wrapper.hpp
//  octopus
//
//  Created by Daniel Cooke on 02/08/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef mappable_reference_wrapper_hpp
#define mappable_reference_wrapper_hpp

#include <functional>
#include <type_traits>

#include <concepts/mappable.hpp>

namespace octopus {

/**
 MappableReferenceWrapper is a wrapper for Mappable types that is also Mappable. They
 can therefore be used in containers and mappable algorithms.
 */
template <typename T>
class MappableReferenceWrapper : public Mappable<MappableReferenceWrapper<T>>
{
public:
    static_assert(is_mappable<std::remove_cv_t<T>>, "not a Mappable");
    
    using type = T;
    
    MappableReferenceWrapper() = delete;
    
    MappableReferenceWrapper(T& mappable) : mappable_ {mappable} {}
    
    MappableReferenceWrapper(T&&) = delete;
    
    MappableReferenceWrapper(const MappableReferenceWrapper&)            = default;
    MappableReferenceWrapper& operator=(const MappableReferenceWrapper&) = default;
    MappableReferenceWrapper(MappableReferenceWrapper&&)                 = default;
    MappableReferenceWrapper& operator=(MappableReferenceWrapper&&)      = default;
    
    operator T&() const noexcept { return mappable_.get(); }
    
    T& get() const noexcept { return mappable_.get(); }
    
    decltype(auto) mapped_region() const noexcept { return mappable_.get().mapped_region(); }
    
private:
    std::reference_wrapper<T> mappable_;
};

} // namespace octopus

namespace std {
    template <typename T>
    struct hash<octopus::MappableReferenceWrapper<T>>
    {
        size_t operator()(const octopus::MappableReferenceWrapper<T> ref) const
        {
            return hash<T>()(ref);
        }
    };
} // namespace std

#endif
