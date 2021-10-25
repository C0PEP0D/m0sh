#ifndef M0SH_REGULAR_H
#define M0SH_REGULAR_H
#pragma once

// std includes
#include <cmath>
// module includes
#include "m0sh/structured.h"

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class Regular : public Structured<TypeVector, TypeRef, TypeContainer> {
    public:
        Regular(const TypeContainer<std::size_t>& n, const TypeContainer<double>& length, const TypeVector& origin);
    public:
        using Structured<TypeVector, TypeRef, TypeContainer>::x;
        TypeVector x(const TypeContainer<int>& ijk) const override;
        using Structured<TypeVector, TypeRef, TypeContainer>::ijk;
        TypeContainer<int> ijk(const TypeRef<const TypeVector>& x) const override;
    public:
        TypeContainer<double> spacing;
        // inherited
        using Structured<TypeVector, TypeRef, TypeContainer>::n;
        using Structured<TypeVector, TypeRef, TypeContainer>::length;
        using Structured<TypeVector, TypeRef, TypeContainer>::origin;
};

// Regular

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
Regular<TypeVector, TypeRef, TypeContainer>::Regular(const TypeContainer<std::size_t>& p_n, const TypeContainer<double>& p_length, const TypeVector& p_origin) : Structured<TypeVector, TypeRef, TypeContainer>::Structured(p_n, p_length, p_origin), spacing(n.size()) {
    for(std::size_t i = 0; i < n.size(); i++) {
        spacing[i] = length[i] / n[i];
    }
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeVector Regular<TypeVector, TypeRef, TypeContainer>::x(const TypeContainer<int>& ijk) const {
    // Compute cell center
    TypeVector x;
    for(std::size_t i = 0; i < spacing.size(); i++) {
        x[i] = origin[i] + ijk[i] * spacing[i] + 0.5 * spacing[i];
    }
    return x;
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<int> Regular<TypeVector, TypeRef, TypeContainer>::ijk(const TypeRef<const TypeVector>& x) const {
    // Compute index
    TypeContainer<int> ijk_(n.size());
    for(std::size_t i = 0; i < n.size(); i++) {
        ijk_[i] = std::floor((x[i] - origin[i]) / spacing[i]);
    }
    return ijk_;
}

}

#endif
