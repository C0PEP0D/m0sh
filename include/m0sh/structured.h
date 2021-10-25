#ifndef M0SH_STRUCTURED_H
#define M0SH_STRUCTURED_H
#pragma once

// std includes
#include <cmath> // div, ldiv_t
#include <numeric> // accumulate
#include <functional> // multiplies
// module include
#include "m0sh/mesh.h"

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class Structured : public Mesh<TypeVector, TypeRef, TypeContainer> {
    public:
        Structured(const TypeContainer<std::size_t>& n, const TypeContainer<double>& length, const TypeVector& origin);
    public:
        std::size_t size() const override;
        
        std::size_t index(const TypeRef<const TypeVector>& x) const override;
        virtual std::size_t index(const TypeContainer<int>& ijk) const;
        
        virtual TypeContainer<int> ijk(const TypeRef<const TypeVector>& x) const = 0;
        virtual TypeContainer<int> ijk(const std::size_t& index) const;
        
        virtual TypeVector x(const TypeContainer<int>& ijk) const = 0;
        TypeVector x(const std::size_t& index) const override;
    public:
        TypeContainer<std::size_t> n;
        TypeContainer<double> length;
        TypeVector origin;
};

// Structured

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
Structured<TypeVector, TypeRef, TypeContainer>::Structured(const TypeContainer<std::size_t>& p_n, const TypeContainer<double>& p_length, const TypeVector& p_origin) : n(p_n), length(p_length), origin(p_origin) {
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
std::size_t Structured<TypeVector, TypeRef, TypeContainer>::size() const {
    return std::accumulate(n.begin(), n.end(), 1.0, std::multiplies<double>());
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
std::size_t Structured<TypeVector, TypeRef, TypeContainer>::index(const TypeRef<const TypeVector>& x) const {
    return index(ijk(x));
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
std::size_t Structured<TypeVector, TypeRef, TypeContainer>::index(const TypeContainer<int>& ijk) const {
    // Compute ijk periodic
    TypeContainer<int> ijkPeriodic = ijk;
    for(std::size_t i = 0; i < n.size(); i++) {
        ijkPeriodic[i] = std::abs(ijk[i]) % n[i];
        if(ijk[i] < 0 && ijkPeriodic[i] != 0) {
            ijkPeriodic[i] = n[i] - ijkPeriodic[i];
        }
    }
    // Compute index
    std::size_t nTot = 1;
    std::size_t index = ijkPeriodic[0];
    for(std::size_t i = 1; i < ijkPeriodic.size(); i++) {
        nTot *= n[i-1];
        index += nTot * ijkPeriodic[i];
    }
    return index;
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<int> Structured<TypeVector, TypeRef, TypeContainer>::ijk(const std::size_t& index) const {
    // Compute nDiv
    TypeContainer<std::size_t> nDiv(n.size(), 1);
    for(std::size_t k = 1; k < n.size(); k++) {
        nDiv[k] = nDiv[k-1] * n[k-1];
    }
    // Compute ijk
    std::ldiv_t dv = {0, (long int)index};
    TypeContainer<int> ijk_(n.size());
    for(std::size_t k = 0; k < n.size(); k++) {
        dv = std::div(dv.rem, nDiv[n.size() - k - 1]);
        ijk_[n.size() - k - 1] = dv.quot;
    }
    return ijk_;
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeVector Structured<TypeVector, TypeRef, TypeContainer>::x(const std::size_t& index) const {
    return x(ijk(index));
}

}

#endif
