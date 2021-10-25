#ifndef M0SH_REGULAR_H
#define M0SH_REGULAR_H
#pragma once

// std includes
#include <cmath>
// module includes
#include "m0sh/structured.h"

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class UnRegular : public Structured<TypeVector, TypeRef, TypeContainer> {
    public:
        UnRegular(const TypeContainer<TypeContainer<double>>& grid);
        UnRegular(const std::size_t dim, const TypeContainer<double>& axis);
        UnRegular(const TypeContainer<double>& axis);
    public:
        using Structured<TypeVector, TypeRef, TypeContainer>::x;
        TypeVector x(const TypeContainer<int>& ijk) const override;
        using Structured<TypeVector, TypeRef, TypeContainer>::ijk;
        TypeContainer<int> ijk(const TypeRef<const TypeVector>& x) const override;
    public:
        TypeContainer<TypeContainer<double>> grid;
        // inherited
        using Structured<TypeVector, TypeRef, TypeContainer>::n;
        using Structured<TypeVector, TypeRef, TypeContainer>::length;
        using Structured<TypeVector, TypeRef, TypeContainer>::origin;
};

// UnRegular

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
UnRegular<TypeVector, TypeRef, TypeContainer>::UnRegular(const TypeContainer<TypeContainer<double>>& p_grid) : Structured<TypeVector, TypeRef, TypeContainer>::Structured(TypeContainer<std::size_t>(p_grid.size()), TypeContainer<double>(p_grid.size()), TypeVector::Zero()), grid(p_grid) {
    for(std::size_t i = 0; i < grid.size(); i++) {
        n[i] = grid[i].size() - 1;
        length[i] = *(grid[i].end() - 1) - *(grid[i].begin());
        origin[i] = grid[i].front();
    }
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
UnRegular<TypeVector, TypeRef, TypeContainer>::UnRegular(const std::size_t dim, const TypeContainer<double>& axis) : UnRegular(TypeContainer<TypeContainer<double>>(dim, axis)) {
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
UnRegular<TypeVector, TypeRef, TypeContainer>::UnRegular(const TypeContainer<double>& axis) : UnRegular(1, axis) {
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeVector UnRegular<TypeVector, TypeRef, TypeContainer>::x(const TypeContainer<int>& ijk) const {
    // Compute ijkPeriodic
    TypeContainer<int> ijkPeriodic(ijk.size());
    TypeContainer<int> ijkPeriodicP(ijk.size());
    for(std::size_t i = 0; i < n.size(); i++) {
        ijkPeriodic[i] = ijk[i] % n[i];
        if(ijkPeriodic[i] < 0) {
            ijkPeriodic[i] += n[i];
        }
        ijkPeriodicP[i] = ijkPeriodic[i] + 1;
    }
    // Compute index
    TypeVector x;
    for(std::size_t i = 0; i < n.size(); i++) {
        x[i] = 0.5 * (grid[i][ijkPeriodic[i]] + grid[i][ijkPeriodicP[i]]);
    }
    return x;
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<int> UnRegular<TypeVector, TypeRef, TypeContainer>::ijk(const TypeRef<const TypeVector>& x) const {
    // Compute xPeriodic
    TypeVector xPeriodic = x;
    for(std::size_t i = 0; i < n.size(); i++) {
        xPeriodic[i] = std::fmod(x[i] - origin[i], length[i]);
        if(xPeriodic[i] < 0.0) {
            xPeriodic[i] += length[i];
        }
        xPeriodic[i] += origin[i];
    }
    // Compute index
    TypeContainer<int> ijk_(n.size());
    for(std::size_t i = 0; i < n.size(); i++) {
        ijk_[i] = std::distance(grid[i].begin(), std::lower_bound(grid[i].begin(), grid[i].end(), xPeriodic[i])) - 1;
    }
    return ijk_;
}

}

#endif
