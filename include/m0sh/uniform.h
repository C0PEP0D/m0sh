#ifndef M0SH_UNIFORM_H
#define M0SH_UNIFORM_H
#pragma once

// std includes
#include <cmath>
// module includes
#include "m0sh/structured.h"

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class Uniform : public Structured<TypeVector, TypeRef, TypeContainer> {
    public:
        using TypeInherited = Structured<TypeVector, TypeRef, TypeContainer>;
    public:
        Uniform(const TypeContainer<std::size_t>& p_nCells, const TypeContainer<double>& p_length, const TypeVector& p_origin, const TypeContainer<bool>& p_periodic) {
            // copy
            nCells = p_nCells;
            length = p_length;
            origin = p_origin;
            periodic = p_periodic;
            // nPoints
            nPoints.resize(nCells.size());
            for(unsigned int i = 0; i < nCells.size(); i++) {
                if(periodic[i]) {
                    nPoints[i] = nCells[i];
                } else {
                    nPoints[i] = nCells[i] + 1;
                } 
            }
            // spacing between points
            spacing.resize(nCells.size());
            for(std::size_t i = 0; i < nCells.size(); i++) {
                spacing[i] = length[i] / nCells[i];
            }
        };
    public:
        using TypeInherited::positionPoint;
        using TypeInherited::positionCell;
        TypeVector positionPoint(const TypeContainer<int>& ijk) const override {
            TypeVector position;
            for(std::size_t i = 0; i < spacing.size(); i++) {
                position[i] = ijk[i] * spacing[i];
                //if(periodic[i]) {
                //    if(position[i] < 0.0) {
                //        position[i] = length[i] - std::fmod(std::abs(position[i]), length[i]);
                //    } else {
                //        position[i] = std::fmod(std::abs(position[i]), length[i]);
                //    }
                //}
                position[i] += origin[i];
            }
            return position;
        };
        TypeVector positionCell(const TypeContainer<int>& ijk) const override {
            TypeVector position;
            for(std::size_t i = 0; i < spacing.size(); i++) {
                position[i] = ijk[i] * spacing[i] + 0.5 * spacing[i];
                //if(periodic[i]) {
                //    if(position[i] < 0.0) {
                //        position[i] = length[i] - std::fmod(std::abs(position[i]), length[i]);
                //    } else {
                //        position[i] = std::fmod(std::abs(position[i]), length[i]);
                //    }
                //}
                position[i] += origin[i];
            }
            return position;
        };
        using TypeInherited::ijkPoint;
        using TypeInherited::ijkCell;
        TypeContainer<int> ijkPoint(const TypeRef<const TypeVector>& position) const override {
            TypeContainer<int> ijk_(nPoints.size());
            for(std::size_t i = 0; i < nPoints.size(); i++) {
                ijk_[i] = std::floor((position[i] - origin[i] - 0.5 * spacing[i]) / spacing[i]);
            }
            return ijk_;
        };
        TypeContainer<int> ijkCell(const TypeRef<const TypeVector>& position) const override {
            TypeContainer<int> ijk_(nCells.size());
            for(std::size_t i = 0; i < nCells.size(); i++) {
                ijk_[i] = std::floor((position[i] - origin[i]) / spacing[i]);
            }
            return ijk_;
        };
    public:
        TypeContainer<double> spacing;
        // inherited
        using Structured<TypeVector, TypeRef, TypeContainer>::nPoints;
        using Structured<TypeVector, TypeRef, TypeContainer>::nCells;
        using Structured<TypeVector, TypeRef, TypeContainer>::length;
        using Structured<TypeVector, TypeRef, TypeContainer>::origin;
        using Structured<TypeVector, TypeRef, TypeContainer>::periodic;
};

}

#endif
