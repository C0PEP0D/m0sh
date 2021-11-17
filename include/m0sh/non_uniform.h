#ifndef M0SH_NON_UNIFORM_H
#define M0SH_NON_UNIFORM_H
#pragma once

// std includes
#include <cmath>
// module includes
#include "m0sh/structured.h"

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class NonUniform : public Structured<TypeVector, TypeRef, TypeContainer> {
    public:
        using TypeInherited = Structured<TypeVector, TypeRef, TypeContainer>;
    public:
        NonUniform(const TypeContainer<TypeContainer<double>>& p_gridPoints, const TypeContainer<bool>& p_periodic) : gridPoints(p_gridPoints) {
            // copy
            periodic = p_periodic;
            // grid
            gridCells.resize(gridPoints.size());
            for(unsigned int i = 0; i < gridCells.size(); i++) {
                gridCells[i].resize(gridPoints[i].size() - 1);
                for(unsigned int j = 0; j < gridCells[i].size(); j++) {
                    gridCells[i][j] = 0.5 * (gridPoints[i][j] + gridPoints[i][j + 1]);
                }
            }
            // origin
            for(unsigned int i = 0; i < gridPoints.size(); i++) {
                origin[i] = gridPoints[i][0];
            }
            // length
            length.resize(gridPoints.size());
            for(unsigned int i = 0; i < gridPoints.size(); i++) {
                length[i] = gridPoints[i][gridPoints[i].size() - 1] - gridPoints[i][0];
            }
            // points and cells
            nPoints.resize(gridPoints.size());
            nCells.resize(gridPoints.size());
            for(unsigned int i = 0; i < nPoints.size(); i++) {
                if(periodic[i]) {
                    nPoints[i] = gridPoints[i].size() - 1;
                    nCells[i] = gridPoints[i].size() - 1;
                } else {
                    nPoints[i] = gridPoints[i].size();
                    nCells[i] = gridPoints[i].size() - 1;
                }
            }
        };

        NonUniform(const std::size_t& dim, const TypeContainer<double>& axisPoints, bool p_periodic) : NonUniform(TypeContainer<TypeContainer<double>>(dim, axisPoints), TypeContainer<bool>(dim, p_periodic)) {
        };

        NonUniform(const TypeContainer<double>& axisPoints, bool p_periodic) : NonUniform(1, axisPoints, p_periodic) {
        };
    public:
        using TypeInherited::positionPoint;
        using TypeInherited::positionCell;
        TypeVector positionPoint(const TypeContainer<int>& ijk) const override {
            // Compute ijkPeriodic
            TypeContainer<int> ijkPeriodic(ijk.size());
            TypeContainer<int> ijkDiv(ijk.size());
            for(std::size_t i = 0; i < nPoints.size(); i++) {
                ijkPeriodic[i] = std::abs(ijk[i]) % nPoints[i];
                ijkDiv[i] = std::abs(ijk[i]) / nPoints[i];
                if(ijk[i] < 0) {
                    ijkPeriodic[i] = nPoints[i] - ijkPeriodic[i];
                    ijkDiv[i] = -ijkDiv[i] - 1;
                }
            }
            // Compute position
            TypeVector position;
            for(std::size_t i = 0; i < nPoints.size(); i++) {
                position[i] = gridPoints[i][ijkPeriodic[i]] + ijkDiv[i] * length[i];
            }
            return position;
        };
        TypeVector positionCell(const TypeContainer<int>& ijk) const override {
            // Compute ijkPeriodic
            TypeContainer<int> ijkPeriodic(ijk.size());
            TypeContainer<int> ijkDiv(ijk.size());
            for(std::size_t i = 0; i < nCells.size(); i++) {
                ijkPeriodic[i] = std::abs(ijk[i]) % nCells[i];
                ijkDiv[i] = std::abs(ijk[i]) / nCells[i];
                if(ijk[i] < 0) {
                    ijkPeriodic[i] = nCells[i] - ijkPeriodic[i];
                    ijkDiv[i] = -ijkDiv[i] - 1;
                }
            }
            // Compute position
            TypeVector position;
            for(std::size_t i = 0; i < nCells.size(); i++) {
                position[i] = gridCells[i][ijkPeriodic[i]] + ijkDiv[i] * length[i];
            }
            return position;
        };
        using TypeInherited::ijkPoint;
        using TypeInherited::ijkCell;
        TypeContainer<int> ijkPoint(const TypeRef<const TypeVector>& position) const override {
            // Compute positionPeriodic
            TypeVector positionPeriodic;
            TypeContainer<int> positionDiv(positionPeriodic.size());
            for(std::size_t i = 0; i < nPoints.size(); i++) {
                positionPeriodic[i] = std::fmod(std::abs(position[i] - origin[i]), length[i]);
                positionDiv[i] = std::floor(std::abs(position[i] - origin[i]) / length[i]);
                if(position[i] < origin[i]) {
                    positionPeriodic[i] = length[i] - positionPeriodic[i];
                    positionDiv[i] = -positionDiv[i] - 1;
                }
                positionPeriodic[i] += origin[i];
            }
            // Compute index
            TypeContainer<int> ijk_(nPoints.size());
            for(std::size_t i = 0; i < nPoints.size(); i++) {
                ijk_[i] = std::distance(gridCells[i].begin(), std::lower_bound(gridCells[i].begin(), gridCells[i].end(), positionPeriodic[i])) + positionDiv[i] * nPoints[i];
            }
            return ijk_;
        };
        TypeContainer<int> ijkCell(const TypeRef<const TypeVector>& position) const override {
            // Compute positionPeriodic
            TypeVector positionPeriodic;
            TypeContainer<int> positionDiv(positionPeriodic.size());
            for(std::size_t i = 0; i < nCells.size(); i++) {
                positionPeriodic[i] = std::fmod(std::abs(position[i] - origin[i]), length[i]);
                positionDiv[i] = std::floor(std::abs(position[i] - origin[i]) / length[i]);
                if(position[i] < origin[i]) {
                    positionPeriodic[i] = length[i] - positionPeriodic[i];
                    positionDiv[i] = -positionDiv[i] - 1;
                }
                positionPeriodic[i] += origin[i];
            }
            // Compute index
            TypeContainer<int> ijk_(nCells.size());
            for(std::size_t i = 0; i < nCells.size(); i++) {
                ijk_[i] = std::distance(gridPoints[i].begin(), std::lower_bound(gridPoints[i].begin(), gridPoints[i].end(), positionPeriodic[i])) - 1 + positionDiv[i] * nCells[i];
            }
            return ijk_;
        };
    public:
        TypeContainer<TypeContainer<double>> gridPoints;
        TypeContainer<TypeContainer<double>> gridCells;
        // inherited
        using TypeInherited::nPoints;
        using TypeInherited::nCells;
        using TypeInherited::length;
        using TypeInherited::origin;
        using TypeInherited::periodic;
};

}

#endif
