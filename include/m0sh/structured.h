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
        Structured() {
        };
    public:
        std::size_t nbPoints() const override {
            return std::accumulate(nPoints.begin(), nPoints.end(), 1.0, std::multiplies<double>());
        };
        std::size_t nbCells() const override {
            return std::accumulate(nCells.begin(), nCells.end(), 1.0, std::multiplies<double>());
        };
        
        std::size_t indexPoint(const TypeRef<const TypeVector>& position) const override {
            return indexPoint(ijkPoint(position));
        };
        std::size_t indexCell(const TypeRef<const TypeVector>& position) const override {
            return indexCell(ijkCell(position));
        };
        virtual std::size_t indexPoint(const TypeContainer<int>& ijk) const {
            // Compute ijk periodic
            TypeContainer<int> ijkPeriodic = ijk;
            for(std::size_t i = 0; i < nPoints.size(); i++) {
                ijkPeriodic[i] = std::abs(ijk[i]) % nPoints[i];
                if(ijk[i] < 0 && ijkPeriodic[i] != 0) {
                    ijkPeriodic[i] = nPoints[i] - ijkPeriodic[i];
                }
            }
            // Compute index
            std::size_t nTot = 1;
            std::size_t index = ijkPeriodic[0];
            for(std::size_t i = 1; i < ijkPeriodic.size(); i++) {
                nTot *= nPoints[i-1];
                index += nTot * ijkPeriodic[i];
            }
            return index;
        };
        virtual std::size_t indexCell(const TypeContainer<int>& ijk) const {
            // Compute ijk periodic
            TypeContainer<int> ijkPeriodic = ijk;
            for(std::size_t i = 0; i < nCells.size(); i++) {
                ijkPeriodic[i] = std::abs(ijk[i]) % nCells[i];
                if(ijk[i] < 0 && ijkPeriodic[i] != 0) {
                    ijkPeriodic[i] = nCells[i] - ijkPeriodic[i];
                }
            }
            // Compute index
            std::size_t nTot = 1;
            std::size_t index = ijkPeriodic[0];
            for(std::size_t i = 1; i < ijkPeriodic.size(); i++) {
                nTot *= nCells[i-1];
                index += nTot * ijkPeriodic[i];
            }
            return index;
        };
        
        virtual TypeContainer<int> ijkPoint(const TypeRef<const TypeVector>& position) const = 0;
        virtual TypeContainer<int> ijkCell(const TypeRef<const TypeVector>& position) const = 0;
        virtual TypeContainer<int> ijkPoint(const std::size_t& index) const {
            // Compute nDiv
            TypeContainer<std::size_t> nDiv(nPoints.size(), 1);
            for(std::size_t k = 1; k < nPoints.size(); k++) {
                nDiv[k] = nDiv[k-1] * nPoints[k-1];
            }
            // Compute ijk
            std::ldiv_t dv = {0, (long int)index};
            TypeContainer<int> ijk_(nPoints.size());
            for(std::size_t k = 0; k < nPoints.size(); k++) {
                dv = std::div(dv.rem, nDiv[nPoints.size() - k - 1]);
                ijk_[nPoints.size() - k - 1] = dv.quot;
            }
            return ijk_;
        };
        virtual TypeContainer<int> ijkCell(const std::size_t& index) const {
            // Compute nDiv
            TypeContainer<std::size_t> nDiv(nCells.size(), 1);
            for(std::size_t k = 1; k < nCells.size(); k++) {
                nDiv[k] = nDiv[k-1] * nCells[k-1];
            }
            // Compute ijk
            std::ldiv_t dv = {0, (long int)index};
            TypeContainer<int> ijk_(nCells.size());
            for(std::size_t k = 0; k < nCells.size(); k++) {
                dv = std::div(dv.rem, nDiv[nCells.size() - k - 1]);
                ijk_[nCells.size() - k - 1] = dv.quot;
            }
            return ijk_;
        };
        
        virtual TypeVector positionPoint(const TypeContainer<int>& ijk) const = 0;
        virtual TypeVector positionCell(const TypeContainer<int>& ijk) const = 0;
        TypeVector positionPoint(const std::size_t& index) const override {
            return positionPoint(ijkPoint(index));
        };
        TypeVector positionCell(const std::size_t& index) const override {
            return positionCell(ijkCell(index));
        };
    public:
        TypeContainer<std::size_t> nPoints;
        TypeContainer<std::size_t> nCells; // if periodic, nCells[i] = nPoints[i] else nCells[i] =  nPoints[i] - 1
        TypeContainer<double> length;
        TypeVector origin;
        TypeContainer<bool> periodic;
};

}

#endif
