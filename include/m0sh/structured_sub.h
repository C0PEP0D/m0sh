#ifndef M0SH_STRUCTURED_SUB_H
#define M0SH_STRUCTURED_SUB_H
#pragma once

// std includes
#include <memory> // shared_ptr
// module include
#include "m0sh/structured.h"

// TODO assert

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class StructuredSub : public Structured<TypeVector, TypeRef, TypeContainer> {
    public:
        using TypeInherited = Structured<TypeVector, TypeRef, TypeContainer>;
    public:
        StructuredSub(const TypeContainer<std::size_t>& p_nCells, const TypeContainer<int>& p_offsetCells, const std::shared_ptr<TypeInherited>& p_sMesh) : offsetCells(p_offsetCells), sMesh(p_sMesh) {
            // copy
            nCells = p_nCells;
            periodic = TypeContainer<bool>(nCells.size(), true);
            // points
            nPoints.resize(nCells.size());
            for(unsigned int i = 0; i < nCells.size(); i++) {
                if(periodic[i]) {
                    nPoints[i] = nCells[i];
                } else {
                    nPoints[i] = nCells[i] + 1;
                } 
            }
            offsetPoints = offsetCells;
            // length and origin
            std::vector<int> offsetEndPoints(offsetPoints.size());
            for(std::size_t i = 0; i < nCells.size(); i++) {
                if(periodic[i]) {
                    offsetEndPoints[i] = offsetPoints[i] + nPoints[i];
                } else {
                    offsetEndPoints[i] = offsetPoints[i] + nPoints[i] - 1;
                }
            }
            length.resize(offsetCells.size());
            TypeVector _length = (sMesh->positionPoint(offsetEndPoints) - sMesh->positionPoint(offsetPoints));
            for(std::size_t i = 0; i < nCells.size(); i++) {
                length[i] = _length[i];
            }
            origin = sMesh->positionPoint(offsetPoints);
            // compute offset periodic
            offsetPointsPeriodic = TypeContainer<int>(offsetPoints.size());
            for(std::size_t i = 0; i < offsetPoints.size(); i++) {
                offsetPointsPeriodic[i] = std::abs(offsetPoints[i]) % sMesh->nPoints[i];
                if(offsetPoints[i] < 0) {
                    offsetPointsPeriodic[i] = sMesh->nPoints[i] - offsetPointsPeriodic[i];
                }
            }
            offsetCellsPeriodic = TypeContainer<int>(offsetCells.size());
            for(std::size_t i = 0; i < offsetCells.size(); i++) {
                offsetCellsPeriodic[i] = std::abs(offsetCells[i]) % sMesh->nCells[i];
                if(offsetCells[i] < 0) {
                    offsetCellsPeriodic[i] = sMesh->nCells[i] - offsetCellsPeriodic[i];
                }
            }
        }
    public:
        using TypeInherited::indexPoint;
        using TypeInherited::indexCell;
        std::size_t indexPoint(const TypeContainer<int>& ijk) const override {
            TypeContainer<int> ijk_ = ijk;
            std::transform(ijk_.begin(), ijk_.end(), offsetPoints.begin(), ijk_.begin(), std::plus<int>());
            return sMesh->indexPoint(ijk_);
        };
        std::size_t indexCell(const TypeContainer<int>& ijk) const override {
            TypeContainer<int> ijk_ = ijk;
            std::transform(ijk_.begin(), ijk_.end(), offsetCells.begin(), ijk_.begin(), std::plus<int>());
            return sMesh->indexCell(ijk_);
        };
        using TypeInherited::ijkPoint;
        using TypeInherited::ijkCell;
        TypeContainer<int> ijkPoint(const std::size_t& index) const override {
            // compute ijk
            TypeContainer<int> ijk_ = sMesh->ijkPoint(index);
            for(std::size_t i = 0; i < ijk_.size(); i++) {
                ijk_[i] -= offsetPointsPeriodic[i];
                if(ijk_[i] < 0) {
                    ijk_[i] = ijk_[i] + sMesh->nPoints[i];
                } else if (ijk_[i] >= nPoints[i]) {
                    ijk_[i] = ijk_[i] - sMesh->nPoints[i];
                }
            }
            return ijk_;
        };
        TypeContainer<int> ijkCell(const std::size_t& index) const override {
            // compute ijk
            TypeContainer<int> ijk_ = sMesh->ijkCell(index);
            for(std::size_t i = 0; i < ijk_.size(); i++) {
                ijk_[i] -= offsetCellsPeriodic[i];
                if(ijk_[i] < 0) {
                    ijk_[i] = ijk_[i] + sMesh->nCells[i];
                } else if (ijk_[i] >= nCells[i]) {
                    ijk_[i] = ijk_[i] - sMesh->nCells[i];
                }
            }
            return ijk_;
        };
        TypeContainer<int> ijkPoint(const TypeRef<const TypeVector>& position) const override {
            TypeContainer<int> ijk_ = sMesh->ijkPoint(position);
            std::transform(ijk_.begin(), ijk_.end(), offsetPoints.begin(), ijk_.begin(), std::minus<int>());
            return ijk_;
        };
        TypeContainer<int> ijkCell(const TypeRef<const TypeVector>& position) const override {
            TypeContainer<int> ijk_ = sMesh->ijkCell(position);
            std::transform(ijk_.begin(), ijk_.end(), offsetCells.begin(), ijk_.begin(), std::minus<int>());
            return ijk_;
        };
        using TypeInherited::positionPoint;
        using TypeInherited::positionCell;
        TypeVector positionPoint(const TypeContainer<int>& ijk) const override {
            TypeContainer<int> ijk_ = ijk;
            std::transform(ijk_.begin(), ijk_.end(), offsetPoints.begin(), ijk_.begin(), std::plus<int>());
            return sMesh->positionPoint(ijk_);
        };
        TypeVector positionCell(const TypeContainer<int>& ijk) const override {
            TypeContainer<int> ijk_ = ijk;
            std::transform(ijk_.begin(), ijk_.end(), offsetCells.begin(), ijk_.begin(), std::plus<int>());
            return sMesh->positionCell(ijk_);
        };
        
        TypeContainer<std::size_t> indexPoints() const override {
            TypeContainer<std::size_t> indexs = TypeInherited::indexPoints();
            for(auto& i : indexs) {
                i = indexPoint(TypeInherited::ijkPoint(i));
            }
            return indexs;
        };
        TypeContainer<std::size_t> indexCells() const override {
            TypeContainer<std::size_t> indexs = TypeInherited::indexCells();
            for(auto& i : indexs) {
                i = indexCell(TypeInherited::ijkCell(i));
            }
            return indexs;
        };
    public:
        TypeContainer<int> offsetPoints;
        TypeContainer<int> offsetCells;
        TypeContainer<int> offsetPointsPeriodic;
        TypeContainer<int> offsetCellsPeriodic;
        std::shared_ptr<TypeInherited> sMesh;
        using TypeInherited::nPoints;
        using TypeInherited::nCells;
        using TypeInherited::length;
        using TypeInherited::origin;
        using TypeInherited::periodic;
};

}

#endif
