#ifndef M0SH_MESH_H
#define M0SH_MESH_H
#pragma once

// std includes
#include <memory> // shared_ptr

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class Mesh {
    public:
        Mesh() {
        };
    public:
        // individual
        virtual TypeVector positionPoint(const std::size_t& index) const = 0;
        virtual TypeVector positionCell(const std::size_t& index) const = 0;
        virtual std::size_t indexPoint(const TypeRef<const TypeVector>& position) const = 0;
        virtual std::size_t indexCell(const TypeRef<const TypeVector>& position) const = 0;
        // nb
        virtual std::size_t nbPoints() const = 0;
        virtual std::size_t nbCells() const = 0;
        // all
        virtual TypeContainer<std::size_t> indexPoints() const {
            TypeContainer<std::size_t> indexs(nbPoints());
            std::iota(indexs.begin(), indexs.end(), 0);
            return indexs;
        };
        virtual TypeContainer<std::size_t> indexCells() const {
            TypeContainer<std::size_t> indexs(nbCells());
            std::iota(indexs.begin(), indexs.end(), 0);
            return indexs;
        };
        virtual TypeContainer<TypeVector> positionPoints() const {
            TypeContainer<TypeVector> positions(nbPoints());
            for(unsigned int i = 0; i < nbPoints(); i++) {
                positions[i] = positionPoint(i);
            }
            return positions;
        };
        virtual TypeContainer<TypeVector> positionCells() const {
            TypeContainer<TypeVector> positions(nbCells());
            for(unsigned int i = 0; i < nbCells(); i++) {
                positions[i] = positionCell(i);
            }
            return positions;
        };
};

}

#endif
