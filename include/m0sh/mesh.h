#ifndef M0SH_MESH_H
#define M0SH_MESH_H
#pragma once

// std includes
#include <memory> // shared_ptr

namespace m0sh {

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
class Mesh {
    public:
        Mesh();
    public:
        virtual TypeVector x(const std::size_t& index) const = 0;
        virtual std::size_t index(const TypeRef<const TypeVector>& x) const = 0;
        virtual std::size_t size() const = 0;
        virtual TypeContainer<std::size_t> indexs() const;
        virtual TypeContainer<TypeVector> xs() const;
};

// Mesh

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
Mesh<TypeVector, TypeRef, TypeContainer>::Mesh() {
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<std::size_t> Mesh<TypeVector, TypeRef, TypeContainer>::indexs() const {
    TypeContainer<std::size_t> ids(size());
    std::iota(ids.begin(), ids.end(), 0);
    return ids;
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<TypeVector> Mesh<TypeVector, TypeRef, TypeContainer>::xs() const {
    TypeContainer<std::size_t> ids = indexs();
    TypeContainer<TypeVector> positions(ids.size());
    for(unsigned int i = 0; i < ids.size(); i++) {
        positions[i] = x(i);
    }
    return positions;
}

}

#endif
