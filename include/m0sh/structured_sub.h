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
        StructuredSub(const TypeContainer<std::size_t>& n, const TypeContainer<int>& offset, const std::shared_ptr<Structured<TypeVector, TypeRef, TypeContainer>>& sMesh);
    public:
        using Structured<TypeVector, TypeRef, TypeContainer>::index;
        std::size_t index(const TypeContainer<int>& ijk) const override;
        using Structured<TypeVector, TypeRef, TypeContainer>::ijk;
        TypeContainer<int> ijk(const std::size_t& index) const override;
        TypeContainer<int> ijk(const TypeRef<const TypeVector>& x) const override;
        using Structured<TypeVector, TypeRef, TypeContainer>::x;
        TypeVector x(const TypeContainer<int>& ijk) const override;
        
        TypeContainer<std::size_t> indexs() const override;
    public:
        TypeContainer<int> offset;
        TypeContainer<int> offsetPeriodic;
        std::shared_ptr<Structured<TypeVector, TypeRef, TypeContainer>> sMesh;
        using Structured<TypeVector, TypeRef, TypeContainer>::n;
        using Structured<TypeVector, TypeRef, TypeContainer>::length;
        using Structured<TypeVector, TypeRef, TypeContainer>::origin;
};

// StructuredSub

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
StructuredSub<TypeVector, TypeRef, TypeContainer>::StructuredSub(const TypeContainer<std::size_t>& p_n, const TypeContainer<int>& p_offset, const std::shared_ptr<Structured<TypeVector, TypeRef, TypeContainer>>& p_sMesh) : Structured<TypeVector, TypeRef, TypeContainer>(p_n, TypeContainer<double>(p_n.size()), TypeVector()), offset(p_offset), sMesh(p_sMesh) {
    // Compute length and origin
    TypeContainer<int> sup(n.size(), 0);
    TypeContainer<int> inf(n.size(), 0);
    TypeContainer<int> supInf(n.size(), 0);
    TypeContainer<int> infSup(n.size(), 0);
    for(std::size_t i = 0; i < n.size(); i++) {
        sup[i] = offset[i] + n[i];
        supInf[i] = offset[i] + n[i] - 1;
        inf[i] = offset[i];
        infSup[i] = offset[i] + 1;
        length[i] = (sMesh->x(sup) - sMesh->x(inf))[i] + 0.5 * (sMesh->x(sup) - sMesh->x(supInf))[i] + 0.5 * (sMesh->x(infSup) - sMesh->x(inf))[i]; // not perfect
    }
    origin =  sMesh->x(offset);
    // compute offset periodic
    offsetPeriodic = TypeContainer<int>(offset.size());
    for(std::size_t i = 0; i < offset.size(); i++) {
        offsetPeriodic[i] = std::abs(offset[i]) % sMesh->n[i];
        if(offset[i] < 0) {
            offsetPeriodic[i] = sMesh->n[i] - offsetPeriodic[i];
        }
    }
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
std::size_t StructuredSub<TypeVector, TypeRef, TypeContainer>::index(const TypeContainer<int>& ijk) const {
    TypeContainer<int> ijk_ = ijk;
    std::transform(ijk_.begin(), ijk_.end(), offset.begin(), ijk_.begin(), std::plus<int>());
    return sMesh->index(ijk_);
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<int> StructuredSub<TypeVector, TypeRef, TypeContainer>::ijk(const std::size_t& index) const {
    // compute ijk
    TypeContainer<int> ijk_ = sMesh->ijk(index);
    for(std::size_t i = 0; i < ijk_.size(); i++) {
        ijk_[i] -= offsetPeriodic[i];
        if(ijk_[i] < 0) {
            ijk_[i] = ijk_[i] + sMesh->n[i];
        } else if (ijk_[i] >= n[i]) {
            ijk_[i] = ijk_[i] - sMesh->n[i];
        }
    }
    return ijk_;
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<int> StructuredSub<TypeVector, TypeRef, TypeContainer>::ijk(const TypeRef<const TypeVector>& x) const {
    TypeContainer<int> ijk_ = sMesh->ijk(x);
    std::transform(ijk_.begin(), ijk_.end(), offset.begin(), ijk_.begin(), std::minus<int>());
    return ijk_;
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeVector StructuredSub<TypeVector, TypeRef, TypeContainer>::x(const TypeContainer<int>& ijk) const {
    TypeContainer<int> ijk_ = ijk;
    std::transform(ijk_.begin(), ijk_.end(), offset.begin(), ijk_.begin(), std::plus<int>());
    return sMesh->x(ijk_);
}

template<typename TypeVector, template<typename...> class TypeRef, template<typename...> class TypeContainer>
TypeContainer<std::size_t> StructuredSub<TypeVector, TypeRef, TypeContainer>::indexs() const {
    TypeContainer<std::size_t> ids = Structured<TypeVector, TypeRef, TypeContainer>::indexs();
    for(auto& i : ids) {
        i = index(Structured<TypeVector, TypeRef, TypeContainer>::ijk(i));
    }
    return ids;
}

}

#endif
