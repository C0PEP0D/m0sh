// std includes
#include <iostream>
#include <vector>
#include <random>
// thirdparties includes
#include <Eigen/Dense>
// lib includes
#include "m0sh/uniform.h"
#include "m0sh/structured_sub.h"

using TypeScalar = double;
// Space
const unsigned int DIM = 3;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
// Mesh
template<typename ...Args>
using TypeContainer = std::vector<Args...>;
using TypeStructured = m0sh::Structured<TypeVector, TypeRef, TypeContainer>;
using TypeMesh = m0sh::Uniform<TypeVector, TypeRef, TypeContainer>;
using TypeSub = m0sh::StructuredSub<TypeVector, TypeRef, TypeContainer>;
// Data
// // Uniform mesh
const std::size_t n = 100;
const double length = 1.0;
const double origin = 0.5;
// // 
const std::size_t nSub = 10;
const TypeContainer<int> offSub = {0, 0, 0};

void print(const std::shared_ptr<TypeSub>& sMesh, std::uniform_int_distribution<int>& uniform, std::default_random_engine& e) {
    TypeContainer<int> ijk;
    TypeVector x;
    std::size_t index;
    // Print info
    ijk = {uniform(e), uniform(e), uniform(e)};
    x = sMesh->positionPoint(ijk);
    index = sMesh->indexPoint(ijk);
    std::cout << "i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << "\nindex: " << index << "\nx: " << x.transpose() << std::endl;
    ijk = sMesh->ijkPoint(x.data());
    std::cout << "xReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    ijk = sMesh->ijkPoint(index);
    std::cout << "indexReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    std::cout << std::endl;
}

int main () {
    std::shared_ptr<TypeMesh> sMesh = std::make_shared<TypeMesh>(TypeContainer<std::size_t>(DIM, n), TypeContainer<double>(DIM, length), TypeVector::Constant(origin), TypeContainer<bool>(DIM, false));
    std::shared_ptr<TypeSub> sSubMesh = std::make_shared<TypeSub>(TypeContainer<std::size_t>(DIM, nSub), offSub, sMesh);
    // Random setup
    std::random_device r;
    std::default_random_engine e(r());
    std::uniform_int_distribution<int> uniform(-int(nSub), 2*nSub-1);
    // Print
    print(sSubMesh, uniform, e);
    print(sSubMesh, uniform, e);
    print(sSubMesh, uniform, e);
}
