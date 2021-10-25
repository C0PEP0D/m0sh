// std includes
#include <iostream>
#include <vector>
#include <random>
// thirdparties includes
#include <Eigen/Dense>
// lib includes
#include "m0sh/regular.h"

using TypeScalar = double;
// Space
const unsigned int DIM = 3;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
// Mesh
template<typename ...Args>
using TypeContainer = std::vector<Args...>;
using TypeMesh = m0sh::Regular<TypeVector, TypeRef, TypeContainer>;
// Data
const int n = 100;
const double length = 1.0;
const double origin = 0.5;

void print(const TypeMesh& mesh, std::uniform_int_distribution<int>& uniform, std::default_random_engine& e) {
    TypeContainer<int> ijk;
    TypeVector x;
    std::size_t index;
    // Print info
    ijk = {uniform(e), uniform(e), uniform(e)};
    x = mesh.x(ijk);
    index = mesh.index(ijk);
    std::cout << "i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << " index: " << index << " x: \n" << x << std::endl;
    ijk = mesh.ijk(x);
    std::cout << "xReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    ijk = mesh.ijk(index);
    std::cout << "indexReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    std::cout << std::endl;
}

int main () { 
    TypeMesh mesh(TypeContainer<std::size_t>(DIM, n), TypeContainer<double>(DIM, length), TypeVector::Constant(origin));
    // Random setup
    std::random_device r;
    std::default_random_engine e(r());
    std::uniform_int_distribution<int> uniform(-n, 2*n);
    // Print
    print(mesh, uniform, e);
    print(mesh, uniform, e);
    print(mesh, uniform, e);
}
