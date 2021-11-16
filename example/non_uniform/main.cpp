// std includes
#include <iostream> // cout, endl
#include <vector>
#include <random> // random_device, default_random_engine, uniform_int_distribution
#include <algorithm> // generate
// thirdparties includes
#include <Eigen/Dense>
// lib includes
#include "m0sh/non_uniform.h"

using TypeScalar = double;
// Space
const unsigned int DIM = 3;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
// Mesh
template<typename ...Args>
using TypeContainer = std::vector<Args...>;
using TypeMesh = m0sh::NonUniform<TypeVector, TypeRef, TypeContainer>;
// Data
const std::size_t n = 10;
const double dx = 0.1;

void print(const TypeMesh& mesh, std::uniform_int_distribution<int>& uniform, std::default_random_engine& e) {
    TypeContainer<int> ijk;
    TypeVector x;
    std::size_t index;
    // Print info
    ijk = {uniform(e), uniform(e), uniform(e)};
    x = mesh.positionCell(ijk);
    index = mesh.indexCell(ijk);
    std::cout << "i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << "\nindex: " << index << "\nx: " << x.transpose() << std::endl;
    ijk = mesh.ijkCell(x);
    std::cout << "xReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    ijk = mesh.ijkCell(index);
    std::cout << "indexReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    std::cout << std::endl;
}

int main () {
    // Build axis
    TypeContainer<double> axis(n+1);
    std::generate(axis.begin(), axis.end(), [x = 0.0] () mutable { return x += dx; });
    // Build grid
    TypeContainer<TypeContainer<double>> grid(DIM, axis);
    // Build mesh finally
    TypeMesh mesh(grid, TypeContainer<bool>(true));
    // Random setup
    std::random_device r;
    std::default_random_engine e(r());
    std::uniform_int_distribution<int> uniform(0, n);
    // Print
    print(mesh, uniform, e);
    print(mesh, uniform, e);
    print(mesh, uniform, e);
}
