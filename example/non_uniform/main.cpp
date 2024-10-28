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
using tSpaceVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args> using tView = Eigen::Map<Args...>;
// Mesh
using tMesh = m0sh::NonUniform<tSpaceVector, tView, DIM>;
// Mesh parameters
const std::size_t n = 10;
const double dx = 1.0;
// Random setup
std::random_device r;
std::default_random_engine e(r());
std::uniform_int_distribution<int> uniform(0, n);

void print(const std::vector<double>* pAxes) {
    std::vector<int> ijk = {uniform(e), uniform(e), uniform(e)};
    const tSpaceVector x = tMesh::positionPoint(pAxes, ijk.data());
    std::size_t index = tMesh::indexPoint(pAxes, ijk.data());
    // Print info
    std::cout << "i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << "\nindex: " << index << "\nx: " << x.transpose() << std::endl;
    ijk = tMesh::ijkCell(pAxes, x.data());
    std::cout << "xReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    ijk = tMesh::ijkPoint(pAxes, index);
    std::cout << "indexReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    std::cout << std::endl;
}

int main () {
    // Build axis
    std::vector<double> axis(n+1);
    std::generate(axis.begin(), axis.end(), [x = -dx] () mutable { return x += dx; });
    // Build grid
    std::vector<std::vector<double>> axes(DIM, axis);
    // Print
    print(axes.data());
    print(axes.data());
    print(axes.data());
}
