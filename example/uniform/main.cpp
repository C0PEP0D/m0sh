// std includes
#include <iostream>
#include <vector>
#include <random>
// thirdparties includes
#include <Eigen/Dense>
// lib includes
#include "m0sh/uniform.h"

using tScalar = double;
// Space
const unsigned int DIM = 3;
using tSpaceVector = Eigen::Matrix<tScalar, DIM, 1>;
template<typename ...Args> using tView = Eigen::Map<Args...>;
// Mesh
using tMesh = m0sh::Uniform<tSpaceVector, tView, DIM>;
// Random setup
const int n = 100;
std::random_device r;
std::default_random_engine e(r());
std::uniform_int_distribution<int> uniform(-n, 2*n);

void print(const double* pOrigin, const double* pSpacing, const unsigned int* pNbPointsPerAxis) {
    std::vector<int> ijk = {uniform(e), uniform(e), uniform(e)};
    tSpaceVector x = tMesh::positionPoint(pOrigin, pSpacing, ijk.data());
    unsigned int index = tMesh::indexPoint(pNbPointsPerAxis, ijk.data());
    // Print info
    std::cout << "i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << "\nindex: " << index << "\nx: " << x.transpose() << std::endl;
    ijk = tMesh::ijkCell(pOrigin, pSpacing, x.data());
    std::cout << "xReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    ijk = tMesh::ijkPoint(pNbPointsPerAxis, index);
    std::cout << "indexReverse: " << " i: " << ijk[0] << " j: " << ijk[1] << " k: " << ijk[2] << std::endl;
    std::cout << std::endl;
}

int main () {
	// Mesh
	const tSpaceVector origin = tSpaceVector::Constant(0.0);
	const tSpaceVector spacing = tSpaceVector::Constant(1.0);
	const std::vector<unsigned int> nbPointsPerAxis = {n, n, n};
    // Print
    print(origin.data(), spacing.data(), nbPointsPerAxis.data());
    print(origin.data(), spacing.data(), nbPointsPerAxis.data());
    print(origin.data(), spacing.data(), nbPointsPerAxis.data());
}
