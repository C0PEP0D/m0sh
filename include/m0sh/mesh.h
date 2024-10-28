#ifndef M0SH_MESH_H
#define M0SH_MESH_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>

namespace m0sh {

template<typename _tSpaceVector, template<typename...> class _tView>
class Mesh {
	public:
		using tSpaceVector = _tSpaceVector;
		template<typename... Args> using tView = _tView<Args...>;
    public:
        Mesh() {
        };
    public:
        // individual
        // virtual static tSpaceVector positionPoint(const std::size_t& index) = 0;
        // virtual static tSpaceVector positionCell(const std::size_t& index) = 0;
        // virtual static std::size_t indexPoint(const double* pPosition) = 0;
        // virtual static std::size_t indexCell(const double* pPosition) = 0;
        // nb
        // virtual std::size_t nbPoints() const = 0;
        // virtual std::size_t nbCells() const = 0;
        // all
        static std::vector<unsigned int> indexs(unsigned int nb) {
            std::vector<unsigned int> indexs(nb);
            std::iota(indexs.begin(), indexs.end(), 0);
            return indexs;
        };
        // static std::vector<tSpaceVector> positionPoints() {
            // std::vector<tSpaceVector> positions(nbPoints());
            // for(unsigned int i = 0; i < nbPoints(); i++) {
                // positions[i] = positionPoint(i);
            // }
            // return positions;
        // };
        // static std::vector<tSpaceVector> positionCells() {
            // std::vector<tSpaceVector> positions(nbCells());
            // for(unsigned int i = 0; i < nbCells(); i++) {
                // positions[i] = positionCell(i);
            // }
            // return positions;
        // };
};

}

#endif
