#ifndef M0SH_STRUCTURED_H
#define M0SH_STRUCTURED_H
#pragma once

// std includes
#include <cmath> // div, ldiv_t
#include <numeric> // accumulate
#include <functional> // multiplies
#include <vector>
#include <cassert>
// module include
#include "m0sh/mesh.h"

namespace m0sh {

template<typename _tSpaceVector, template<typename...> class _tView, unsigned int _Dim>
class Structured : public Mesh<_tSpaceVector, _tView> {
	public:
		using tBase = Mesh<_tSpaceVector, _tView>;
	public:
		using tSpaceVector = typename tBase::tSpaceVector;
		template<typename... Args> using tView = typename tBase::tView<Args...>;
	public:
		constexpr static unsigned int Dim = _Dim;
    public:
        Structured() {
        };
    public:
		static std::vector<unsigned int> nbCellsPerAxis(const unsigned int* pNbPointsPerAxis) {
			std::vector<unsigned int> output(Dim);
			for(unsigned int i = 0; i < Dim; i++) {
				output[i] = pNbPointsPerAxis[i] - 1;
			}
            return output;
        };
    
        static unsigned int nb(const unsigned int* pNbPerAxis) {
            return std::accumulate(pNbPerAxis, pNbPerAxis + Dim, 1.0, std::multiplies<double>());
        };

        static std::vector<unsigned int> indexs(const unsigned int* pNbPerAxis) {
            return tBase::indexs(nb(pNbPerAxis));
        };

        static unsigned int indexPoint(const unsigned int* pNbPointsPerAxis, const int* pIjk) {
            // Compute index
            unsigned int nTot = 1;
            unsigned int index = pIjk[0];
           	assert((pIjk[0] > -1 && pIjk[0] < pNbPointsPerAxis[0]) && "ijk is not in range");
            for(unsigned int i = 1; i < Dim; i++) {
            	assert((pIjk[i] > -1 && pIjk[i] < pNbPointsPerAxis[i]) && "ijk is not in range");
                nTot *= pNbPointsPerAxis[i-1];
                index += nTot * pIjk[i];
            }
            return index;
        };
        
        static unsigned int indexPointPeriodic(const unsigned int* pNbPointsPerAxis, const int* pIjk) {
            // Compute ijk periodic
            std::vector<int> ijkPeriodic(Dim);
            for(unsigned int i = 0; i < Dim; i++) {
                ijkPeriodic[i] = std::abs(pIjk[i]) % (pNbPointsPerAxis[i] - 1);
                if(pIjk[i] < 0 && ijkPeriodic[i] != 0) {
                    ijkPeriodic[i] = (pNbPointsPerAxis[i] - 1) - ijkPeriodic[i];
                }
            }
            // Compute index
            unsigned int nTot = 1;
            unsigned int index = ijkPeriodic[0];
            for(unsigned int i = 1; i < Dim; i++) {
                nTot *= (pNbPointsPerAxis[i-1] - 1);
                index += nTot * ijkPeriodic[i];
            }
            return index;
        };

        static unsigned int indexCell(const unsigned int* pNbPointsPerAxis, const int* pIjk) {
            return indexPointPeriodic(pNbPointsPerAxis, pIjk);
        };
        
        static std::vector<int> ijk(const unsigned int* pNbPerAxis, const unsigned int index) {
            // Compute nDiv
            std::vector<unsigned int> nDiv(Dim, 1);
            for(unsigned int k = 1; k < Dim; k++) {
                nDiv[k] = nDiv[k-1] * pNbPerAxis[k-1];
            }
            // Compute ijk
            std::ldiv_t dv = {0, (long int)index};
            std::vector<int> ijk_(Dim);
            for(unsigned int k = 0; k < Dim; k++) {
                dv = std::div(dv.rem, nDiv[Dim - k - 1]);
                ijk_[Dim - k - 1] = dv.quot;
            }
            return ijk_;
        };

        static std::vector<unsigned int> subIndexsPoints(const unsigned int* pNbPointsPerAxis, const unsigned int* pSubNbPointsPerAxis, const int* pOffset) {
			std::vector<unsigned int> output(nb(pSubNbPointsPerAxis));
            for (const unsigned int subIndex : indexs(pSubNbPointsPerAxis)) {
               	std::vector<int> _ijk = ijk(pSubNbPointsPerAxis, subIndex);
               	for(unsigned int i = 0; i < Dim; i++){
               		_ijk[i] += pOffset[i];
               	}
               	output[subIndex] = indexPoint(pNbPointsPerAxis, _ijk.data());
            }
            return output;
        };

        static std::vector<unsigned int> subIndexsPointsPeriodic(const unsigned int* pNbPointsPerAxis, const unsigned int* pSubNbPointsPerAxis, const int* pOffset) {
			std::vector<unsigned int> output(nb(pSubNbPointsPerAxis));
            for (const unsigned int subIndex : indexs(pSubNbPointsPerAxis)) {
               	std::vector<int> _ijk = ijk(pSubNbPointsPerAxis, subIndex);
               	for(unsigned int i = 0; i < Dim; i++){
               		_ijk[i] += pOffset[i];
               	}
               	output[subIndex] = indexPointPeriodic(pNbPointsPerAxis, _ijk.data());
            }
            return output;
        };

        static std::vector<unsigned int> subIndexsCells(const unsigned int* pNbPointsPerAxis, const unsigned int* pSubNbPointsPerAxis, const int* pOffset) {
            return subIndexsPointsPeriodic(pNbPointsPerAxis, pSubNbPointsPerAxis, pOffset);
        };
        
        // virtual tSpaceVector positionPoint(const tContainer<int>& ijk) const = 0;
        // virtual tSpaceVector positionCell(const tContainer<int>& ijk) const = 0;
        // tSpaceVector positionPoint(const unsigned int& index) const override {
            // return positionPoint(ijkPoint(index));
        // };
        // tSpaceVector positionCell(const unsigned int& index) const override {
            // return positionCell(ijkCell(index));
        // };
};

}

#endif
