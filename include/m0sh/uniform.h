#ifndef M0SH_UNIFORM_H
#define M0SH_UNIFORM_H
#pragma once

// std includes
#include <cmath>
// module includes
#include "m0sh/structured.h"

namespace m0sh {

template<typename _tSpaceVector, template<typename...> class _tView, unsigned int _Dim>
class Uniform : public Structured<_tSpaceVector, _tView, _Dim> {
    public:
        using tBase = Structured<_tSpaceVector, _tView, _Dim>;
   		using tSpaceVector = typename tBase::tSpaceVector;
   		template<typename... Args> using tView = typename tBase::tView<Args...>;
   	public:
   		using tBase::Dim;
    public:
    	Uniform(){
    		
    	};
   	public:
    	static unsigned int indexPoint(const unsigned int* pNbPointsPerAxis, const int* pIjk) {
            return tBase::index(pNbPointsPerAxis, pIjk);
        };
        static unsigned int indexCell(const unsigned int* pNbPointsPerAxis, const int* pIjk) {
			const std::vector<unsigned int> _nbCellsPerAxis = tBase::nbCellsPerAxis(pNbPointsPerAxis);
           	return tBase::index(_nbCellsPerAxis.data(), pIjk);
       	};
    public:
        static tSpaceVector positionPoint(const double* pOrigin, const double* pSpacing, const int* pIjk) {
            tSpaceVector position;
            for(std::size_t i = 0; i < Dim; i++) {
                position[i] = pOrigin[i] + pIjk[i] * pSpacing[i];
            }
            return position;
        };
        
        static tSpaceVector positionCell(const double* pOrigin, const double* pSpacing, const int* pIjk) {
            tSpaceVector position;
            for(std::size_t i = 0; i < Dim; i++) {
                position[i] = pOrigin[i] + pIjk[i] * pSpacing[i] + 0.5 * pSpacing[i];
            }
            return position;
        };

        static tSpaceVector positionPoint(const double* pOrigin, const double* pSpacing, const unsigned int* pNbPointsPerAxis, const unsigned int p_index) {
            return positionPoint(pOrigin, pSpacing, tBase::ijk(pNbPointsPerAxis, p_index));
        };
        
        static tSpaceVector positionCell(const double* pOrigin, const double* pSpacing, const unsigned int* pNbPointsPerAxis, const unsigned int p_index) {
        	std::vector<int> nbCellsPerAxis = tBase::nbCellsPerAxis(pNbPointsPerAxis);
            return positionCell(pOrigin, pSpacing, tBase::ijk(nbCellsPerAxis.data(), p_index));
        };
        
   public:
        static std::vector<int> ijkCell(const double* pOrigin, const double* pSpacing, const double* pPosition) {
            std::vector<int> ijk_(Dim);
            for(std::size_t i = 0; i < Dim; i++) {
                ijk_[i] = std::floor((pPosition[i] - pOrigin[i]) / pSpacing[i]);
            }
            return ijk_;
        };

        static unsigned int indexCell(const double* pOrigin, const double* pSpacing, const unsigned int* pNbPointsPerAxis, const double* pPosition) {
        	const std::vector<int> _ijkCell = ijkCell(pOrigin, pSpacing, pPosition);
        	const std::vector<unsigned int> _nbCellsPerAxis = tBase::nbCellsPerAxis(pNbPointsPerAxis);
            return tBase::index(_nbCellsPerAxis.data(), _ijkCell.data());
        };
};

}

#endif
