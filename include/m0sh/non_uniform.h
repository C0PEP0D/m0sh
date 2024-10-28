#ifndef M0SH_NON_UNIFORM_H
#define M0SH_NON_UNIFORM_H
#pragma once

// std includes
#include <cmath>
#include <vector>
// module includes
#include "m0sh/structured.h"

namespace m0sh {

template<typename _tSpaceVector, template<typename...> class _tView, unsigned int _Dim>
class NonUniform : public Structured<_tSpaceVector, _tView, _Dim> {
	public:
		using tBase = Structured<_tSpaceVector, _tView, _Dim>;
		using tSpaceVector = typename tBase::tSpaceVector;
		template<typename... Args> using tView = typename tBase::tView<Args...>;
	public:
		using tBase::Dim;
	public:
		NonUniform(){
			
		};
	public:
		static std::vector<unsigned int> nbPointsPerAxis(const std::vector<double>* pAxesPoints) {
			std::vector<unsigned int> output(Dim);
			for(unsigned int i = 0; i < Dim; i++) {
				output[i] = pAxesPoints[i].size();
			}
            return output;
        };
    public:
    	static std::vector<int> ijkPoint(const std::vector<double>* pAxesPoints, const unsigned int index) {
			std::vector<unsigned int> _nbPointsPerAxis = nbPointsPerAxis(pAxesPoints);
			return tBase::ijk(_nbPointsPerAxis.data(), index);
        };
        static std::vector<int> ijkPointPeriodic(const std::vector<double>* pAxesPoints, const unsigned int index) {
			std::vector<unsigned int> _nbPointsPerAxis = nbPointsPerAxis(pAxesPoints);
			std::vector<unsigned int> nbCellsPerAxis = tBase::nbCellsPerAxis(_nbPointsPerAxis.data());
			return tBase::ijk(nbCellsPerAxis.data(), index);
        };
        static std::vector<int> ijkCell(const std::vector<double>* pAxesPoints, const unsigned int index) {
			return ijkPointPeriodic(pAxesPoints, index);
        };
    public:
		static tSpaceVector positionPoint(const std::vector<double>* pAxesPoints, const int* pIjk) {
			// Compute ijkPeriodic
			std::vector<int> ijkPeriodic(Dim);
			std::vector<int> ijkDiv(Dim);
			for(std::size_t i = 0; i < Dim; i++) {
				ijkPeriodic[i] = std::abs(pIjk[i]) % (pAxesPoints[i].size() - 1);
				ijkDiv[i] = std::abs(pIjk[i]) / (pAxesPoints[i].size() - 1);
				if(pIjk[i] < 0) {
					ijkPeriodic[i] = pAxesPoints[i].size() - 1 - ijkPeriodic[i];
					ijkDiv[i] = -ijkDiv[i] - 1;
				}
			}
			// Compute position
			tSpaceVector position;
			for(std::size_t i = 0; i < Dim; i++) {
				const double axisLength = pAxesPoints[i][pAxesPoints[i].size() - 1] - pAxesPoints[i][0];
				position[i] = pAxesPoints[i][ijkPeriodic[i]] + ijkDiv[i] * axisLength;
			}
			return position;
		};
		
		static tSpaceVector positionCell(const std::vector<double>* pAxesPoints, const int* pIjk) {
			std::vector<int> ijkPlusOne(Dim);
			for(std::size_t i = 0; i < Dim; i++) {
				ijkPlusOne[i] = pIjk[i] + 1;
			}
			const tSpaceVector positionPointIjk = positionPoint(pAxesPoints, pIjk);
			const tSpaceVector positionPointIjkPlusOne = positionPoint(pAxesPoints, ijkPlusOne.data());
			return 0.5 * (positionPointIjk + positionPointIjkPlusOne);
		};

		static tSpaceVector positionPoint(const std::vector<double>* pAxesPoints, const unsigned int p_index) {
			const std::vector<unsigned int> _nbPointsPerAxis = nbPointsPerAxis(pAxesPoints);
			const std::vector<int> _ijk = tBase::ijk(_nbPointsPerAxis.data(), p_index);
			return positionPoint(pAxesPoints, _ijk.data());
		};

		static tSpaceVector positionPointPeriodic(const std::vector<double>* pAxesPoints, const unsigned int p_index) {
			const std::vector<unsigned int> _nbPointsPerAxis = nbPointsPerAxis(pAxesPoints);
			const std::vector<unsigned int> _nbCellsPerAxis = tBase::nbCellsPerAxis(_nbPointsPerAxis.data());
			const std::vector<int> _ijk = tBase::ijk(_nbCellsPerAxis.data(), p_index);
			return positionPoint(pAxesPoints, _ijk.data());
		};

		static tSpaceVector positionCell(const std::vector<double>* pAxesPoints, const unsigned int p_index) {
			const std::vector<unsigned int> _nbPointsPerAxis = nbPointsPerAxis(pAxesPoints);
			const std::vector<unsigned int> _nbCellsPerAxis = tBase::nbCellsPerAxis(_nbPointsPerAxis.data());
			const std::vector<int> _ijk = tBase::ijk(_nbCellsPerAxis.data(), p_index);
			return positionCell(pAxesPoints, _ijk.data());
		};

	public:
		static std::vector<int> ijkCell(const std::vector<double>* pAxesPoints, const double* pPosition) {
			// Compute positionPeriodic
			tSpaceVector positionPeriodic;
			std::vector<int> positionDiv(Dim);
			for(std::size_t i = 0; i < Dim; i++) {
				const double axisOrigin = pAxesPoints[i][0];
				const double axisLength = pAxesPoints[i][pAxesPoints[i].size() - 1] - pAxesPoints[i][0];
				positionPeriodic[i] = std::fmod(std::abs(pPosition[i] - axisOrigin), axisLength);
				positionDiv[i] = std::floor(std::abs(pPosition[i] - axisOrigin) / axisLength);
				if(pPosition[i] < axisOrigin) {
					positionPeriodic[i] = axisLength - positionPeriodic[i];
					positionDiv[i] = -positionDiv[i] - 1;
				}
				positionPeriodic[i] += axisOrigin;
			}
			// Compute index
			std::vector<int> ijk_(Dim);
			for(std::size_t i = 0; i < Dim; i++) {
				ijk_[i] = std::distance(pAxesPoints[i].begin(), std::upper_bound(pAxesPoints[i].begin(), pAxesPoints[i].end(), positionPeriodic[i])) - 1 + positionDiv[i] * (pAxesPoints[i].size() - 1);
			}
			return ijk_;
		};

		static unsigned int indexCell(const std::vector<double>* pAxesPoints, const double* pPosition) {
        	const std::vector<unsigned int> _nbPointsPerAxis = nbPointsPerAxis(pAxesPoints);
        	const std::vector<int> _ijkCell = ijkCell(pAxesPoints, pPosition);
            return tBase::indexCell(_nbPointsPerAxis.data(), _ijkCell.data());
        };
};

}

#endif
