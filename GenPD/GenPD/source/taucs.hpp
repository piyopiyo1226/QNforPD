#pragma once
#include <iostream>
#include <vector>
#include <ccomplex>
//included by the simulation
#include "Sparce.hpp"

//namespace arcsim {
	

	//----------------------------------------------------------

	std::vector<double> taucs_linear_solve(const arcsim::SpMat<double> &A, const std::vector<double> &b);

	template <int m>
	std::vector<arcsim::Vec<m> > taucs_linear_solve(const arcsim::SpMat<arcsim::Mat<m, m> > &A, const std::vector<arcsim::Vec<m> > &b);


	////-------------------------------------------1.2------------------------------------------------
	//taucs_ccs_matrix *sparse_to_taucs(const arcsim::SpMat<double> &As);
	////------------------------------------------1,1--------------------------------------------

	//std::vector<double> alglib_linear_solve(const arcsim::SpMat<double> &A, const std::vector<double> &b);


	
//
//}
//-----------------------------------------------------------------------------------------------


