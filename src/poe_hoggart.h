#ifndef ___EM___
#define ___EM___

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Eigen/Core"
#include "Eigen/Dense"

namespace PLNK {
	double pChisq(double x, double df);
	double pT(double x, double df);
}

namespace HOGGART {
	
	// chiara, from R code
	// with use of fastlm std pheno or not
	int nm_hoggart(
		const Eigen::VectorXi &snp,
		const Eigen::VectorXd &y,
		const Eigen::VectorXd &snpAB,
		const bool standardize,
		double &beta_hat,
		double &se_hat,
		double &p_hat);
		
	// from quicktest impl.
	// always standardize phenotype
	int quicktest_hoggard (
		std::vector<double> pAA,
		std::vector<double> pAB,
		std::vector<double> pBB,
		std::vector<double> y,
		double n, double &beta, double &se, double &pval, double &varAA, double &varAB, double &varBB);
	
	// from quicktest impl. with use of eigen::array
	// by default standardize phenotype
	int eigen_hoggard (
		const Eigen::ArrayXd &pAA,
		const Eigen::ArrayXd &pAB,
		const Eigen::ArrayXd &pBB,
		const Eigen::ArrayXd &y,
		double &beta_hat, double &se_hat, double &p_hat, double &varAA, double &varAB, double &varBB,
		const bool standardize=true);

	double fuzzy_median (std::vector<double> y, std::vector<double> p);	
}


#endif