#include "nm_em.h"
#include "eigenMedian.h"
#include "lm/fastLm.h"
#include "corr/eigenCorr.h"
#include "norm/normdens.h"

#include <climits>
#include <cfloat>
extern "C" {
#include "dcdflib/cdflib.h"
}
using namespace Eigen;

namespace PLNK {

	// check netstat chisq distrib working
	// R test:
	//	> pchisq(.285, 1, lower.tail=FALSE)
	//[1] 0.5934426
	//> 1-pchisq(.285, 1, lower.tail=FALSE)
	//[1] 0.4065574
	//> pchisq(.285, 1)
	//[1] 0.4065574
	//MyType chi = PLNK::pChisq(0.285, 1.0);

	double PLNK::pChisq(double x, double df)
	{

		//if ( ! realnum(x) ) return -9;

		double p, q;
		int st = 0; // error variable
		int w = 1; // function variable
		double bnd = 1; // boundary function

		// NCP is set to 0
		cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

		// Check status
		if (st != 0) return -9;

		// Return p-value
		//return q;
		return p;

	}

	double PLNK::pT(double T, double df)
	{

		//if ( ! realnum(T) ) return -9; 

		//T = fabs(T); ---> return 1-p

		double p, q;
		int st = 0;      // error variable
		int w = 1;       // function variable
		double bnd = 1;  // boundary function

		// NCP is set to 0
		cdft(&w, &p, &q, &T, &df, &st, &bnd);

		// Check status
		if (st != 0) return -9;

		// Return two-sided p-value
		//return 2*q;
		return p;

	}
}


namespace HOGGART {
	// implementation as described in paper: <ref>, with use of Eigen lib  
	int nm_hoggart(
		const Eigen::VectorXi &snp,
		const Eigen::VectorXd &y,
		const Eigen::VectorXd &snpAB,
		double &beta_hat,
		double &se_hat,
		double &p_hat) {

		int N = snp.size();
		MyType med[3];

		std::vector<MyType> groups[3];
		for (auto i = 0; i < N; i++) {
			int pos = snp(i);
			groups[pos].push_back(y(i));
		}
		for (auto i = 0; i < 3; i++) {
			med[i] = median(groups[i]);
		}

		Eigen::ArrayXd medY(y.size());
		for (auto i = 0; i < N; i++) {
			int pos = snp(i);
			medY(i) = med[pos];
		}
		medY = (y.array() - medY).abs();
		MyMatrix X(N, 2);
		X.setConstant(N, 2, 1.0);
		X.col(1) = snpAB;
		lmsol::lmres r1 = lmsol::fastLm(X, medY, lmsol::ColPivQR_t);
		VectorXd se = std::get<lmsol::se>(r1);
		VectorXd coef = std::get<lmsol::coefficients>(r1);

		beta_hat = coef(1);
		p_hat = 1.0;
		se_hat = se(1);
		if (fabs(se_hat) > 0.000000001)
			p_hat = PLNK::pNorm01(beta_hat / se_hat, false);
		return 1;
	}

	// wrapper function with standardization flag
	int nm_hoggart(
		const Eigen::VectorXi &snp,
		const Eigen::VectorXd &y,
		const Eigen::VectorXd &snpAB,
		const bool standardize,
		double &beta_hat,
		double &se_hat,
		double &p_hat) {

		if (standardize) {
			double meanY = y.mean();
			double stdY = stdDev<double>(y.array());
			return nm_hoggart(snp, (y.array() - meanY) / stdY, snpAB, beta_hat, se_hat, p_hat);
		}

		return nm_hoggart(snp, y, snpAB, beta_hat, se_hat, p_hat);
	}

	//-----------------------------------------------------------------------------------------------------
	// from here code ported or adapted from quicktest
	//		
#include "sortMatlab.h"

// test fuzzy_median
//double myy[] = { 100.0, 110.0, 9.0, 8.0, 10.0, 80.0, 90.0 };
//double myp[] = { 0.0,  0.0, 1.0, 1.0, 1.0 , 0.0, 0.0 };
//std::vector<double> y(myy, myy + sizeof(myy) / sizeof(double));
//std::vector<double> p(myp, myp + sizeof(myp) / sizeof(double));
//double med = HOGGART::fuzzy_median(y, p);

	double fuzzy_median(std::vector<double> y, std::vector<double> p) {

		int n = p.size();
		std::vector<double> sp(n);
		std::vector<double> csp(n);
		std::vector<double> icsp(n);
		std::vector<double> v(n);
		std::vector<double> w(n);
		//  vector<double> sy (n);
		//  vector<int> syix (n);

		std::vector<double> tmpy;
		tmpy.resize(n);
		for (int i = 0; i < n; ++i) {
			tmpy[i] = y[i];
		}
		std::vector<size_t> syix;
		std::vector<double> sy;

		sort(tmpy, sy, syix);
		for (int i = 0; i < n; ++i) { // order p and y
			sp[i] = p[syix[i]];
		}

		csp[0] = sp[0];
		for (int i = 1; i < n; ++i) { // create cumulative sum of sorted p
			csp[i] = csp[i - 1] + sp[i];
		}

		icsp[n - 1] = sp[n - 1];
		for (int i = n - 2; i >= 0; --i) { // create backward cumulative sum of sorted p
			icsp[i] = icsp[i + 1] + sp[i];
		}

		double mv = 2; // can only be 1 at max
		for (int i = 0; i < n; ++i) {
			if (sp[i] == 0) {
				v[i] = 2;
			}
			else {
				v[i] = fabs(csp[i] - icsp[i]);
			}
			mv = std::min(mv, v[i]);
		}

		double sw = 0;
		for (int i = 0; i < n; ++i) {
			if (v[i] == mv) {
				w[i] = 1;
				sw += 1;
			}
			else {
				w[i] = 0;
			}
		}

		double med = 0;
		for (int i = 0; i < n; ++i) {
			// double tmp = w[i]/sw;
			// cout << i << " " << sy[i] << " " << tmp << endl;
			med += sy[i] * w[i] / sw;
		}

		return(med);
	}

	double fuzzy_mean(std::vector<double> y, std::vector<double> p) {

		int n = p.size();

		double sumyp = 0;
		double sump = 0;
		for (int i = 0; i < n; ++i) {
			sumyp += p[i] * y[i];
			sump += p[i];
		}

		double med = sumyp / sump;

		return(med);
	}

	double fuzzy_var(std::vector<double> y, std::vector<double> p) {

		int n = p.size();

		double sumyp = 0;
		double sump = 0;
		for (int i = 0; i < n; ++i) {
			sumyp += p[i] * y[i];
			sump += p[i];
		}

		double med = sumyp / sump;

		double sumy2p = 0;
		for (int i = 0; i < n; ++i) {
			sumy2p += p[i] * pow(y[i] - med, 2);
		}

		double var = sumy2p / (sump - 1);

		return(var);
	}

	// implementation derived from quicktest
	int quicktest_hoggart(std::vector<double> pAA, std::vector<double> pAB, std::vector<double> pBB, std::vector<double> y,
		double n, double &beta_hat, double &se_hat, double &p_hat, double &varAA, double &varAB, double &varBB) {

		int nInt = int(n);
		double df = nInt - 1;
		double meanY = 0;
		for (int i = 0; i < n; ++i) {
			meanY += y[i] / n;
		}
		double stdY = 0;
		for (int i = 0; i < n; ++i) {
			stdY += pow(y[i] - meanY, 2);
		}
		stdY = sqrt(stdY / (n - 1));
		for (int i = 0; i < n; ++i) {
			y[i] = (y[i] - meanY) / stdY;
		}

		double yBB = fuzzy_median(y, pBB);
		double yAB = fuzzy_median(y, pAB);
		double yAA = fuzzy_median(y, pAA);

		varAA = fuzzy_var(y, pAA);
		varAB = fuzzy_var(y, pAB);
		varBB = fuzzy_var(y, pBB);

		std::vector<double> y0(nInt);
		std::vector<double> y0sq(nInt);
		double my;
		double sum_y0sq = 0;
		for (int i = 0; i < n; ++i) {
			my = yAA * pAA[i] + yAB * pAB[i] + yBB * pBB[i];
			y0[i] = fabs(y[i] - my);
			y0sq[i] = pow(y0[i], 2);
			sum_y0sq += y0sq[i];
		}

		double n0 = 0;
		double n1 = 0;
		double n2 = 0;
		double z02 = 0;
		double z1 = 0;
		double zz = 0;

		for (int i = 0; i < n; ++i) {
			n0 += pAA[i];
			n1 += pAB[i];
			n2 += pBB[i];
		}
		for (int i = 0; i < n; ++i) {
			z02 += (pAA[i] + pBB[i])*y0[i] / (n0 + n2);
			z1 += pAB[i] * y0[i] / n1;
			zz += y0[i] / n;
		}

		double w_bar = n1 / n;
		double yTg = -(n0 + n2)*w_bar*(z02 - zz) + n1 * (1 - w_bar)*(z1 - zz);
		double gTg = (n1 - n * pow(w_bar, 2));

		beta_hat = yTg / gTg;
		double rss = (sum_y0sq - n * pow(zz, 2)) - pow(beta_hat, 2)*gTg;
		se_hat = sqrt((rss / (n - 1)) / gTg);
		//double p_hat = pt(-beta_hat/se_hat, df, 1, 0);
		p_hat = PLNK::pT(-beta_hat / se_hat, df);

		//  beta = beta_hat;
		//  se = se_hat;
		//  pval = p_hat;

		return (1);
	}

	// implementation derived from quicktest with use of Eigen::lib
	int eigen_hoggart(const Eigen::ArrayXd &pAA, const Eigen::ArrayXd &pAB, const Eigen::ArrayXd &pBB, const Eigen::ArrayXd &yy,
		double &beta_hat, double &se_hat, double &p_hat, double &varAA, double &varAB, double &varBB,
		const bool standardize) {
		const int AA = 0;
		const int AB = 1;
		const int BB = 2;

		Eigen::ArrayXd y = yy;
		double n = double(yy.size());
		int nInt = int(n);
		double df = nInt - 1;

		if (standardize) {
			// set new y centered on 0
			double meanY = yy.sum() / n;
			double stdY = (yy - meanY).square().sum();
			stdY = sqrt(stdY / df);
			y -= meanY;
			y /= stdY;
		}
		//
		std::vector<MyType> groups[3];
		// groups[].reserve(n);
		for (auto i = 0; i < n; i++) {
			int pos = -1;
			if (pAA[i])
				pos = AA;
			else if (pAB[i])
				pos = AB;
			else if (pBB[i])
				pos = BB;

			if (pos != -1) {
				groups[pos].push_back(y(i));
			}
		}

		double yAA = median(groups[AA]); //fuzzy_median(y, pBB);
		double yAB = median(groups[AB]); //fuzzy_median(y, pAB);
		double yBB = median(groups[BB]); //fuzzy_median(y, pAA);

		varAA = group_var(y, pAA);
		varAB = group_var(y, pAB);
		varBB = group_var(y, pBB);

		Eigen::ArrayXd y0(nInt);   // Eigen::ArrayXd y0sq (nInt);
		double sum_y0sq = 0;
		Eigen::ArrayXd my = (yAA*pAA + yAB * pAB + yBB * pBB);
		y0 = (y - my).abs();						// y0sq = y0.square();
		sum_y0sq = y0.square().sum();

		double n0 = pAA.sum();	//groups[AA].size()
		double n1 = pAB.sum();
		double n2 = pBB.sum();
		double z02 = ((pAA + pBB)*y0).sum() / (n0 + n2);
		double z1 = (pAB*y0).sum() / n1;
		double zz = y0.sum() / n;

		double w_bar = n1 / n;
		double yTg = -(n0 + n2)*w_bar*(z02 - zz) + n1 * (1 - w_bar)*(z1 - zz);
		double gTg = (n1 - n * pow(w_bar, 2));

		beta_hat = yTg / gTg;
		double rss = (sum_y0sq - n * pow(zz, 2)) - pow(beta_hat, 2)*gTg;
		se_hat = sqrt((rss / (n - 1)) / gTg);
		//p_hat = pt(-beta_hat/se_hat, df, 1, 0);
		p_hat = PLNK::pT(-beta_hat / se_hat, df);

		return (1);
	}

} // hoggart namespace end