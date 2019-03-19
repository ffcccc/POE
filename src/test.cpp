#include <iostream>
// data loaders
#include "fast-cpp-csv-parser-master/csv.h"
// EM
#include "nm_em.h"
#include "numDeriv/numderiv.h"
#include "lm/fastLm.h"
#include "eigenMedian.h"
#include "corr/eigenCorr.h"
#include "norm/normdens.h"
#include "time/timer11.h"


int testtime()
{
    Timer tmr;
    double t = tmr.elapsed();
    std::cout << t << std::endl;

    tmr.reset();
    t = tmr.elapsed();
    std::cout << t << std::endl;
    return 0;
}

using namespace std;

// 
int main()
{
	// read from csv
	int N = 1;
	int NOrderMagn = 1;
	int NTimes = 1;
	int NDataset = 1;
	std::vector<double> gY (N, 0.0); 
	std::vector<double> gAA(N, 0.0);
	std::vector<double> gAB(N, 0.0);
	std::vector<double> gBB(N, 0.0);
	std::vector<double> gSNPd(N, 0);
	std::vector<int> gSNP(N, 0);
	
	// test poe 
	for(int k=0; k<NDataset; k++){
		// file 1k
		string fname("../../data/dataset_1_1_1000");
		N = 1000;
		
		// file 100k
		if(k==1){
			fname = "../../data/dataset_1_1_100000";
			N = 100000;
		}

		// read 4-column csv
		io::CSVReader<4> in(fname);
		gY.resize (N, 0.0); 
		gAA.resize(N, 0.0);
		gAB.resize(N, 0.0);
		gBB.resize(N, 0.0);
		gSNP.resize(N, 0.0);
		gSNPd.resize(N, 0.0);
		
		in.read_header(io::ignore_extra_column, "snp", "w_ab", "w_bb", "y");
		double val_snp;
		double val_bb;
		double val_ab;
		double val_y;
		int linecount(0);
		while (in.read_row(val_snp, val_ab, val_bb, val_y))
		{
			// do stuff with the data
			gY [linecount] = val_y;
			gAB[linecount]= val_ab;
			gBB[linecount]= val_bb;
			gSNPd[linecount] = val_snp;
			gSNP[linecount] = int(val_snp);
			if(gSNP[linecount] == 0){
				// not present in dataset
				gAA[linecount] = 1.0;
			}
			linecount++;
		}
		Eigen::Map<Eigen::VectorXd> parY(gY.data(), gY.size());
		Eigen::Map<Eigen::VectorXd> parAA(gAA.data(), gAA.size());
		Eigen::Map<Eigen::VectorXd> parAB(gAB.data(), gAB.size());
		Eigen::Map<Eigen::VectorXd> parBB(gBB.data(), gBB.size());
		Eigen::Map<Eigen::VectorXi> parSNP(gSNP.data(), gSNP.size());
		Eigen::Map<Eigen::VectorXd> parSNPd(gSNPd.data(), gSNPd.size());

		double avgAA = group_mean<double>(parY.array(), parAA.array());
		double avgBB = group_mean<double>(parY.array(), parBB.array());

		ValArr vv(4), vv2(4);    
		MyType aic = 0.0;
		
		cout << "Pop. size: " << N << endl;
		cout << "# Hogg test | POE test | POE + Emp. Var. | POE + [Emp.+ Obs.] Var." << endl;
		MyMatrix empVar, empVar2, obsVar, empVarInv, obsVarInv;
		MyMatrix obsVar0, obsVar1, obsVar2, obsVarInv0, obsVarInv1, obsVarInv2;
		
		for(auto j=0; j<NOrderMagn; j++){
			NTimes = pow(10,j);
			cout << NTimes;
			
			Timer tmrHogg;
			for(int i=0; i<NTimes; i++){
				double beta_hat, se_hat, p_hat, varAA, varAB, varBB;
				// hoggart test as described in paper - no standardization
				HOGGART::nm_hoggart(parSNP, parY, parAB, false, beta_hat, se_hat, p_hat);
				MyType res = beta_hat;
				MyType pres = p_hat;
				// hoggart test as described in paper - with standardization
				HOGGART::nm_hoggart(parSNP, parY, parAB, true, beta_hat, se_hat, p_hat);
				MyType res2 = beta_hat;
				MyType pres2 = p_hat;

				// compare hoggart test from quicktest with use of eigen vectors and standardization
				HOGGART::eigen_hoggard(parAA, parAB, parBB, parY, beta_hat, se_hat, p_hat, varAA, varAB, varBB, true);
				double hogg2 = beta_hat;
				double phogg2 = p_hat;
				// hoggart test from quicktest with use eigen vectors but no standardization
				HOGGART::eigen_hoggard(parAA, parAB, parBB, parY, beta_hat, se_hat, p_hat, varAA, varAB, varBB, false);
				double hogg3 = beta_hat;
				double phogg3 = p_hat;
				// reference hoggart test from quicktest with standardization
				HOGGART::quicktest_hoggard(gAA, gAB, gBB, gY, double(gY.size()), beta_hat, se_hat, p_hat, varAA, varAB, varBB);
				double hogg4 = beta_hat;
				double phogg4 = p_hat;
				
			}
			std::cout << " | " << tmrHogg.elapsed()/1000;
			
			//Timer tmr_g;
			//for(int i=0; i<NTimes; i++){
			//	aic = EM::gauss_em(parSNP, parY, outParams,	false, 2, 2, 2);
			//}
			//std::cout << " | " << tmr_g.elapsed()/1000;

		}
	}

	return 0;
}


// valarray helpers x apply func.

double clear(double x)
{
    return 0.0;
}

double print(double x)
{
    std::cout << x;
    return x;
}

bool isNA(double i)
{
    return (i == NA);
}

bool isNA(std::vector<double> x)
{
    bool res = std::any_of(x.begin(), x.end(), [](double i) { return (i == NA); });
    return res;
}