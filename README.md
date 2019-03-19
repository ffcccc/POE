## Welcome to the POE Page

The POE software implements statistical methods for detection of Parent of Origin Effect (POE), in dataset of unrelated individuals. The aim is to suggest hints for model selection in genotype association testing.

Release v0.1 includes different implementations of the method described by C. Hoggarth in <ref>. 

The source code and example files can be browsed and downloaded from the github site [editor on GitHub](https://github.com/ffcccc/POE).

This work in large part follows the code was originally written by Toby Johnson and Zoltán Kutalik in <ref>.

### Usage

POE is a lightweight C++ library. Where feasible It makes use of vectorization to enhance performances, referring to the header-only Eigen::Array library and to the helper library EigenUtils.

### Usage
See the test.cpp file for a working example.

```markdown
				// reference hoggart test from quicktest with standardization
				HOGGART::quicktest_hoggard(gAA, gAB, gBB, gY, double(gY.size()), beta_hat, se_hat, p_hat, varAA, varAB, varBB);
				double hogg4 = beta_hat;
				double phogg4 = p_hat;ck
```

### References
-For more details see [Quicktest](https://www2.unil.ch/cbg/index.php?title=QuickTest).
-For more details see [Paper](https://doi.org/10.1371/journal.pgen.1004508).
