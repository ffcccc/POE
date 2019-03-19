## Welcome to the POE Page

The POE software implements statistical methods for detection of Parent of Origin Effect (POE), in dataset of unrelated individuals. The aim is to suggest hints for model selection in genotype association testing.

Release v0.1 includes different implementations of the method described by C. Hoggarth published in PLOS Genetics (Novel Approach Identifies SNPs in SLC2A10 and KCNK9 with Evidence for Parent-of-Origin Effect on Body Mass Index). 

The source code and example files can be browsed and downloaded from [github](https://github.com/ffcccc/POE).

The quicktest_hoggard() and eigen_hoggard() algorithms in large part follow the flux originally conceived by Toby Johnson and Zoltán Kutalik in the Quicktest sw: the first implementation is a C++ porting to simplify the interface to STL, the second is hacked to overcome some memory bottlenecks and exploit a little bit of the power of the Eigen framework.
Finally, the code of nm_hoggart() follows the approach of C. Sacco published in Statistical Applications in Genetics and Molecular Biology (A statistical test for detecting parent-of-origin effects when parental information is missing).
### Usage

POE is a lightweight C++ library. Where feasible It makes use of vectorization to enhance performances, referring to the header-only Eigen::Array library and to the helper library EigenUtils.

### Usage
See the test.cpp file for a working example.

```markdown
    ...
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
    ...
```

### References
-[Novel Approach Identifies SNPs in SLC2A10 and KCNK9 with Evidence for Parent-of-Origin Effect on Body Mass Index](https://doi.org/10.1371/journal.pgen.1004508)

-[Quicktest](https://www2.unil.ch/cbg/index.php?title=QuickTest)

-[A statistical test for detecting parent-of-origin effects when parental information is missing](https://doi.org/10.1515/sagmb-2017-0007).

