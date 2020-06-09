Here is code verifying the computations and theorems of the manuscript 
                                          
                                          "Sporadic Cubic Torsion" 
                                          
             By Maarten Derickx, Anastassia Etropolski, Mark van Hoeij, Jackson S. Morrow, and David Zureick-Brown

In order to utilize this code, one needs to do the following first:

0. Upload the Magma files (The two latter .m files come Ozman--Siksek "Quadratic points on modular curves")
  - functions.m
  - modEqns.m 
  - quadPts.m 

1. Upload the following Python files from the folder mdsage (These files come from Maarten Derickx's Github https://koffie.github.io/mdsage/doc/html/index.html)
  - cuspidal_classgroup.py
  - maartens_sage_functions.py
  - modular_unit_divisors.py

Once the above files are loaded, you will be able to run the below code verifying claims found in our paper.

The code verifying the ranks of modular Jacobians computations from Section 3 can be found at master-ranks.m

The code verifying the cuspidal torsion computations from Section 4 can be found at torsionComputations.py

For each specific case, code verifiying the claims in the paper can be found at master-N.m (E.g., master-22.m corresponds to code verifying the case of X1(22))

There are two additional files 
  - master-45-precomputations.m (Computes a model for the intermediate curve XH(45) from Section 6.7)
  - master-45-verifications.m   (Verifies that our model for XH(45) is isomorphic to a model constructed from David Zywina's code from "Computing actions on cusp forms" http://pi.math.cornell.edu/~zywina/papers/AtkinLehner/ComputeAtkinLehner.txt)


