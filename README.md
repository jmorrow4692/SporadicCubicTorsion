Here is code verifying the computations and theorems of the manuscript 

                                          "Sporadic Cubic Torsion" 
                                          
    By Maarten Derickx, Anastassia Etropolski, Mark van Hoeij, Jackson S. Morrow, and David Zureick-Brown

There are 3 folders in this repository. One corresponding to each computer algebra program we use: Magma, Sage, and Maple.

In order to utilize this code, one needs to do the following:

0. Upload the Magma files (The two latter .m files come Ozman--Siksek "Quadratic points on modular curves")
  - functions.m
  - modEqns.m 
  - quadPts.m 

1. Upload the following Python files from the folder mdsage (These files come from Maarten Derickx's Github https://koffie.github.io/mdsage/doc/html/index.html)
  - cuspidal_classgroup.py
  - maartens_sage_functions.py
  - modular_unit_divisors.py
(Alternatively, one could download the .zip file mdsage.zip and unzip this in your local directory).

Once the above files are loaded, you will be able to run the below code verifying claims found in our paper.

In the Magma folder, one finds:
  - master-ranks.m (Magma code verifying the ranks of modular Jacobians computations from Section 3)
  - master-N.m (Magma code for the computations in Sections 6 and 7 of the paper can be found at. For example, the file master-22.m corresponds to code verifying       the case of X1(22)).

In the Sage folder, one finds:
  - torsionComputations.py (Sage code verifying the cuspidal torsion computations from Section 4).
  
In the Maple folder, one finds:
  - Cuspidal_part_of_J1N_output (Output of algorithm CuspidalClassGroupStructure from http://www.math.fsu.edu/~hoeij/files/X1N/cusp_divisors_program which gives the     structure of the subgroup of J1(N) generated by the Galois orbits of cusps for N <= 100.)
  - FindModel_XH_45_output (Gives a singular model for the modular curve XH(45) from Subsection 6.7)
  - Lift_XY_back_to_X1_45_output (Gives a way to lift a cubic point on XH(45) back to X1(45) to check whether it gives a cubic point on X1(45))
  - Additionally, the input files of these.


