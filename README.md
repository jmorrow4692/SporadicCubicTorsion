Here is code verifying the computations and theorems of the manuscript 

                                          "Sporadic Cubic Torsion" 
                                          
    By Maarten Derickx, Anastassia Etropolski, Mark van Hoeij, Jackson S. Morrow, and David Zureick-Brown

There are 3 folders in this repository. One corresponding to each computer algebra program we use: Magma, Sage, and Maple.

In order to utilize this code, one needs to do the following:

0. Download the following Magma files (The two latter .m files come Ozman--Siksek "Quadratic points on modular curves" https://arxiv.org/abs/1806.08192)
    - functions.m
    - modEqns.m 
    - quadPts.m 

1. Download the following Python files (These files come from Maarten Derickx's Github https://koffie.github.io/mdsage/doc/html/index.html)
   - cuspidal_classgroup.py
   - maartens_sage_functions.py
   - modular_unit_divisors.py

Once the above files are loaded, you will be able to run the below code verifying claims found in our paper.

In the Magma folder, one finds:
   - master-ranks.m (Magma code verifying the ranks of modular Jacobians computations from Section 3)
   - master-N.m (Magma code for the computations in Sections 6 and 7 of the paper can be found at. For example, the file master-22.m corresponds to code verifying       the case of X1(22)).
   - ModelForXH45.txt (A text file with the model for the modular curve XH(45) from Subsection 6.7 and a way to lift a cubic point on XH(45) back to X1(45))

In the Sage folder, one finds:
  - torsionComputations.py (Sage code verifying the cuspidal torsion computations from Section 4).
  
In the Maple folder, one finds:
  - Cuspidal_part_of_J1N_output (Output of algorithm CuspidalClassGroupStructure from   
    http://www.math.fsu.edu/~hoeij/files/X1N/cusp_divisors_program which gives the structure of the subgroup of J1(N)(Q) generated by the Galois orbits of cusps)
  - Cuspidal_part_of_J2_2N (Output of a Maple implementation (not yet ready for distribution) which gives the structure of the subgroup of J1(2,2N)(Q) generated by 
    the Galois orbits of cusps. These results agree with the output from the Sage program in torsionComputations.py)
  - FindModel_XH_45_output (Gives a singular model for the modular curve XH(45) from Subsection 6.7)
  - Lift_XY_back_to_X1_45_output (Gives a way to lift a cubic point on XH(45) back to X1(45) to check whether it gives a cubic point on X1(45))
  - Model_X1_54_index_2 (Determines the index 2 subfield of Q(X1(54)), lists generators M1, M2 of this subfield, gives their algebraic relation, and gives a model 
    for the index 2, genus 22 quotient of X1(54))
  - Additionally, the input files of these.



