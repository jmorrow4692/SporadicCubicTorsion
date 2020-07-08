//////////////////////////////////////////////////////////////////////////////////////////
//  X1(65) --> X0(65)
//////////////////////////////////////////////////////////////////////////////////////////

load "modEqns.m"; 
 N := 65; 

/****************************************************************************** 
 Here is a summary of the argument.

 X := X_0(65) is a genus 5 curve, not trigonal (via canonical model)
 J := J_0(65) is isogenous to E x A1 x A2, 

 X has 4 cusps (0, infinity, 1/5, 1/13). 
 X_0(65)(F_3) has 4 points. Two are rational cusps. 

A cubic point of X_1(65)(Q) has bad reduction at 3, since there are no elliptic curves over 
F_3^n (n = 1..3) with a rational 65 torsion point by the Hasse bound trick. Thus a cubic point 
reduces to a divisor that is a sum of cusps. 

The map X_0(65)^{(3)}  \to J_e(65) to the winding quotient is not a formal immersion at all 
cuspidal divisors. Fortunately, by Lemma 7.4, we only need to verify the formal immersion 
criterion at particular sums of cusps on $X_0(65)$: those of the form 3*[infinity], 3*[0], 
3*[1/5], and 3*[1/13].

******************************************************************************/

  //////////////////////////////////////////////////////////////////////////////////////////
  // Check J_0(65) decomposition and ranks
  //////////////////////////////////////////////////////////////////////////////////////////

  N := 65;
  dec := Decomposition(JZero(N));
  for i in [1..#dec] do
    Dimension(dec[i]),  IsZeroAt(LSeries(dec[i]),1);
  end for;  
  
  Rank(EllipticCurve(dec[1]));  // 1

/*  
    Modular abelian variety 65A of dimension 1, level 5*13 and conductor 5*13
    over Q,
    Modular abelian variety 65B of dimension 2, level 5*13 and conductor
    5^2*13^2 over Q,
    Modular abelian variety 65C of dimension 2, level 5*13 and conductor
    5^2*13^2 over Q

*/ 
     

  //////////////////////////////////////////////////////////////////////////////////////////
  // Equations (from Samir Siksek) for the canonical model of X_0(65). Its not trigonal.
  //////////////////////////////////////////////////////////////////////////////////////////  

  X := modeqns(N,13); // 3 quadratics in P^4, so not trigonal
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Verify the formal immersion criterion
  //////////////////////////////////////////////////////////////////////////////////////////   
  
  M := CuspForms(N);
    B := Basis(M);
  R<q> := PowerSeriesRing(Integers());

  
  // Find simultaneous eigenspaces
  // Find a p such that T_p acts with rational eigenvalues
  for p in PrimesUpTo(200) do 
     <p, Coefficient(Newforms(M)[1][1],p), 
         Coefficient(Newforms(M)[2][1],p), 
         Coefficient(Newforms(M)[3][1],p)>; 
  end for; 

  // p = 199 works
assert (IsIntegral(Coefficient(Newforms(M)[1][1],199)) eq true);
assert (IsIntegral(Coefficient(Newforms(M)[2][1],199)) eq true);
assert (IsIntegral(Coefficient(Newforms(M)[3][1],199)) eq true);
  
  // Diagonalize
  A := HeckeOperator(M,199);
  Eigenvalues(A); // { <20, 2>, <-16, 1>, <4, 2> }
  E1 := Basis(Eigenspace(A,-16));
  E2 := Basis(Eigenspace(A,20));
  E3 := Basis(Eigenspace(A,4)); 
  
  betterBasis := 
  [
      &+[E1[1][i]*B[i] : i in [1..#B]],
      &+[E2[1][i]*B[i] : i in [1..#B]],
      &+[E2[2][i]*B[i] : i in [1..#B]],
      &+[E3[1][i]*B[i] : i in [1..#B]],
      &+[E3[2][i]*B[i] : i in [1..#B]]
  ];
  
  
  // E2 and E3 correspond to rank 0 factors. Use these in the formal immersion criterion
  tups := [tup : tup in CartesianPower([0..3],4) | &+[tup[i] : i in [1..4] ] eq 3];

  for tup in tups do
      "************************************";
      mat := Matrix([          
          [Coefficient(f,k) : k in [1..tup[1]]] cat 
          [Coefficient(AtkinLehnerOperator(5, f),k) : k in [1..tup[2]]] cat 
          [Coefficient(AtkinLehnerOperator(13,f),k) : k in [1..tup[3]]] cat 
          [Coefficient(AtkinLehnerOperator(65,f),k) : k in [1..tup[4]]]
          where f is betterBasis[i]: 
          i in [2..5]
      ]);
      Rank(RMatrixSpace(GF(3),4,3)!mat), tup, mat;
//      Rank(RMatrixSpace(Rationals(),4,3)!mat), tup, mat;      
  end for;

  // By Lemma 7.4, we are done as the matrices above have rank 3 for the quadruples <3,0,0,0>, <0,3,0,0>, <0,0,3,0>, and <0,0,0,3> 
      
  // Reality check: the map Sym^3 X_0(65)--> J_0(65) is an immersion..
  tups := [tup : tup in CartesianPower([0..3],4) | &+[tup[i] : i in [1..4] ] eq 3];

  for tup in tups do
      mat := Matrix([          
          [Coefficient(f,k) : k in [1..tup[1]]] cat 
          [Coefficient(AtkinLehnerOperator(5, f),k) : k in [1..tup[2]]] cat 
          [Coefficient(AtkinLehnerOperator(13,f),k) : k in [1..tup[3]]] cat 
          [Coefficient(AtkinLehnerOperator(65,f),k) : k in [1..tup[4]]]
          where f is betterBasis[i]: 
          i in [1..5]
      ]);
      Rank(RMatrixSpace(GF(3),5,3)!mat), tup, mat;
//      Rank(RMatrixSpace(Rationals(),5,3)!mat);
  end for;

  // We get that all of the matrices have rank 3 as the map Sym^3 X_0(65)--> J_0(65) is a formal immersion because X_0(65) is not trigonal
      
