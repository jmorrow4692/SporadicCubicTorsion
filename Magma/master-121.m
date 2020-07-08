//////////////////////////////////////////////////////////////////////////////////////////
//  X1(121) --> X0(121)
//////////////////////////////////////////////////////////////////////////////////////////
 
 N := 121; 

/****************************************************************************** 
 Here is a summary of the argument.

 X := X_0(121) is a genus 6 curve, not trigonal (via canonical model computation)
 J := J_0(121) is isogenous to a product of 6 elliptic curves (see below)

 X has 12 cusps, 2 rational (0 and infinity), and one Galois orbit of size 10, defined over Q(zeta_11)
 X_0(121)(F_5) has 2 rational cusps  (see Remark 7.2 why we use F_5 instead of F_3).

   Quick argument: 0 - infty is torsion, so is non-zero mod 5 since torsion injects;
   the prime 5 splits in Q(zeta_11) as 2 primes with inertia degree 5, so those cusps are only defined over F_5^5
   
 A cubic point of X_1(121)(Q) has bad reduction at 5, since there are no elliptic curves 
 over F_5^n (n = 1..3) with a rational 121 torsion point (by brute force search).
 Thus a cubic point reduces to a divisor that is a sum of cusps. 
 Because the inertia degree of 5 in Q(zeta_11) is 5, the cusps are a combination of 0 and infinity. 
 
 We can thus verify the formal immersion criterion at the degree 3 effective divisors supported at 
 infinity and 0 (using the Atkin-Lehner involution to compute the q-expansion at 0). 
******************************************************************************/

  //////////////////////////////////////////////////////////////////////////////////////////
  // Check J_0(121) decomposition and ranks
  //////////////////////////////////////////////////////////////////////////////////////////

  N := 121;
  dec := Decomposition(JZero(N));
  for i in [1..#dec] do
    Dimension(dec[i]),  IsZeroAt(LSeries(dec[i]),1);
  end for;  

  // Double check
  for e in dec do 
    Rank(EllipticCurve(e)); 
  end for;
  // 1,0,0,0,0,0       
 
/*  
    Modular abelian variety 65A of dimension 1, level 11^2 and conductor 11^2
    over Q,
    Modular abelian variety 65B of dimension 1, level 11^2 and conductor 11^2
    over Q,
    Modular abelian variety 65C of dimension 1, level 11^2 and conductor 11^2
    over Q,
    Modular abelian variety 65D of dimension 1, level 11^2 and conductor 11^2
    over Q,
    Modular abelian variety N(11,65,1)(11A) of dimension 1, level 11^2 and
    conductor 11 over Q,
    Modular abelian variety N(11,65,11)(11A) of dimension 1, level 11^2 and
    conductor 11 over Q  
*/ 
 


  //////////////////////////////////////////////////////////////////////////////////////////
  // Equations (via Ozman--Siksek's code) for the canonical model of X_0(121). Its not trigonal.
  //////////////////////////////////////////////////////////////////////////////////////////  

  F := Rationals();
  F := FiniteField(3);
  P<[x]> := PolynomialRing(F,6);
  eqns := [
  x[1]*x[3] - x[2]^2 - 41*x[2]*x[6] - x[3]^2 - 4*x[3]*x[5] - 12*x[3]*x[6] + 2*x[4]^2 - 8*x[4]*x[5] +
    32*x[4]*x[6] - 15*x[5]^2 + 5*x[5]*x[6] - 20*x[6]^2,
x[1]*x[4] - x[2]*x[3] - 24*x[2]*x[6] - x[3]^2 + x[3]*x[4] - 3*x[3]*x[5] - 8*x[3]*x[6] + x[4]^2 - 4*x[4]*x[5] +
    20*x[4]*x[6] - 9*x[5]^2 + 2*x[5]*x[6] - 11*x[6]^2,
x[1]*x[5] - 20*x[2]*x[6] - x[3]*x[4] - 2*x[3]*x[5] - 8*x[3]*x[6] + 2*x[4]^2 - 3*x[4]*x[5] + 18*x[4]*x[6] -
    8*x[5]^2 + 3*x[5]*x[6] - 2*x[6]^2,
x[1]*x[6] - 2*x[2]*x[6] - x[3]*x[6] + 2*x[4]*x[6] - x[5]^2,
x[2]*x[4] - 37*x[2]*x[6] - x[3]^2 - x[3]*x[4] - 5*x[3]*x[5] - 11*x[3]*x[6] + 3*x[4]^2 - 6*x[4]*x[5] +
    29*x[4]*x[6] - 14*x[5]^2 + 4*x[5]*x[6] - 8*x[6]^2,
x[2]*x[5] + 34*x[2]*x[6] + 3*x[3]*x[5] + 14*x[3]*x[6] - x[4]^2 + 5*x[4]*x[5] - 31*x[4]*x[6] + 13*x[5]^2 -
    3*x[5]*x[6] + 12*x[6]^2];
  Xcan := Curve(ProjectiveSpace(P), eqns);
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Verify the formal immersion criterion
  //////////////////////////////////////////////////////////////////////////////////////////   
  
  M := CuspForms(N);
  R<q> := PowerSeriesRing(Integers());

    dec := Decomposition(JZero(N));
    for i in [2..5] do 
	Rank(EllipticCurve(dec[i])); 
    end for;

  for tup in [[3-j, j] : j in [0..3]] do
      mat := Matrix([          
          [Coefficient(f,k) : k in [1..tup[1]]] cat 
          [Coefficient(AtkinLehnerOperator(121,f),k) : k in [1..tup[2]]]
          where f is M!(R!Newform(dec[i])) : 
          i in [2..5]
      ]);
      mat;
      Rank(RMatrixSpace(GF(5),4,3)!mat);
  end for;
      
      
