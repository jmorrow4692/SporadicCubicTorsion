////////////////////////////////////////////////////////////////////////////////////////
// This is a complete determination of the cubic points on X1(45)
//////////////////////////////////////////////////////////////////////////////////////////

 N := 45;


/*****************************************************************************************
 Here is a summary of the argument.

 There is a quotient XH of genus 5, where H is the index 2 subgroup of B(45) with upper
 left entry congruent to 1,4,11, or 14 mod 15. 

 See the Maple file FindModel_XH_45_input for a computation of a model, and the paper
 for an explanation of the computation.

 There are 4 rational points and 6 quadratic points (all cuspidal), and
 there are actually 8 cubic points on XH. This gives 52 cubic divisors.

 The cubic points are not all cusps, so we check that the cubic points do not lift to
 X_1(45). 

 The cusps generate a subgroup of J_H isomorphic to  [ 2, 4, 48, 0 ]
 The upper bound we get from local methods is  [ 2, 4, 24, 48, 0 ]
 The Hecke bound gives [2,4,48]

*****************************************************************************************/


  //////////////////////////////////////////////////////////////////////////////////////////
  // Input the homebrewed functions
  //////////////////////////////////////////////////////////////////////////////////////////

  load "functions.m";

  F := Rationals();

  //////////////////////////////////////////////////////////////////////////////////////////
  // Precomputed equations for X_H(45) see file FindModel_XH_45_input
  //////////////////////////////////////////////////////////////////////////////////////////
  
  A2<x,y> := AffineSpace(F,2);
  
  fH45 := x^4*y^2+3*x^3*y^3+x^2*y^4+2*x^4*y+13*x^3*y^2+13*x^2*y^3+2*x*y^4+x^4+17*x^3*y+39*x^2*y^2+17*x*y^3+y^4+6*x^3+39*x^2*y+39*x*y^2+6*y^3+9*x^2+27*x*y+9*y^2;
  XH := ProjectiveClosure(Curve(A2,fH45));

  
  //////////////////////////////////////////////////////////////////////
  // Compute a bound on the torsion.
  //////////////////////////////////////////////////////////////////////

  torsionData := {@@};
//  for p in [q : q in PrimesUpTo(200) | not q in PrimeDivisors(2*N) cat [11]] do
  for p in [7,17] do      
      invs := Invariants(ClassGroup(Curve(Reduction(XH,p))));
      torsionData := Include(torsionData, invs);
      <p,torsBound(torsionData)>;
  end for;

  "The rational torsion subgroup is a subgroup of", torsBound(torsionData); // [ 48, 8, 4, 2 ]  


  //////////////////////////////////////////////////////////////////////////////////////////
  // Construct known small degree points.
  //////////////////////////////////////////////////////////////////////////////////////////

  pts := PointSearch(XH,1000);  // 5 pts, one is singular, all cusps
    pts := [D : D in &cat([Places(D) : D in pts ]) | Degree(D) eq 1];
  eqns2 := divisorSearch(XH,2,5); // 6, all cusps
  eqns3 := divisorSearch(XH,3,3); // 6, sporadic, need 2 more

  // find the other degree 3 points via automorhpisms
time  A := AutomorphismGroup(XH);
  divsH3   := {@XH meet Scheme(XH,I) : I in eqns3 @};
    divsH3 := {@Divisor(XH, XH meet D) : D in divsH3@}; 
    divsH3 := divsH3 join {@Pullback(a,D) : D in divsH3, a in {A.1} @};
    assert #divsH3 eq 8; // found them

  divsH2   := {@XH meet Scheme(XH,I) : I in eqns2 @};

  assert {Degree(D) : D in divsH2} eq {2};
  assert {IsIrreducible(D) : D in divsH2} eq {true};

  
  //////////////////////////////////////////////////////////////////////////////////////////
  // See what subgroup the known divisors generate  
  // (Easier mod a prime, and the quadratic points suffice)
  //////////////////////////////////////////////////////////////////////////////////////////  

  p := 7;
  XHp<[X]> := Curve(Reduction(XH,p));
  Rp := Parent(X[1]);
  basePt  := Divisor(XHp![-3,0,1]);
    pic, mPic := ClassGroup(XHp);    
    idealsHp := {@[Rp!e : e in eqns2[i] ] : i in [2,4,5]@};
    divsHp   := {@Divisor(XHp,Scheme(XHp,I)) : I in idealsHp@};
    
  global, mGlobal := sub<pic | [(Inverse(mPic))(D - Degree(D)*basePt) : D in divsHp]>;
     Invariants(global);    
  
    //-------------------------------------------------------------------------------------
    // returns [ 48, 4, 2 ]
    //-------------------------------------------------------------------------------------

  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Hard code the generators for torsion.
  //////////////////////////////////////////////////////////////////////////////////////////  
   
  ideals := 
  {@
      [
          X[2]^2 - X[2]*X[3] - 3*X[3]^2,
          X[1] + X[2]
      ],
      [
          X[1]^2 + 3*X[1]*X[2] + X[2]^2,
          X[3]
      ],
      [
          X[2]^2 + 3*X[3]^2,
          X[1] + 1/2*X[2] + 3/2*X[3]
      ]
  @};

  idealsHp := {@[Rp!e : e in I ] : I in idealsHp@};
  divsHp   := {@Divisor(XHp,Scheme(XHp,I)) : I in idealsHp@};
      
  global, mGlobal := sub<pic | [(Inverse(mPic))(D - Degree(D)*basePt) : D in divsHp]>;
     Invariants(global); // double check   


  //////////////////////////////////////////////////////////////////////
  // Compute the image of Abel--Jacobi mod 7
  //////////////////////////////////////////////////////////////////////
  
  "There are", [#Places(XHp,i) : i in [1..3]], "places of degree 1, 2, and 3 over F_7"; // [8, 28, 112];  

  pl1 := Places(XHp,1);
  cubicDivisors :=
           {@ &+[Divisor(pl1[i])*tup[i] : i in [1..#pl1]] : tup in CartesianPower([0..3], #pl1) | &+[i : i in tup] eq 3@}
      join {@Divisor(D1) + Divisor(D2) : D1 in Places(XHp,1), D2 in Places(XHp,2)@}
      join {@Divisor(D) : D in Places(XHp,3)@};

  validCubicImages := {@@};
  for D0 in cubicDivisors do
      D := D0 - Degree(D0)*basePt;
      if Inverse(mPic)(D) in global then
        validCubicImages :=
        validCubicImages join {@Inverse(mPic)(D)@};
      end if;
  end for;
  #validCubicImages, "map to reductions of torsion divisors"; // 52
  "There are", Binomial(#pts + 3 - 1, 3) + #pts * #eqns2 + #divsH3, "rational cubic divisors, so we have found all of the points!";

  
  ////////////////////////////////////////////////////////////////////////
  // Verify that the cubic points do not lift to cubic points of X1(45)
  ////////////////////////////////////////////////////////////////////////       

  for i in [2..#divsH3] do
    D := divsH3[i];
    pt := RepresentativePoint(Support(D)[1]);
    K := Parent(pt[1]);
    Q<Z> := PolynomialRing(K);

    // Computed in the maple file Lift_XY_back_to_X1_45_input; given a cubic point with coordinates X,Y on XH(45),
    // the roots of this are the x coordinates of the preimage in X1_45

    f := Z^6 - (X^4+3*X^3*Y-2*X^2*Y^2-X*Y^3+7*X^3+8*X^2*Y-5*X*Y^2+15*X^2+9*X*Y+9*X+12*Y)/(X^2*Y+X^2+2*X*Y+X+2*Y)*Z^5
     + (X^5*Y+3*X^4*Y^2+X^3*Y^3+2*X^5+17*X^4*Y+24*X^3*Y^2+6*X^2*Y^3+23*X^4+94*X^3*Y+68*X^2*Y^2+9*X*Y^3
       +100*X^3+216*X^2*Y+75*X*Y^2+4*Y^3+192*X^2+201*X*Y+24*Y^2+135*X+51*Y)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z^4
     - (4*X^7*Y+12*X^6*Y^2+2*X^5*Y^3+8*X^7+74*X^6*Y+110*X^5*Y^2+17*X^4*Y^3+95*X^6+454*X^5*Y+418*X^4*Y^2+51*X^3*Y^3
       +460*X^5+1368*X^4*Y+836*X^3*Y^2+71*X^2*Y^3+1163*X^4+2267*X^3*Y+912*X^2*Y^2+42*X*Y^3+1622*X^3+2051*X^2*Y
       +492*X*Y^2+5*Y^3+1176*X^2+879*X*Y+90*Y^2+336*X+105*Y)/(X+1)/(X*Y+3*X+3*Y+3)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z^3
     - (X^5*Y+3*X^4*Y^2+X^3*Y^3+9*X^4*Y+18*X^3*Y^2+4*X^2*Y^3+X^4+36*X^3*Y+43*X^2*Y^2+6*X*Y^3+9*X^3+71*X^2*Y
       +47*X*Y^2+3*Y^3+27*X^2+63*X*Y+18*Y^2+27*X+12*Y)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z^2
     + (X^5*Y+3*X^4*Y^2+X^3*Y^3+X^5+12*X^4*Y+18*X^3*Y^2+4*X^2*Y^3+10*X^4+50*X^3*Y+41*X^2*Y^2+5*X*Y^3+37*X^3
       +91*X^2*Y+39*X*Y^2+2*Y^3+60*X^2+69*X*Y+12*Y^2+36*X+12*Y)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z
     + 1
	where X is pt[1] where Y is pt[2];
   assert IsIrreducible(f);
  end for;

  // We need to handle i = 1 separately; when we evaluate the coefficient of Z^3 at divsH3[1] we get 0/0
  // If we evaluate the Y coordinate first, then the X coordinate, then this resolves the issue
  
  ppt := RepresentativePoint(Support(divsH3[1])[1]);
  K := Parent(ppt[1]);
  Q<X,Y> := PolynomialRing(K,2);
  _<Z> := PolynomialRing(FieldOfFractions(Q));
    f := Z^6 - (X^4+3*X^3*Y-2*X^2*Y^2-X*Y^3+7*X^3+8*X^2*Y-5*X*Y^2+15*X^2+9*X*Y+9*X+12*Y)/(X^2*Y+X^2+2*X*Y+X+2*Y)*Z^5
    + (X^5*Y+3*X^4*Y^2+X^3*Y^3+2*X^5+17*X^4*Y+24*X^3*Y^2+6*X^2*Y^3+23*X^4+94*X^3*Y+68*X^2*Y^2+9*X*Y^3
       +100*X^3+216*X^2*Y+75*X*Y^2+4*Y^3+192*X^2+201*X*Y+24*Y^2+135*X+51*Y)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z^4
    - (4*X^7*Y+12*X^6*Y^2+2*X^5*Y^3+8*X^7+74*X^6*Y+110*X^5*Y^2+17*X^4*Y^3+95*X^6+454*X^5*Y+418*X^4*Y^2+51*X^3*Y^3
       +460*X^5+1368*X^4*Y+836*X^3*Y^2+71*X^2*Y^3+1163*X^4+2267*X^3*Y+912*X^2*Y^2+42*X*Y^3+1622*X^3+2051*X^2*Y
       +492*X*Y^2+5*Y^3+1176*X^2+879*X*Y+90*Y^2+336*X+105*Y)/(X+1)/(X*Y+3*X+3*Y+3)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z^3
    - (X^5*Y+3*X^4*Y^2+X^3*Y^3+9*X^4*Y+18*X^3*Y^2+4*X^2*Y^3+X^4+36*X^3*Y+43*X^2*Y^2+6*X*Y^3+9*X^3+71*X^2*Y
       +47*X*Y^2+3*Y^3+27*X^2+63*X*Y+18*Y^2+27*X+12*Y)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z^2
    + (X^5*Y+3*X^4*Y^2+X^3*Y^3+X^5+12*X^4*Y+18*X^3*Y^2+4*X^2*Y^3+10*X^4+50*X^3*Y+41*X^2*Y^2+5*X*Y^3+37*X^3
       +91*X^2*Y+39*X*Y^2+2*Y^3+60*X^2+69*X*Y+12*Y^2+36*X+12*Y)/(X^2*Y+X^2+2*X*Y+2*X+Y)*Z
    + 1;
  coeffsFirst  := [Evaluate(c,[X,ppt[2]]) : c in Coefficients(f)];  
  coeffsSecond := [Evaluate(c,[ppt[1],0]) : c in coeffsFirst];  
  _<T> := PolynomialRing(K);
  f := &+[T^(i-1)*coeffsSecond[i] : i in [1..#coeffsSecond]];
  assert IsIrreducible(f);
