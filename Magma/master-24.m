/////////////////////////////////////////////////////////////////////////////////////////
// This is a complete determination of the cubic points on X1(24)
//////////////////////////////////////////////////////////////////////////////////////////
 
/****************************************************************************** 
 Here is a summary of the argument.
  - X1(24) has genus 5
  - The local torsion bound is [ 2, 2, 6, 120, 0 ]
  - The Hecke bound is [2, 2, 2, 120]
  - The additional argument from Theorem 4.13 proves that the torsion is [2,2,120]
  - There are 6 rational points and 5 quadratic points, which gives
    6*5 + Binomial(8,3) = 86 cubic divisors
  - Computing (mod 5) the intersection of the image of Abel--Jacobi 
  - and the reduction of the rational torsion gives 86 cubic divisors.
******************************************************************************/

  N := 24;

  //////////////////////////////////////////////////////////////////////////////////////////
  // Input the homebrewed functions 
  //////////////////////////////////////////////////////////////////////////////////////////
  
  load "functions.m";
  
  F := Rationals();     
  A2<x,y> := AffineSpace(F,2);
  
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Equations for X1(24), from Derickx--Sutherland
  //////////////////////////////////////////////////////////////////////////////////////////
      
    A2<x,y> := AffineSpace(F,2);
    X:=Curve(A2,(x^2 - 3)*y^4 + x*(x^4 + 2*x^3 + 2*x^2 + 2*x - 3)*y^2 - x^4 - x^2);
    Xp := ProjectiveClosure(X);           

  
  //////////////////////////////////////////////////////////////////////
  // Get the canonical model.
  // This presents it as a complete intersection of 3 quadratics in P4.
  //////////////////////////////////////////////////////////////////////
  
    phi := CanonicalMap(Xp);
    Xsm := CanonicalImage(Domain(phi),phi); 
     P<[T]> := AmbientSpace(Xsm);
    phi := map<Xp -> Xsm | DefiningEquations(phi)>;
  
  //////////////////////////////////////////////////////////////////////
  // Hard coded
  //////////////////////////////////////////////////////////////////////

      P<[T]> := ProjectiveSpace(F,4);    
      Xsm := Curve(P,
      [
      T[1]*T[2] + 6*T[2]^2 + 11*T[2]*T[3] + 4*T[3]^2 + T[4]^2 - 4*T[2]*T[5] - 15*T[3]*T[5],
      -T[2]^2 + T[1]*T[3] - 2*T[2]*T[3] - 6*T[3]^2 + T[2]*T[5] + 4*T[3]*T[5],
      T[2]^2 + 4*T[2]*T[3] + T[3]^2 + T[1]*T[5] + 2*T[2]*T[5] - 5*T[3]*T[5]    
      ]);

  //////////////////////////////////////////////////////////////////////
  // Compute the local torsion bound
  //////////////////////////////////////////////////////////////////////   

  torsData := {@@};      
//  for p in [q : q in PrimesUpTo(40) | not q in PrimeDivisors(2*N) ] do
  for p in [5 ,7] do     
      invs := Invariants(ClassGroup(Curve(Reduction(X,p))));
      torsData := torsData join {@invs@};
      <p,invs>;
  end for;
   
    "The rational torsion subgroup is a subgroup of", torsBound(torsData); ; // [ 120, 6, 2, 2 ]
    "The Hecke bound is",  [ 120, 2, 2, 2 ]; // [ 120, 2, 2, 2 ], computed in sageComputations
    
    "The additional argument from Theorem 4.13 gives a bound of", [120,2,2]; // [ 120, 2, 2, 2 ], computed in sageComputations

    
  ////////////////////////////////////////////////////////////////////////
  // Compute the known small degree points
  // (Hard coded)
  ////////////////////////////////////////////////////////////////////////

  C := Xsm;
  basePt := [ 1 , 0 , 0 ,  0 , 0];

  
  // Degree 1  
  pts := PointSearch(C,20);
  pts := {@ [ 0 , 0 , 0 ,  0 , 1], 
            [-3 , 1 , 0 , -1 , 1],
            [-3 , 1 , 0 ,  1 , 1], 
            [-14, 4, -1,  -1, 1],
            [-14, 4, -1,   1, 1],
            [ 1 , 0 , 0 ,  0 , 0] @};
 
  divs1 := {@Divisor(Xsm!pts[i]) - Divisor(Xsm!basePt) : i in [1..#pts]@};

  
  // Degree 2
  // pts2 := divisorSearch(Xsm,2,3); // takes 88 seconds
 
  P4<[Y]> := AmbientSpace(C);
  eqns2 := [
      [
          Y[2]^2 - 4*Y[2]*Y[5] + 17/4*Y[5]^2,
          Y[1] + 4*Y[2] - 3/2*Y[5],
          Y[3] + 1/2*Y[5],
          Y[4]
      ],
      [
          Y[3]^2 + Y[3]*Y[5] + 1/3*Y[5]^2,
          Y[1] - 6*Y[3] + 4*Y[5],
          Y[2] + 2*Y[3] - Y[5],
          Y[4]
      ],
      [
          Y[2]^2 + 4*Y[2]*Y[3] + Y[3]^2,
          Y[1] + 2*Y[2] - 5*Y[3],
          Y[4],
          Y[5]
      ],
      [
          Y[3]^2 + Y[3]*Y[5] + 1/3*Y[5]^2,
          Y[1] - 15*Y[3] - 3*Y[5],
          Y[2] + 5*Y[3] + Y[5],
          Y[4] - Y[5]
      ],
      [
          Y[3]^2 + Y[3]*Y[5] + 1/3*Y[5]^2,
          Y[1] - 15*Y[3] - 3*Y[5],
          Y[2] + 5*Y[3] + Y[5],
          Y[4] + Y[5]
      ]  
          ];      
  
  divs2 := {@Divisor(C, C meet Scheme(C,e)) : e in eqns2@};
  "The quadratic cusps are defined over fields with discriminant", 
  [ Discriminant(Integers(Parent(RepresentativePoint(Support(pt)[1])[1]))) : pt in divs2];  
  divs2 := {@D - Degree(D)*Divisor(Xsm!basePt) : D in divs2@};
 
  // degree 4
//  pts4 := divisorSearch(Xsm,4,2);   // 21 seconds
  eqns4 := [
      [
	  T[3]^2 + T[3]*T[5] + 1/3*T[5]^2,
	  T[4]^2,
	  T[1] - 6*T[3] + 4*T[5],
	  T[2] + 2*T[3] - T[5]
      ],
      [
	  T[3]^2 + T[3]*T[5] + 1/2*T[5]^2,
	  T[4]^2 + 2*T[3]*T[5] + T[5]^2,
	  T[1] - 9*T[3] + T[5],
	  T[2] + 3*T[3]
      ],
      [
	  T[2]^2 + 4*T[2]*T[3] + T[3]^2,
	  T[4]^2,
	  T[1] + 2*T[2] - 5*T[3],
	  T[5]
      ],
      [
	  T[3]^2 + T[3]*T[5] + 1/6*T[5]^2,
	  T[4]^2 + 6*T[3]*T[5] + T[5]^2,
	  T[1] - 21*T[3] - T[5],
	  T[2] + 5*T[3]
      ]
  ];

  divs4 := {@Divisor(C, C meet Scheme(C,e)) : e in eqns4@};
  divs4 := {@D - Degree(D)*Divisor(Xsm!basePt) : D in divs4@};
    
  divs := divs1 join divs2 join divs4;

    
  //////////////////////////////////////////////////////////////////////
  // See what these divisors generate
  //////////////////////////////////////////////////////////////////////   
 
  p := 5;
  Cp<[T]> := Curve(Reduction(Xsm,p));
    pic,mPic := ClassGroup(Cp);    
  basePtp := Divisor(Cp!basePt);

  divs1p := {@Divisor(Cp!pts[i]) - basePtp : i in [1..#pts]@};
  divs2p := {@Divisor(Cp, Cp meet Scheme(Cp,[Parent(T[1])!e :e in I])) : I in eqns2@};
    extraTorsion := Support(divs2p[1])[1] - Support(divs2p[1])[2];
    divs2p := {@D - Degree(D)*basePtp : D in divs2p@};
  divs4p := {@Divisor(Cp, Cp meet Scheme(Cp,[Parent(T[1])!e : e in I])) : I in eqns4@};
    divs4p := {@D - Degree(D)*basePtp : D in divs4p@};

  divsp := divs1p join divs2p join divs4p;

  global, mGlobal := sub<pic | [(Inverse(mPic))(D) : D in divsp]>;
    "The known Q-rational divisors generate a subgroup isomorphic to", Invariants(global); // [ 2, 2, 120 ]
    
  // Did this part this part by hand
  basis := {@
            divsp[1] - divsp[8],
            divsp[3] + divsp[1] - divsp[8] - (3*(divsp[5] - divsp[7]) - divsp[4]),
            3*(divsp[5] - divsp[7]) - divsp[4]
             @};
      assert [Order((Inverse(mPic))(D)) : D in basis] eq  [ 2, 2, 120 ];  
      global, mGlobal := sub<pic | [(Inverse(mPic))(D) : D in basis]>; // [ 2, 2, 120 ]
      assert Invariants(global) eq  [ 2, 2, 120 ];    

           
  //////////////////////////////////////////////////////////////////////
  // "Sieve"
  //////////////////////////////////////////////////////////////////////  
  
  "There are", [#Places(Cp,i) : i in [1..3]], "places of degree 1, 2, and 3 over F_5"; // [8, 16, 32];    

  pl1 := Places(Cp,1);
  cubicDivisors :=
           {@ &+[Divisor(pl1[i])*tup[i] : i in [1..#pl1]] : tup in CartesianPower([0..3], #pl1) | &+[i : i in tup] eq 3@}
      join {@Divisor(D1) + Divisor(D2) : D1 in Places(Cp,1), D2 in Places(Cp,2)@}
      join {@Divisor(D) : D in Places(Cp,3)@};
    
  validCubicImages := {@@};
  for D0 in cubicDivisors do
      D := D0 - Degree(D0)*basePtp;
      if Inverse(mPic)(D) in global then
        validCubicImages := 
        validCubicImages join {@Inverse(mPic)(D)@};
      end if;
  end for;
  #validCubicImages, "of them are in the image of Abel--Jacobi"; // 86  
  "There are", Binomial(#pts + 3 - 1, 3) + #pts * #eqns2, "rational cubic divisors, so we have found all of the points!";
