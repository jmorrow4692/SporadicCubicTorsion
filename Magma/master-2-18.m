//////////////////////////////////////////////////////////////////////////////////////////
// This is a complete determination of the cubic points on X1(2,18)
//////////////////////////////////////////////////////////////////////////////////////////

 N := 18;

/******************************************************************************
 Here is a summary of the argument.

   - X1(2,18) is a genus 7 curve
   - Torsion:
     - Checking locally gives a bound on the torsion of Z/126Z x Z/42Z x Z/6Z
     - We can find a subgroup of rational torsion of Z/126Z x Z/42Z x Z/2Z
     - Hecke bound tells us that this is all torsion
   - Known points
     - there are 9 known rational cusps, 3 quadratic cusps, and 3 cubic cusps
     - (9 + 3 -1 choose 3) + 9*3 + 3 = 195 rational degree 3 divisors	 
   - We compute the intersection of the image of the Abel-Jacobi with our known subgroup mod 5, and this gives 195 cubic divisors
   - At the end, we check that the cubic points are indeed cusps

******************************************************************************/

  //////////////////////////////////////////////////////////////////////////////////////////
  // Input the homebrewed functions
  //////////////////////////////////////////////////////////////////////////////////////////

  load "functions.m";

  //////////////////////////////////////////////////////////////////////////////////////////
  // Model of X(2,18) from Derickx--Sutherland
  ////////////////////////////////////////////////////////////////////////////////////////// 

m := 2;
n := 9;
g := 7;
F := Rationals();
R<u,v>:=PolynomialRing(F,2);
X := 3*v^6 + 6*u^2*v^5 + (4*u^4 - 3)*v^4 + (u^6 - u^4 - 3*u^2 - 1)*v^3 - 3*u^4*v^2 - 2*u^4*v - u^4;
q := 1/u;
t := 4*u*v/(u^2-u^2*v-3*v-1);
E := [0,t^2-2*q*t-2,0,-(t^2-1)*(q*t+1)^2,0];
P := [R!0,0];
Q := [(t+1)*(q*t+1),t*(q*t+1)*(t+1)];
j218 := jInvariant(EllipticCurve(E));

  A2<u,v> := AffineSpace(R);
    X<[T]> := ProjectiveClosure(Curve(A2,X));
    "Genus is", Genus(X); // 7

 //////////////////////////////////////////////////////////////////////
  // Get the canonical model.
  // This presents it as a complete intersection of 3 quadratics in P4.
  //////////////////////////////////////////////////////////////////////
    
    phi := CanonicalMap(X);
    Xsm<[T]> := CanonicalImage(Domain(phi),phi); 
    phi := map<X -> Xsm | DefiningEquations(phi)>;

  //////////////////////////////////////////////////////////////////////
  // Compute the local torsion bound
  //////////////////////////////////////////////////////////////////////   

  torsData := [];
  for p in [5,7] do     
      invs := Invariants(ClassGroup(Curve(Reduction(Xsm,p))));
      torsData := torsData cat [invs];
      <p,invs>;
  end for;

    "Torsion bound is", torsBound(torsData);
    // [ 126, 42, 6 ]; additional primes do not improve the bound


  ////////////////////////////////////////////////////////
  // Create the global torsion
  ////////////////////////////////////////////////////////

  // pts1 := PointSearch(Xsm,100);

  // global rational points, all cusps

  // hard coded
  basePt := [-3,1,0,0,0,0,0];

  pts1 := 
  {@
	[-3,1,0,0,0,0,0],
  	[9,-2,0,0,3,-2,1],
	[-3,2,-1,0,1,0,0],  
	[9,0,3,-1,6,-3,1],
	[14,-2,2,-1,5,-3,1],
        [-4,-1,0,0,-3,1,0],
	[14,-2,-2,1,5,-3,1],
	[9,0,-3,1,6,-3,1],
	[-3,2,1,0,1,0,0] 
  @};

  // Find the ideals of the quadratic cusps
  //  pts2 := divisorSearch(Xsm,2,1);
    
  //  #pts2, "quadratic points"; // 3    

// Hard code the ideals defining the quadratic points
    eqns2 := {@
    [
        T[6]^2 + 4*T[6]*T[7] + 13/3*T[7]^2,
        T[1] + 19/2*T[6] + 17/2*T[7],
        T[2] - T[6] + T[7],
        T[3] - 3/2*T[6] - 1/2*T[7],
        T[4] - T[7],
        T[5] + 5/2*T[6] + 5/2*T[7]
    ],
    [
        T[6]^2 + 4*T[6]*T[7] + 13/3*T[7]^2,
        T[1] + 19/2*T[6] + 17/2*T[7],
        T[2] - T[6] + T[7],
        T[3] + 3/2*T[6] + 1/2*T[7],
        T[4] + T[7],
        T[5] + 5/2*T[6] + 5/2*T[7]
    ],
    [
        T[5]^2 + 5*T[5]*T[6] + 7*T[6]^2,
        T[1] + T[5] + 6*T[6],
        T[2] - T[5] - 2*T[6],
        T[3],
        T[4],
        T[7]
    ]
@};

// Find the ideals of the cubic cusps
//    pts3 := divisorSearch(Xsm,3,1);
//    #pts3; // 1

eqns3 := {@
    [
        T[2]^2 - 4/3*T[2]*T[7] - 89/3*T[6]*T[7] - 81*T[7]^2,
        T[2]*T[6] + 11/3*T[2]*T[7] + 31/3*T[6]*T[7] + 28*T[7]^2,
        T[6]^2 - 1/3*T[2]*T[7] + 7/3*T[6]*T[7] - T[7]^2,
        T[1] + 4*T[2] + 13*T[6] + 29*T[7],
        T[3],
        T[4],
        T[5] + 4*T[6] + 6*T[7]
    ],
    [
    T[4]^2 + 6*T[4]*T[7] + 6*T[6]*T[7] + 15*T[7]^2,
    T[4]*T[6] + T[4]*T[7] - 3*T[6]*T[7] - 7*T[7]^2,
    T[6]^2 + 2/3*T[4]*T[7] + 8*T[6]*T[7] + 13*T[7]^2,
    T[1] - 5*T[4] + 1/2*T[6] - 17/2*T[7],
    T[2] + 2*T[4] + 2*T[6] + 6*T[7],
    T[3] + T[4] - 3/2*T[6] - 9/2*T[7],
    T[5] + 5/2*T[6] + 3/2*T[7]
    ],
    [
    T[4]^2 - 6*T[4]*T[7] + 6*T[6]*T[7] + 15*T[7]^2,
    T[4]*T[6] + T[4]*T[7] + 3*T[6]*T[7] + 7*T[7]^2,
    T[6]^2 - 2/3*T[4]*T[7] + 8*T[6]*T[7] + 13*T[7]^2,
    T[1] + 5*T[4] + 1/2*T[6] - 17/2*T[7],
    T[2] - 2*T[4] + 2*T[6] + 6*T[7],
    T[3] + T[4] + 3/2*T[6] + 9/2*T[7],
    T[5] + 5/2*T[6] + 3/2*T[7]
]
@};

  //////////////////////////////////////////////////////////////////////////////////////////
  // The divisorSearch only returned 1 cubic point
  // But when we computed the image of Abel--Jacobi mod 5, we found 2 more cubic points.
  // To see if they lift to Q, we worked "backwards": we computed which global torsion point D
  // a given cubic point mapped to, and then checked whether D is the image of a cubic point
  // over Q. 
  //////////////////////////////////////////////////////////////////////////////////////////


"There are", Binomial(#pts1 + 3 - 1, 3) + #pts1 * #eqns2 + #eqns3, "rational cubic divisors";    


  //////////////////////////////////////////////////////////////////////
  // Verify that the known points generate [2,42,126]
  //////////////////////////////////////////////////////////////////////    

divs1X := {@ Divisor(Xsm!pt) : pt in pts1@};
divs2X := {@Divisor(Xsm,Scheme(Xsm,eqn)) : eqn in eqns2 @};
divs3X := {@Divisor(Xsm,Scheme(Xsm,eqn)) : eqn in eqns3 @};

      basePtX := Divisor(Xsm!basePt);

    basisX := [
	  (divs1X[4] - basePtX),
	  divs1X[7] + divs1X[8] - 2*basePtX,    	  
	  7*(divs3X[1] - 3*basePtX)	  
      ]; 

/*
    order(basisX[1]); // 126
    order(basisX[2]); // 42
    order(basisX[3]); // 2

  time exists(t){ <a,b,c> : a in [0..125], b in [0..41], c in [0,1] |
      not <a,b,c> eq <0,0,0>
      and IsPrincipal(a*basisX[1] + b*basisX[2] + c*basisX[3]) };
*/

 // We check that this generates torsion mod 5

   //////////////////////////////////////////////////////////////////////////////////////////
  // We compute the intersection of the image of the Abel-Jacobi with our known subgroup
  // mod 5, and this gives 195 cubic divisors.
  //////////////////////////////////////////////////////////////////////////////////////////
  p := 5;    
  data := {@@};  
    Xp<[t]> := Curve(Reduction(Xsm,p));
      A, mA := ClassGroup(Xp);
      divs1Xp := {@ Divisor(Xp![GF(p)!cpt : cpt in pt]) : pt in pts1@};
      divs2Xp := {@Divisor(Xp,Scheme(Xp,[Parent(t[2])!e : e in eqn])) : eqn in eqns2 @};
      divs3Xp := {@Divisor(Xp,Scheme(Xp,[Parent(t[2])!e : e in eqn])) : eqn in eqns3 @};
      basePtXp := Divisor(Xp!basePt);

    global, mGlobal := sub<A | [Inverse(mA)(D - Degree(D)*basePtXp) : D in divs1Xp join divs2Xp join divs3Xp]>;
      Invariants(global);
      // [2, 42, 126]     

   basisXp := [
	  (divs1Xp[4] - basePtXp),
	  divs1Xp[7] + divs1Xp[8] - 2*basePtXp,    	  
	  7*(divs3Xp[1] - 3*basePtXp)	  
      ]; 
      [Order(Inverse(mA)(D)) : D in basisXp];
      //  [ 126, 42, 2 ]

    global, mGlobal := sub<A | [Inverse(mA)(D) : D in basisXp]>;
      Invariants(global);
      // [ 126, 42, 2 ]
 

  //////////////////////////////////////////////////////////////////////
  // "Sieve"
  //////////////////////////////////////////////////////////////////////  

  B := AbelianGroup([126,42,2]);  
  mBtoA := hom<B -> A | [(Inverse(mA))(basisXp[i]) : i in [1..3]]>;
    
   
    tups := CartesianProduct( [[0..Order(A.i)-1] : i in [1..#Generators(A)-1] ] );     
  
    B := AbelianGroup([126,42,2]);  
    mBtoA := hom<B -> A | [(Inverse(mA))(basisXp[i]) : i in [1..#basisXp]]>;  

    img := sub<A | Image(mBtoA)>;

    places1 := Places(Xp,1);

    divsp := 
         {@&+[tup[i]*places1[i] : i in [1..#places1]] : tup in tuplesOfDegree(#places1,3)@}
    join {@Divisor(D1) + Divisor(D2) : D1 in Places(Xp,1), D2 in Places(Xp,2)@}
    join {@Divisor(D) : D in Places(Xp,3)@};

    for D0 in divsp do         
      D := Inverse(mA)(D0 - 3*basePtXp);
      if D in img then 
          data := data join {@ Eltseq( Inverse(mBtoA)(D) ) @};
      end if; 
    end for; 

  "X^(3)(Q) has size at most", #data; // 195


  //////////////////////////////////////////////////////////////////////
  // Check that cubic points are all cusps
  ////////////////////////////////////////////////////////////////////// 

  cubicPts := {@Divisor(Xsm,Xsm meet Scheme(Xsm,I)) : I in eqns3@};
    // convert to places, so that we can use RepresentativePoint
    cubicPts := {@Support(D)[1] : D in cubicPts@};
  

  // Determine curve quotient from Xsm to X1(18)
  A := AutomorphismGroup(Xsm); // Time: 265.380
  X18,m18 := CurveQuotient(AutomorphismGroup(Xsm,[A.1]));
  
  // model from Sutherland for X1(18)
  A2<x,y> := AffineSpace(Rationals(),2);
  X := ProjectiveClosure(Curve(A2,y^2 + (x^3 - 2*x^2 + 3*x + 1)*y + 2*x)); 

  // Determine isomorphism between curve quotient and Sutherlands model
  boo,iso18 :=  IsIsomorphic(X18,X);   
  assert boo;

  // j-map from X1(18) --> X(1) via Sutherland
  r:=(x^2 - x*y - 3*x + 1)/((x-1)^2*(x*y+1));
  s:=(x^2 - 2*x - y)/(x^2 - x*y - 3*x - y^2 - 2*y);
  E:=[s-r*s+1,r*s-r^2*s,r*s-r^2*s,0,0];
  j18 := jInvariant(EllipticCurve(E));
  P1<s,t> := ProjectiveSpace(F,1);
  j := map<X -> P1 | [Numerator(j18),Denominator(j18)]>;
  
  // j-map from Xsm --> X1(18) --> X(1)
  j218 := m18*iso18*j;

  // Now check that all the cubic points are cusps
  j218(RepresentativePoint(cubicPts[1])); // (1:0)
  j218(RepresentativePoint(cubicPts[2])); // (1:0)
  j218(RepresentativePoint(cubicPts[3])); // (1:0)





