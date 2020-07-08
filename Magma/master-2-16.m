////////////////////////////////////////////////////////////////////////////////////////
// This is a complete determination of the cubic points on X1(2,16)
//////////////////////////////////////////////////////////////////////////////////////////
 
 N := 16; 

/****************************************************************************** 
 Here is a summary of the argument.

 We work directly on X := X(2,16); J(2,16) has rank 0.
 X has 8 rational points and 2 quadratic points (all cusps). 
 These give rise to 136 rational points of X^(3).

 The local bound on torsion is   [2, 2, 20, 20]
 The file mdsage/torsionComputations.sage gives a better bound of [2, 20, 20]
 The rational and quadratic cusps generate a subgroup isomorphic to [2, 20, 20]

 To finish, we compute the intersection of Abel--Jacobi with the mod 3 reduction 
 of the global torsion. 
******************************************************************************/

  load "functions.m";

  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Model of X(2,16) from Derickx--Sutherland
  //////////////////////////////////////////////////////////////////////////////////////////  
  
    F := Rationals();
    R<u,v>:=PolynomialRing(F,2);
    X := v^4 + (u^3 + 2*u)*v^3 + (u^4 - 2)*v^2 - (u^3 + 2*u)*v + 1;
    q := (2*u*v+v^2-1)/(v-1)^2;
    t := (1-v)^3*(u*v+u+2*v-2)*(2*u*v+v^2-1)/
         (2*u^4*v^3+6*u^3*v^4-6*u^3*v^2+4*u^2*v^5+8*u^2*v^4-24*u^2*v^3+8*u^2*v^2+4*u^2*v+u*v^6+
	 6*u*v^5-15*u*v^4+15*u*v^2-6*u*v-u+2*v^6-6*v^5+6*v^4-  4*v^3+6*v^2-6*v+2);
    E := [0,t^2-2*q*t-2,0,-(t^2-1)*(q*t+1)^2,0];
    P := [R!0,0];
    Q := [(t+1)*(q*t+1),t*(q*t+1)*(t+1)];
    
    A2<u,v> := AffineSpace(R);
    X<[T]> := ProjectiveClosure(Curve(A2,X));
    "Genus is", Genus(X); // 5

    
  //////////////////////////////////////////////////////////////////////
  // Get the canonical model.
  // This presents it as a complete intersection of 3 quadratics in P4.
  //////////////////////////////////////////////////////////////////////
    
    phi := CanonicalMap(X);
    Xsm<[T]> := CanonicalImage(Domain(phi),phi); 

    
  //////////////////////////////////////////////////////////////////////
  // Compute the local torsion bound
  //////////////////////////////////////////////////////////////////////   

    "Torsion bound is", Invariants(ClassGroup(Curve(Reduction(Xsm,3))));
    // [2, 2, 20, 20]; additional primes do not improve the bound

    // The file mdsage/torsionComputations.sage gives a better bound of [2, 20, 20]

    
    ////////////////////////////////////////////////////////////////////////
    // Compute the known small degree points
    ////////////////////////////////////////////////////////////////////////
    
    pts := PointSearch(Xsm, 1000);    
    #pts, "rational points";
    // returns 8 points
        
    // hard coded
    basePt := [0,1,-1,1,0];
    
    pts := 
    {@
	[ 0, 1, -1, 1, 0 ],
	[ -1, 1, 0, 0, 0 ],
	[ 0, 1, 1, 1, 0 ],
	[ 0, 1, 0, 0, 0 ],
	[ 1, 0, 0, 1, 0 ],
	[ 0, -1, -1, -1, 2 ],
	[ 0, 0, 0, 1, 0 ],
	[ 0, 1, -1, 1, 2 ]
    @};

    // Find the ideals of the quadratic cusps
    pts2 := divisorSearch(Xsm,2,1);
    
    #pts2, "quadratic points"; // 2    
    /*
      [
      Place at (1 : -$.1 - 1 : $.1 : $.1 + 1 : 1),
      Place at (-1 : $.1 + 1 : $.1 : -$.1 - 1 : 1)
      ]
   */

    // Hard code the ideals defining the quadratic points
    eqns2 := {@
    [
        2*T[4]^2 - 2*T[4]*T[5] + T[5]^2,
        T[1] - T[5],
        T[2] + T[4],
        T[3] - T[4] + T[5]
    ],
    [
        2*T[4]^2 + 2*T[4]*T[5] + T[5]^2,
        T[1] + T[5],
        T[2] + T[4],
        T[3] + T[4] + T[5]
    ]
    @};

   "There are", Binomial(#pts + 3 - 1, 3) + #pts * #eqns2, "rational cubic divisors";    

   
  //////////////////////////////////////////////////////////////////////
  // Verify that the known points generate [2,20,20]
  //////////////////////////////////////////////////////////////////////    

    p := 3;    

    Xp<[t]> := Curve(Reduction(Xsm,p));
      pic, mPic := ClassGroup(Xp);
      divs1Xp := {@ Divisor(Xp![GF(p)!cpt : cpt in pt]) : pt in pts@};
      divs2Xp := {@Divisor(Xp,Scheme(Xp,[Parent(t[2])!e : e in eqn])) : eqn in eqns2 @};
      basePtXp := Divisor(Xp!basePt);

    global, mGlobal := sub<pic | [Inverse(mPic)(D - Degree(D)*basePtXp) : D in divs1Xp join divs2Xp]>;
      Invariants(global);
      // [2, 20, 20]
          
    basisXp := [
	  5*(divs1Xp[3] - basePtXp),
	  divs1Xp[7] - basePtXp,
	  divs1Xp[8] - basePtXp	  
      ]; 
      [Order(Inverse(mPic)(D)) : D in basisXp];
      // [2, 20, 20]  

    global, mGlobal := sub<pic | [Inverse(mPic)(D) : D in basisXp]>;
      Invariants(global);
      // [2, 20, 20]
      

  //////////////////////////////////////////////////////////////////////////////////////////
  // We compute the intersection of the image of the Abel-Jacobi with our known subgroup
  // mod 3, and this gives 136 cubic divisors.
  // This takes about 16 seconds
  //////////////////////////////////////////////////////////////////////////////////////////       

  p := 3;    
  data := {@@};
  
   Xp<[t]> := Curve(Reduction(Xsm,p));
      pic, mPic := ClassGroup(Xp);
      basePtXp := Divisor(Xp!basePt);
      divs1Xp := {@ Divisor(Xp![GF(p)!cpt : cpt in pt]) : pt in pts@};      

  basisXp := [
	  5*(divs1Xp[3] - basePtXp),
	  divs1Xp[7] - basePtXp,
	  divs1Xp[8] - basePtXp	  
      ];       
     
    tups := CartesianProduct( [[0..Order(pic.i)-1] : i in [1..#Generators(pic)-1] ] );     
  
    B := AbelianGroup([2,20,20]);  
    mBtoPic := hom<B -> pic | [(Inverse(mPic))(basisXp[i]) : i in [1..#basisXp]]>;  

    img := sub<pic | Image(mBtoPic)>;

    places1 := Places(Xp,1);

    divsp := 
         {@&+[tup[i]*places1[i] : i in [1..#places1]] : tup in tuplesOfDegree(#places1,3)@}
    join {@Divisor(D1) + Divisor(D2) : D1 in Places(Xp,1), D2 in Places(Xp,2)@}
    join {@Divisor(D) : D in Places(Xp,3)@};

    for D0 in divsp do         
      D := Inverse(mPic)(D0 - 3*basePtXp);
      if D in img then 
          data := data join {@ Eltseq( Inverse(mBtoPic)(D) ) @};
      end if; 
    end for; 

  "X^(3)(Q) has size at most", #data; // 136

      
