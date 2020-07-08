//////////////////////////////////////////////////////////////////////////////////////////
// This is a complete determination of the cubic points on X1(30)
//////////////////////////////////////////////////////////////////////////////////////////
 
 N := 30; 

/****************************************************************************** 
 Here is a summary of the argument.
 
 We map to X0(30), which is  a genus 3 hyperelliptic curve. Since the rational
 points of X0(30) are cusps, any non cuspidal cubic point of X1(30) maps 
 to a cubic point of X0(30). The map X0(30)^(3) --> J0(30) is not injective, but
 is still injective on rational cubic points. 

 The rational torsion of J0(30) is Z/2Z x Z/4Z x Z/24Z. The local bound on the torsion is [2,2,4,24];
 Since X_0(30) is hyperelliptic, we rule out the extra two torsion point via the 
 Galois action on the Weierstrass points (via the command twoTorsionRank)
 Computing the preimage of Abel--Jabobi (and a small argument to verify that 
 the cubic points on X0(30) do not lift to X1(30) succeeds.
******************************************************************************/


  //////////////////////////////////////////////////////////////////////////////////////////
  // Input the homebrewed functions 
  //////////////////////////////////////////////////////////////////////////////////////////
  
  load "functions.m";

  
  //////////////////////////////////////////////////////////////////////////////////////////
  // We double checked that this model is isomorphic to the model from the peer reviewed code
  // of Ozman--Siksek
  //////////////////////////////////////////////////////////////////////////////////////////  
  
  X := SmallModularCurve(N);

  
  //////////////////////////////////////////////////////////////////////
  // Compute the local torsion bound
  //////////////////////////////////////////////////////////////////////   

  torsData := {@@};  
//  for p in [q : q in PrimesUpTo(40) | not q in PrimeDivisors(2*N) ] do
  for p in [7,23 ] do     
      invs := Invariants(ClassGroup(Curve(Reduction(X,p))));
      torsData := torsData join {@invs@};
      <p,invs>;
  end for;

  "The rational torsion subgroup is a subgroup of", torsBound(torsData); ; // [24,4,2,2]
  
  "The Hecke bounds does not improve this";     
  
     simpX, m := SimplifiedModel(X);
     mInv := Inverse(m);

  "The F2-rank of J[2] is", twoTorsionRank(HyperellipticPolynomials(simpX)); 

  ////////////////////////////////////////////////////////////////////////
  // Compute the known small degree points
  ////////////////////////////////////////////////////////////////////////
     

     pts := RationalPoints(simpX : Bound := 50);
     D1 :=   (Divisor(mInv(pts[1])) - Divisor(mInv(pts[5])));  
       order(D1); // 24
     D2 := 3*(Divisor(mInv(pts[1])) - Divisor(mInv(pts[2])));
       order(D2); // 4
     D3 := 3*(Divisor(mInv(pts[1])) - Divisor(mInv(pts[4])));
       order(D3); // 2
  
     exists{[a,b,c] : a in [0..order(D1)-1], b in [0..order(D2)-1], c in [0..order(D3)-1]
                    | not [a,b,c] eq [0,0,0] and  IsPrincipal(a*D1 + b*D2 + c*D3)};
  // false
  // so D1, D2, D3 generate as the F2-rank of J[2] is 3


  ////////////////////////////////////////////////////////////////////////
  // Now determine the image of Abel--Jacobi
  ////////////////////////////////////////////////////////////////////////        
     
  basePt := Divisor(X![1,0,0]);      
  cubicPts := {@@};		      
  for a in [0..order(D1)-1], b in [0..order(D2)-1], c in [0..order(D3)-1] do
      D := a*D1 + b*D2 + c*D3;
      R,mR := RiemannRochSpace(D + 3*basePt);
      sup := Support(D + 3*basePt + Divisor(mR(R.1)));
      if 3 in [Degree(D) : D in sup] then
	  cubicPts := cubicPts join {@sup[1]@};
      end if;      
  end for;
  
    
  ////////////////////////////////////////////////////////////////////////
  // Compute the j invariants of the cubic points and check to see if 
  // there is a 30-torsion point on a twist of a corresponding curve.
  ////////////////////////////////////////////////////////////////////////        

  j := jInvariant(X,30);
  
  num := Numerator(j);
  den := Denominator(j);
  P1 := Curve(ProjectiveSpace(Rationals(),1));
  j0 := map<X -> P1 | [num,den]>;


  // 0 in {Norm(j0(RepresentativePoint(pt))[2]) : pt in cubicPts};
  // All non-cuspidal, but we knew that from the explicit description of the cusps
  
  // Check for CM j-invariants
  if not &and( [IsIntegral(j0(RepresentativePoint(pt))[1]) : pt in cubicPts])
    then
      "All j-invariants are non-integral, thus not CM";
  end if;

  // For each j, reduce mod ell until we find an ell such that no twist has a torsion point.
  for pt in cubicPts do   
      jpt := j0(RepresentativePoint(pt))[1]; 
      K := Parent(RepresentativePoint(pt)[1]); 
      hasTwistWithNTorsion(EllipticCurveFromjInvariant(K!jpt),N, 1000);
  end for;
