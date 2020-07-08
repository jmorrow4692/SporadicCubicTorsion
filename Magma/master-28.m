//////////////////////////////////////////////////////////////////////////////////////////
//  This is a complete determination of the cubic points on X1(28)
//////////////////////////////////////////////////////////////////////////////////////////
 
/****************************************************************************** 
 Here is a summary of the argument.
 
 - X_1(28) has genus 10, and rank 0.
 - There are 9 rational cusps, 3 quadratic cusps, 1 cubic cusp, and 2 sextic cusps
 - The torsion subgroup has invariants [2, 4, 12, 936] (from the Hecke bound)
 - Given the complexity of the model, we use modular units to
   find, explicitly, the cusps (and as usual hard code the defining ideals).
 - Mod 3 there are 9 degree 1 places (all reductions of cusps), 
   3 degree 2 places (all reductions of cusps), and 
   5 degree 3 places. We check that the reduction of the cubic cusp is the only 
   degree 3 place which meets the reduction of J(Q)
   
******************************************************************************/

  //////////////////////////////////////////////////////////////////////////////////////////
  // Input the homebrewed functions 
  //////////////////////////////////////////////////////////////////////////////////////////
  
  load "functions.m";

  //////////////////////////////////////////////////////////////////////////////////////////
  // Defining equations from Sutherland, which are the modular units
  //////////////////////////////////////////////////////////////////////////////////////////

A2<x,y> := AffineSpace(Rationals(),2);

Def := <
[],
(-x^9*y^6 + x^8*y^6 + 2*x^8*y^5 - x^7*y^6 - x^7*y^5 - x^7*y^4 + 2*x^6*y^5 - 
    x^6*y^4 - x^5*y^4 + x^5*y^3)/(x^5*y^5 - 3*x^5*y^4 + 3*x^5*y^3 - x^5*y^2 - 
    x^4*y^6 - 4*x^4*y^5 + 24*x^4*y^4 - 38*x^4*y^3 + 25*x^4*y^2 - 6*x^4*y + 
    7*x^3*y^6 - 18*x^3*y^5 + 2*x^3*y^4 + 33*x^3*y^3 - 36*x^3*y^2 + 13*x^3*y - 
    x^3 - 9*x^2*y^6 + 34*x^2*y^5 - 43*x^2*y^4 + 12*x^2*y^3 + 17*x^2*y^2 - 
    14*x^2*y + 3*x^2 + 2*x*y^6 - 7*x*y^5 + 5*x*y^4 + 10*x*y^3 - 20*x*y^2 + 
    13*x*y - 3*x + y^6 - 6*y^5 + 15*y^4 - 20*y^3 + 15*y^2 - 6*y + 1),
(-x^4*y^3 + x^4*y^2 + 3*x^3*y^3 - 4*x^3*y^2 + x^3*y - 4*x^2*y^3 + 7*x^2*y^2 - 
    3*x^2*y + 3*x*y^3 - 7*x*y^2 + 5*x*y - x - y^3 + 3*y^2 - 3*y + 1)/(x^5*y^3 - 
    2*x^4*y^2 + x^3*y),
(-x^2*y^2 + x^2*y + 2*x*y^2 - 3*x*y + x - y^2 + 2*y - 1)/(x^3*y^2 - x^2*y),
(x^3*y^3 - 2*x^3*y^2 + x^3*y - 3*x^2*y^3 + 7*x^2*y^2 - 5*x^2*y + x^2 + 3*x*y^3 -
    8*x*y^2 + 7*x*y - 2*x - y^3 + 3*y^2 - 3*y + 1)/(x^5*y^3 - 2*x^4*y^2 + x^3*y),
(-y + 1)/(x*y),
(-y^2 + 2*y - 1)/(x^2*y^2 - x*y),
(-y^2 + 2*y - 1)/(x^3*y^2 - x^2*y),
(-y^2 + 2*y - 1)/(x^3*y^3 - x^2*y^2),
x-y+1, 
x^2*y-x*y^2+y-1,
x-y,
x^3*y-x^2*y^2-x^2*y+x*y^2-y+1,
x^2*y-x*y^2-x*y+y^2-1,
x^2*y-x^2-x*y^2+x*y-x+y-1,
x^2*y^2-2*x^2*y-x*y^3+2*x*y^2-y+1,
x^4*y-x^3*y^3-x^3*y+x^2*y^4+x^2*y-x^2-x*y^4+x*y^3-x*y^2+x*y+y^3-2*y^2+y,
x^4*y-3*x^3*y^2+x^3*y+2*x^2*y^3+x^2*y^2-2*x^2*y-2*x*y^3+2*x*y^2+y^2-2*y+1,
x^5*y^2-3*x^4*y^3+3*x^3*y^4+2*x^3*y^2-2*x^3*y-x^2*y^5-x^2*y^4-x^2*y^3+x^2*y^2+x^2*y+x*y^5-x*y^4+2*x*y^3-2*x*y^2-x*y+x-y^4+2*y^3-y^2,
x^4*y-2*x^3*y^2-x^3*y+x^2*y^3+x^2*y^2+2*x^2*y-x^2-2*x*y^2+x*y+y-1,
x^5*y^2-3*x^4*y^3-x^4*y^2+3*x^3*y^4+3*x^3*y^3-x^2*y^5-3*x^2*y^4-2*x^2*y^2+3*x^2*y+x*y^5+3*x*y^3-4*x*y^2-y^4+y^3+y-1,
x^5*y^2-x^5*y-2*x^4*y^3+2*x^4*y+x^3*y^4+3*x^3*y^3-4*x^3*y^2+x^3-2*x^2*y^4+x^2*y^3+x^2*y^2-x^2*y+x*y^4-x*y^3+x*y^2-x*y-y^3+2*y^2-y,
x^7*y^4 - 2*x^7*y^3 + x^7*y^2 - 3*x^6*y^5 + 4*x^6*y^4 + 3*x^5*y^6 - 5*x^5*y^4 - 
    3*x^5*y^3 + 3*x^5*y^2 - x^5*y - x^4*y^7 - 4*x^4*y^6 + 5*x^4*y^5 + 6*x^4*y^4 
    - x^4*y^3 - x^4*y^2 - x^4*y + 2*x^3*y^7 - 3*x^3*y^5 - 6*x^3*y^4 + 2*x^3*y^3 
    + 5*x^3*y^2 - x^3*y - x^2*y^7 + 3*x^2*y^5 + 2*x^2*y^4 - x^2*y^3 - 7*x^2*y^2 
    + 4*x^2*y + x*y^6 - 2*x*y^5 - 3*x*y^4 + 8*x*y^3 - 4*x*y^2 + y^3 - 3*y^2 + 
    3*y - 1,
x^5*y^4 - 3*x^5*y^3 + 3*x^5*y^2 - 2*x^4*y^5 + 5*x^4*y^4 - 5*x^4*y^3 + x^3*y^6 - 
    x^3*y^5 + x^3*y^4 - x^3*y^3 + 4*x^3*y^2 - 3*x^3*y - x^2*y^6 + x^2*y^5 - 
    x^2*y^4 + x^2*y^3 - x^2*y^2 + x^2*y + 2*x*y^5 - 5*x*y^4 + 5*x*y^3 - 2*x*y^2 
    - x*y + x - y^4 + 3*y^3 - 3*y^2 + y,
x^8*y^2 + x^7*y^5 - 4*x^7*y^4 - x^7*y^2 - 3*x^6*y^6 + 11*x^6*y^5 - 2*x^6*y^4 + 
    2*x^6*y^2 - 2*x^6*y + 3*x^5*y^7 - 9*x^5*y^6 - 4*x^5*y^5 + 4*x^5*y^4 - 
    2*x^5*y^3 + 4*x^5*y^2 - x^4*y^8 + x^4*y^7 + 9*x^4*y^6 - 5*x^4*y^5 + 
    6*x^4*y^4 - 15*x^4*y^3 + 7*x^4*y^2 - 2*x^4*y + x^4 + x^3*y^8 - 4*x^3*y^7 + 
    3*x^3*y^6 - 14*x^3*y^5 + 26*x^3*y^4 - 15*x^3*y^3 + 7*x^3*y^2 - 5*x^3*y + x^3
    - x^2*y^7 + 8*x^2*y^6 - 13*x^2*y^5 + 12*x^2*y^4 - 18*x^2*y^3 + 19*x^2*y^2 - 
    8*x^2*y + x^2 - 4*x*y^5 + 13*x*y^4 - 16*x*y^3 + 10*x*y^2 - 4*x*y + x + y^4 -
    4*y^3 + 6*y^2 - 4*y + 1,
x^6*y^2 - 5*x^5*y^3 + 2*x^5*y^2 + x^4*y^5 + 5*x^4*y^4 - 2*x^4*y^2 - x^4*y - 
    2*x^3*y^6 - 4*x^3*y^4 + 9*x^3*y^2 - 4*x^3*y + x^2*y^7 + 5*x^2*y^4 - 
    11*x^2*y^3 + 4*x^2*y^2 + x^2*y - x*y^7 + 2*x*y^6 - 4*x*y^5 + 6*x*y^4 - 
    3*x*y^3 + x*y^2 - 2*x*y + x + y^6 - 3*y^5 + 3*y^4 - y^3,
x^10*y^3 - 6*x^9*y^4 + 15*x^8*y^5 + 3*x^8*y^3 - 3*x^8*y^2 - 19*x^7*y^6 - 
    6*x^7*y^5 - 3*x^7*y^4 + 5*x^7*y^3 + 3*x^7*y^2 + 12*x^6*y^7 + 18*x^6*y^6 - 
    6*x^6*y^5 + 3*x^6*y^4 - 12*x^6*y^3 - 3*x^6*y^2 + 3*x^6*y - 3*x^5*y^8 - 
    18*x^5*y^7 + 6*x^5*y^6 - 12*x^5*y^5 + 27*x^5*y^4 - 3*x^5*y^3 - x^5*y^2 - 
    2*x^5*y + 6*x^4*y^8 + 3*x^4*y^7 + 10*x^4*y^6 - 24*x^4*y^5 + 3*x^4*y^4 - 
    2*x^4*y^3 + 5*x^4*y^2 + x^4*y - x^4 - 3*x^3*y^8 - 3*x^3*y^7 + 15*x^3*y^5 - 
    6*x^3*y^4 - 3*x^3*y^3 - x^3*y^2 + x^3*y + 6*x^2*y^7 - 12*x^2*y^6 + 
    10*x^2*y^5 - 13*x^2*y^4 + 15*x^2*y^3 - 7*x^2*y^2 + x^2*y - 4*x*y^6 + 
    13*x*y^5 - 16*x*y^4 + 10*x*y^3 - 4*x*y^2 + x*y + y^5 - 4*y^4 + 6*y^3 - 4*y^2
    + y,
x^8*y^2 - 5*x^7*y^3 + 10*x^6*y^4 + 2*x^6*y^2 - x^6*y - 10*x^5*y^5 - 10*x^5*y^3 +
    8*x^5*y^2 - x^5*y + 5*x^4*y^6 + x^4*y^5 + 13*x^4*y^4 - 9*x^4*y^3 - x^4*y^2 -
    x^4*y - x^3*y^7 - 2*x^3*y^6 - 4*x^3*y^5 - 4*x^3*y^4 + 8*x^3*y^3 + 2*x^3*y^2 
    - x^3*y + x^2*y^7 - x^2*y^6 + 8*x^2*y^5 - 11*x^2*y^4 + 8*x^2*y^3 - 
    10*x^2*y^2 + 5*x^2*y - 2*x*y^6 + 5*x*y^5 - 10*x*y^4 + 14*x*y^3 - 8*x*y^2 + 
    x*y + y^5 - 3*y^4 + 4*y^3 - 4*y^2 + 3*y - 1
>;

 //////////////////////////////////////////////////////////////////////
  // Example X1(28)
  //////////////////////////////////////////////////////////////////////


  X28 := Curve(A2,Def[28]);
  F28 := FunctionField(X28);
  Xp<[Z]> := ProjectiveClosure(X28);
  pts := PointSearch(Xp,100);
  rationalPlaces := {@pl: pl in Places(pt), pt in pts@};

  torsionData := {@@}; 
  for p in [3,5] do
      invs := Invariants(ClassGroup(Curve(Reduction(Xp,p))));
      torsionData := Include(torsionData, invs);
      <p,torsBound(torsionData)>;      
  end for;                

  "The rational torsion subgroup is a subgroup of", torsBound(torsionData); // [ 936, 24, 4, 4 ]
      
  PX<[Z]>  := AmbientSpace(Xp);
  CX<[Z]> := CoordinateRing(PX);

//////////////////////////////////////////////////////////////////////
// quadratic cusps
//////////////////////////////////////////////////////////////////////


//CC1 := Zeros(F28!(Def[7]))[2]; 
//D21 := Divisor(CC1);
I21 := Ideal([
    Z[2]*Z[3]^2 - Z[3]^3,
    Z[1]
]);


//CC2 := Zeros(F28!(Def[7]))[4];
//D22 := Divisor(CC2);
I22 := Ideal([
    Z[1]^2 - 2*Z[1]*Z[3] + 2*Z[3]^2,
    Z[1]*Z[2] - Z[1]*Z[3],
    Z[2]*Z[3] - Z[3]^2
]);

//CC3 := Poles(F28!(Def[7]))[2];
//D23 := Divisor(CC3);
I23 := Ideal([
    Z[2]^2*Z[3]^2 + Z[3]^4,
    Z[1]
]);

//////////////////////////////////////////////////////////////////////
// cubic cusp
//////////////////////////////////////////////////////////////////////

//CC4 := (Poles(F28!(Def[7]))[3]);
//D34 := Divisor(CC4);
I34 := Ideal([
    Z[2]^2*Z[3] - Z[1]*Z[3]^2 + Z[2]*Z[3]^2 - 2*Z[3]^3,
    Z[1]^2 + 2*Z[1]*Z[3] - Z[2]*Z[3] - Z[3]^2,
    Z[1]*Z[2] - Z[3]^2
]);

//////////////////////////////////////////////////////////////////////
// sextic cusps
//////////////////////////////////////////////////////////////////////

//CC5 := (Poles(F28!(Def[2]))[4]);
//D65 := Divisor(CC5);
I65 := Ideal([
    Z[2]^3*Z[3] - 6*Z[1]^2*Z[3]^2 + 34*Z[1]*Z[2]*Z[3]^2 - 
        21*Z[2]^2*Z[3]^2 - 8*Z[1]*Z[3]^3 - 39*Z[2]*Z[3]^3 + 
        31*Z[3]^4,
    Z[1]^3 - 3*Z[1]^2*Z[3] - Z[1]*Z[2]*Z[3] + Z[1]*Z[3]^2 + 
        2*Z[3]^3,
    Z[1]^2*Z[2] - 2*Z[1]^2*Z[3] + Z[1]*Z[2]*Z[3] - 2*Z[2]^2*Z[3] + 
        Z[1]*Z[3]^2 - 4*Z[2]*Z[3]^2 + 4*Z[3]^3,
    Z[1]*Z[2]^2 - 3*Z[1]^2*Z[3] + 11*Z[1]*Z[2]*Z[3] - 8*Z[2]^2*Z[3] 
        - 2*Z[1]*Z[3]^2 - 14*Z[2]*Z[3]^2 + 12*Z[3]^3
]);

//CC6 := (Zeros(F28!(Def[3]))[6]); 
//D66 := Divisor(CC6);
I66 := Ideal([
    Z[2]^4 - 5*Z[2]^3*Z[3] + Z[1]^2*Z[3]^2 + 11*Z[2]^2*Z[3]^2 + 
        Z[1]*Z[3]^3 - 13*Z[2]*Z[3]^3 + 7*Z[3]^4,
    Z[1]^3 - Z[2]^3 + 5*Z[2]^2*Z[3] + 6*Z[1]*Z[3]^2 - 11*Z[2]*Z[3]^2
        + 6*Z[3]^3,
    Z[1]*Z[2] - Z[2]*Z[3] + Z[3]^2
]);

//////////////////////////////////////////////////////////////////////
// Reducing everything modulo 3
//////////////////////////////////////////////////////////////////////

  p := 3;
  Cp<[T]> := Curve(Reduction(Xp,p));
  kCp := FunctionField(Cp);  
  pic,mPic := ClassGroup(Cp);  
  "There are", [#Places(Cp,i) : i in [1..3]], "places of degree 1, 2, and 3 over F_3"; // [9, 3, 5];  

  rationalPlaces := [pl : pl in Places(Cp,1)];

// Reducing the cusps

P<[Z]>  := AmbientSpace(Cp);
CC<[Z]> := CoordinateRing(P);

//////////////////////////////////////////////////////////////////////
// reducing quadratic cusps
//////////////////////////////////////////////////////////////////////

I1p,m1p := Reduction(I21,p);
I1p := Ideal([CC!Basis(I1p)[1],CC!Basis(I1p)[2]]);
D1p:= Divisor(Cp,Cp meet Scheme(P,I1p));
DD1p := Support(D1p)[3];
assert (Degree(DD1p) eq 2);

I2p,m2p := Reduction(I22,p);
I2p := Ideal([CC!Basis(I2p)[1],CC!Basis(I2p)[2],CC!Basis(I2p)[3]]);
D2p:= Divisor(Cp,Cp meet Scheme(P,I2p));
DD2p := Support(D2p)[3];
assert (Degree(DD2p) eq 2);

I3p,m3p := Reduction(I23,p);
I3p := Ideal([CC!Basis(I3p)[1],CC!Basis(I3p)[2]]);
D3p:= Divisor(Cp,Cp meet Scheme(P,I3p));
DD3p := Support(D3p)[2];
assert (Degree(DD3p) eq 2);

//////////////////////////////////////////////////////////////////////
// reducing cubic cusps
//////////////////////////////////////////////////////////////////////

I4p,m4p := Reduction(I34,p);
I4p := Ideal([CC!Basis(I4p)[1],CC!Basis(I4p)[2],CC!Basis(I4p)[3]]);
D4p:= Divisor(Cp,Cp meet Scheme(P,I4p));
DD4p := Support(D4p)[3];
assert (Degree(DD4p) eq 3);

//////////////////////////////////////////////////////////////////////
// reducing sextic cusps
//////////////////////////////////////////////////////////////////////

I5p,m5p := Reduction(I65,p);
I5p := Ideal([CC!Basis(I5p)[1],CC!Basis(I5p)[2],CC!Basis(I5p)[3],CC!Basis(I5p)[4]]);
D5p:= Divisor(Cp,Cp meet Scheme(P,I5p));

DD51p := Support(D5p)[3];
assert (Degree(DD51p) eq 3);
DD52p := Support(D5p)[4];
assert (Degree(DD52p) eq 3);

I6p,m6p := Reduction(I66,p);
I6p := Ideal([CC!Basis(I6p)[1],CC!Basis(I6p)[2]]);
D6p:= Divisor(Cp,Cp meet Scheme(P,I6p));

DD6p := (Support(D6p)[2]);
assert (Degree(DD6p) eq 6);




// To summarize, 
// the degree 1 places come from reductions of rational cusps
// the degree 2 places come from reductions of quadratic cusps
// 3 of the degree 3 places come from reduction of cubic and sextic cusp

  basePt := Places(Cp,1)[1];
  divs := {@pl - Degree(pl)*basePt : pl in Places(Cp,1)@} join {@pl - Degree(pl)*basePt : pl in Places(Cp,2)@} join {DD4p - Degree(DD4p)* basePt}; 

  global, mGlobal := 
     sub<pic | [(Inverse(mPic))(divs[i]) : i in [1..#divs] ]>;  
  assert Invariants(global) eq [ 2,4,12,936];
   
  validCubicImages := {@@};
  for pl in Places(Cp,3) do
      D := Divisor(pl) - Degree(pl)*basePt;
      if Inverse(mPic)(D) in global then
        validCubicImages := 
        validCubicImages join {@Inverse(mPic)(D)@};
      end if;
  end for;
  #validCubicImages, "of them are in the image of Abel--Jacobi"; // 1

  // This tells us that the only cubic point on X1(28) is the known cubic cusp. 
