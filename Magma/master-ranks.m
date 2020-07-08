//////////////////////////////////////////////////////////////////////////////////////////
// This code verifies the computation of the rank of Jzero(N) and JOne(N)
//////////////////////////////////////////////////////////////////////////////////////////

  genus0Gamma0 :=  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 25];
  genus0Gamma1 := [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];
  
  JZeroRankZero := 
  [1..36] cat  [38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 59, 60, 62, 63, 64, 66, 68, 69, 70, 71, 72, 75, 76, 78, 80, 81, 84, 87, 90  , 94, 95, 96, 98, 100, 104, 105, 108, 110, 119, 120, 126, 132, 140, 144, 150, 168, 180 ];
  
  JZeroPositiveRank := [ 57, 58, 65, 77, 82, 85, 88, 91, 92, 93, 99, 102, 112, 115, 117, 118, 121, 123, 124, 125, 128, 133, 135, 136, 138, 141, 142, 143, 145, 147,   152, 153, 155, 156, 160, 161, 162, 165, 169, 175, 177, 187, 188, 189, 190, 192, 196, 200, 203, 205, 207, 208, 209, 210, 213, 216, 217, 220, 221, 225, 235,   238, 240, 243, 245, 247, 252, 253, 261, 275, 280, 287, 288, 289, 295, 299, 300, 315, 319, 323, 329, 341, 343, 355, 357, 360, 361, 377, 391, 403, 413, 437,
  451, 475, 493, 497, 517, 527, 529, 533, 551, 589, 611, 649, 667, 697, 713, 767, 779, 781, 799, 833, 841, 893, 899, 923, 943, 961, 1003, 1081, 1121, 1189,
  1207, 1271, 1349, 1357, 1363, 1457, 1633, 1681, 1711, 1829, 1927, 2059, 2201, 2209, 2419, 2773, 2911, 3337, 3481, 4189 ];

  allInds := {N : N in JZeroRankZero} join {N : N in JZeroPositiveRank};
  
  JOnePositiveRank:= { 63, 80, 95, 104, 105, 126, 144 };

  JOneRankZero := [N : N in JZeroRankZero | not N in JOnePositiveRank];
  
  monsterPrimes := {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71};
  
  
  ////////////////////////////////////////////////////////////////////////////////
  // Check that we've covered every case
  ////////////////////////////////////////////////////////////////////////////////
  
  assert IsEmpty({<p,N> : p in monsterPrimes, N in JZeroRankZero | not exists{d : d in Divisors(p*N) | d in  allInds} });
 {d : d in Divisors(N), N in JZeroPositiveRank};

 
  ////////////////////////////////////////////////////////////////////////////////
  // Only need to check maximal (with respect to divisibility) N for rank zero,
  // and minimal N for positive rank
  ////////////////////////////////////////////////////////////////////////////////

  JZeroRankZeroMinimal := [N : N in JZeroRankZero | not exists{M : M in JZeroRankZero | M gt N and IsDivisibleBy(M,N)}];
  JOneRankZeroMinimal  := [N : N in JOneRankZero  | not exists{M : M in JOneRankZero  | M gt N and IsDivisibleBy(M,N)}]; 

  // this next, in progress
  JZeroRankOneMinimal := [N : N in JZeroRankOne | not exists{M : M in JZeroRankOne | M lt N and IsDivisibleBy(N,M)}];
  
 
  ////////////////////////////////////////////////////////////////////////////////
  // Helper functions
  ////////////////////////////////////////////////////////////////////////////////
  
  mu := function(N : gamma := 0)
    return (gamma eq 1 select EulerPhi(N) else 1 )*N*&*[1 + 1/p : p in PrimeDivisors(N)];
  end function;
  
  sturm := function(k,N : gamma := 0)
    return Ceiling(k*mu(N : gamma := gamma)/12);
  end function;
  
  windingQuotientDimension := function(N : gamma := 0);
    M := gamma eq 1 select ModularSymbols(Gamma1(N)) else ModularSymbols(N);
      S := CuspidalSubspace(M);
      dimS := Dimension(S)/2;
      pi := ProjectionMap(S);     
    base := M!<1,[Cusps()|Infinity(),0]>;
    V,m := VectorSpace(M);
    minv := Inverse(m);
    B := NextPrime(sturm(2,N  : gamma := gamma));
      gens := [base];
      G := sub<V | [minv(base)]>;
      // precompute Hecke operators
      for p in PrimesUpTo(NextPrime(B)) do  T := HeckeOperator(M,p); end for;
    bool_fullSpace := false;
      p := 1;
    while p le B and not bool_fullSpace do
      p := NextPrime(p);
      bool_newGens := true;
      while bool_newGens do
        gens := [ m(v): v in Basis(G)];  
        G := sub<V | [minv(g*HeckeOperator(M,p^i)) : g in gens, i in [0,1] ] >;
        bool_newGens := Dimension(G) gt #gens;
        dimPr := Dimension(sub<V | [minv(pi(g)) : g in gens] >);
        bool_fullSpace := dimPr eq dimS;
      end while;
    end while;
    return dimPr,dimS, bool_fullSpace;
  end function;
  
  
  // Modular symbol computations over finite fields are only implemented for Gamma0
  windingQuotientDimensionWithField := function(N : gamma := 0, char := 0);
    F := char eq 0 select Rationals() else FiniteField(char);
    M := ModularSymbols(N,2,F);
      S := CuspidalSubspace(M);
      dimS := Dimension(S)/2;
      pi := ProjectionMap(S);     
    base := M!<1,[Cusps()|Infinity(),0]>;
    V,m := VectorSpace(M);
    minv := Inverse(m);
    B := NextPrime(sturm(2,N  : gamma := gamma));
      gens := [base];
      G := sub<V | [minv(base)]>;
    bool_fullSpace := false;
      p := 1;
    while p le B and not bool_fullSpace do
      p := NextPrime(p);
      bool_newGens := true;
      while bool_newGens do
        gens := [ m(v): v in Basis(G)];  
        G := sub<V | [minv(g*HeckeOperator(M,p^i)) : g in gens, i in [0,1] ] >;
        bool_newGens := Dimension(G) gt #gens;
        dimPr := Dimension(sub<V | [minv(pi(g)) : g in gens] >);
        bool_fullSpace := dimPr eq dimS;
      end while;
    end while;
    return dimPr,dimS, bool_fullSpace;
  end function;


  ////////////////////////////////////////////////////////////////////////////////
  // Main computation
  ////////////////////////////////////////////////////////////////////////////////
  
   for N in JZeroRankZero cat JZeroPositiveRank do  
      if not N in genus0Gamma0 then
      N, "JZero(N)", windingQuotientDimensionWithField(N : gamma := 0, char := NextPrime(10000));
      else
      N, "JZero(N)", "genus 0";
    end if;
   end for;
  
       
  // This one takes quite awhile     
   for N in JZeroRankZero do  
      if not N in genus0Gamma1 then
      N, "JOne(N)", windingQuotientDimension(N : gamma := 1);
      else
      N, "JOne(N)", "genus 0";
    end if;
   end for;
  
  
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // This computes the N for which J_1(2,2N) has rank 0
  //////////////////////////////////////////////////////////////////////////////////////////     
  
    /* Theory: 
        - there is a birational isomorphism X_Delta(4N) --> X1(2,2N)
          - and thus a surjection X_1(4N) --> X1(2,2N). 
          - so if X_1(4N) has rank 0, then so does X1(2,2N)
        - on the other hand, there are maps
          - X1(2,2N) --> X0(4N)
          - X1(2,2N) --> X1(2N)
          -- so if either X0(4N) or X1(2N) have positive rank, then so does X1(2,2N)
        - this resolves everything except for N = 20,26,36. These we do by hand.            
    */
        
    // First eliminate N such that X_0(4N) or X1(2N) have positive rank
    rankMaybeZero:= [N : N in [1..90] | 4*N in JZeroRankZero and not 2*N in JOnePositiveRank];
    
    // These definitely have rank zero
    rankZeroEasy := [N : N in rankMaybeZero | not 4*N in JOnePositiveRank];
    
    rankPossiblyPositive := [N : N in rankMaybeZero | not N in rankZeroEasy];
    rankPossiblyPositive;
    
    for N in rankPossiblyPositive do
      delta := [i*(j*2*N + 1) : i in [1,-1], j in [0,1]];
      delta;
      dec := Decomposition(JH(4*N,delta));
      for A in dec do
        Dimension(A), IsZeroAt(LSeries(A),1);
      end for;
    end for;
  
    // Delta is not cyclic in this example, so we can't 
    dec20 := Decomposition(JOne(4*20));
    for A in dec20 do
        Dimension(A), IsZeroAt(LSeries(A),1);
    end for;
    // dec20[18] is the only positive rank one.
    dec20_2 := Decomposition(JH(80,16)); // surjects onto X_delta
    dec20[18] in dec20_2;                // false
    
    // additional check: the characters attached to factors of X_1(2,2N) have order 2.
    Order(DirichletCharacter(dec20[18])); // 4
            
    dec26 := Decomposition(JH(4*26,4));
    for A in dec26 do
        IsZeroAt(LSeries(A),1);
    end for;
        
    dec36 := Decomposition(JH(4*36,4));
    for A in dec36 do
        IsZeroAt(LSeries(A),1);
    end for;
  
     
  //////////////////////////////////////////////////////////////////////////////////////////
  // From Maarten Derickx; same idea, via modular symbols.
  // We use this to double check our calculations
  //////////////////////////////////////////////////////////////////////////////////////////
  
    function TwistOfSimpleModularSymbolsSpace(S,chi);
  //On input of a irreducible new modular symbols space S and a character chi of prime power 
  //modulus output the modular symbols space corresponding to the twist of M by the primitive character associated to chi
  //i.e. if f is the modular form associated to S then outputs the modular symbols space
  //corresponding to the modular form f_chi
  //the sign of M should be 1 or -1, and the sign of the output will be the same
    if Conductor(chi) eq 1 then;
      return S;
    end if;
    chi := AssociatedPrimitiveCharacter(chi);
    chi_S := DirichletCharacter(S);
    m := Modulus(chi_S);
    n := Modulus(chi);
    assert IsPrimePower(n);
    p := PrimeDivisors(n)[1];
    sign := Sign(S);
    k := Weight(S);
    chi_t := Extend(chi_S*chi^2,m*n^2);
    Mt := ModularSymbols(chi_t,k,sign);
    St := CuspidalSubspace(Mt);
    for Si in NewformDecomposition(St) do;
      Snew := AssociatedNewSpace(Si);
      tf,chi_i := IsTwist(S,Snew,p);
      if tf and AssociatedPrimitiveCharacter(chi_i) eq chi then;
        return Snew;
      end if;
    end for;
    print "Did not find a twist while we should have!!!!";
    assert false;
  end function;
    
  function PostiveRankNewFactors(m,n);
    //Let m,n be two integer and let 
    //G be the congruence subgroup given by the matrices of the form
    //
    //   [a b]
    //   [c d]
    //
    //with a,d congruent to 1 modulo mn
    //and c congruent to 0 modulo m^2n
    //Then this function returns one modular symbols space
    //for every isogeny class of simple abelian varieties that occurs as an isogeny factor
    //of J(G) and obtains positive rank over Q(zeta_m)
    pr_new_factors := [];
    D := FullDirichletGroup(m*n);
    Chi := FullDirichletGroup(m);
    for chi in Elements(Chi) do;
      for d in GaloisConjugacyRepresentatives(D) do;
        d1 := Extend(d,m^2*n);
        //if IsOdd(d) then;
        //  continue;
        //end if;
  
        M := ModularSymbols(d1,2,1);
        S := CuspidalSubspace(M);
        for Si in NewformDecomposition(S) do;
          Snew := AssociatedNewSpace(Si);
          St := TwistOfSimpleModularSymbolsSpace(Snew,chi);
          if Dimension(St) ne Dimension(WindingSubmodule(St)) then;
            Append(~pr_new_factors,Snew);
          end if;
        end for;
      end for;
    end for;
    return pr_new_factors;
  end function;
  
  function IsX1mnRankZero(m,n);
    return #PostiveRankNewFactors(m,n) eq 0;
  end function;

   
  // double check; takes about 200 seconds
  assert {N : N in JZeroRankZero} diff JOnePositiveRank eq {N : N in   {N : N in JZeroRankZero} diff JOnePositiveRank|  IsX1mnRankZero(1,N)};
  assert JOnePositiveRank eq {N : N in JOnePositiveRank| not IsX1mnRankZero(1,N)};
  assert {N : N in {20,26,36} | not IsX1mnRankZero(2,N)} eq {36};
              
  // This part is in progress, maybe delete     
   for N in JOneRankZeroMinimal do  
      if not N in genus0Gamma1 then
      N, "JOne(N)", windingQuotientDimension(N : gamma := 1);
      else
      N, "JOne(N)", "genus 0";
    end if;
   end for;
