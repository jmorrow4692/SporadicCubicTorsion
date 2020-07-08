//////////////////////////////////////////////////////////////////////////////////////////
// A few homebrewed functions
//////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////
  // This pins down the global class group, given ClassGroup data at several primes
  //////////////////////////////////////////////////////////////////////////////////////////
  
  torsBound := function(data)
    n := Minimum([#d : d in data]);
    return [ GCD( [ d[#d-i]: d in data])  :  i in [1..n-1] ];
  end function;  

  ////////////////////////////////////////////////////////////////////
  // Tuples of Degree n
  //
  // This returns all lists of length m 
  // with positive entries
  // that sum to n
  ////////////////////////////////////////////////////////////////////  

  function tuplesOfDegree(m,n)
    if m eq 0 then 
      return {@@};
    elif n eq 0 then 
      return {@[0 : i in [1..m]]@};
    elif n eq 1 then 
      return {@[i eq j select 1 else 0 : i in [1..m]] : j in [1..m]@};
    elif m eq 1  then 
      return {@[n]@};    
    else 
      return &join{@
             {@ tup cat [n-i]: tup in $$(m-1,i)@}
	     : i in [0..n]
	     @};
    end if;
  end function;

  
  /////////////////////////////////////////////////////////////////////////
  // divisorSearch
  //
  // Input: a smooth projective curve C
  //        a positive integer deg
  //        a positive integer bound
  //        an optional boolean verbose
  // 
  // Output: the defining ideals of any points of degree deg on C
  //        which arise as the intersection of L and C, 
  //        where L is a hyperplane with coefficients of size at most bound
  //
  ///////////////////////////////////////////////////////////////////////// 

  function divisorSearch(C,deg,bound: verbose := false)
    P<[X]> := AmbientSpace(C);
    B := -1;
    tups := {};
    d := Dimension(P);
    ideals := {@@};
      while B lt bound do
	  tupsOld := tups;
          B := B + 1;  
  	if verbose then B; end if;
          tups := [P![tup[i] : i in [1..#tup]] : 
                   tup in CartesianPower([-B..B],d+1) |
                   not tup eq <0 : i in [1..d+1]>
                   and not P![tup[i] : i in [1..#tup]] in tupsOld];
          tups := [P![tup[i] : i in [1..#tup]] : 
                   tup in CartesianPower([-B..B],d+1) |
                   not tup eq <0 : i in [1..d+1]>
                   and not P![tup[i] : i in [1..#tup]] in tupsOld];
          for c in tups  do  
            L := Scheme(P,&+[c[i]*X[i] : i in [1..Dimension(AmbientSpace(C)) + 1]]);
              irr := IrreducibleComponents(L meet C);
              degs := [Degree(ReducedSubscheme(i)) : i in irr];
              if  deg in degs then
    	      if verbose then c; end if;
    	    for cpt in irr do  
    	      if  Degree(ReducedSubscheme(cpt)) eq deg 
    	          and not Basis(Ideal(ReducedSubscheme(cpt))) in ideals then	      
    	            ideals := ideals join {@Basis(Ideal(ReducedSubscheme(cpt)))@};	   
  		    if verbose then ideals; end if;
    	      end if;
    	    end for;
              end if;  
          end for;
       end while;
     return ideals;
  end function;  
 

  
  /////////////////////////////////////////////////////////////////////////
  // hasTwistWithNTorsion
  //
  // Input: an Elliptic Curve E
  //        a positive integer N
  //        a positive integer bound
  // 
  // Output: bool, p
  //        bool is false if there is a prime p less than bound such that for all primes
  //        pp over p, E mod pp has no quadratic twist with an N-torsion point
  //
  //        In this case, no quadratic twist of E has an N torsion point
  /////////////////////////////////////////////////////////////////////////

  function hasTwistWithNTorsion(E,N, bound)
    E := IntegralModel(E);
    K := BaseRing(E);
    OK := Integers(K);
    badPrimes := PrimeDivisors(Norm(Conductor(E))) cat
		 PrimeDivisors(Discriminant(OK)) cat
		 PrimeDivisors(N);
    bool := true;
    p := 1;
    while bool do
      p := NextPrime(p);
      if not p in badPrimes then
        bool := exists{pp : pp in Factorization(p*OK) |
                0 in [o mod N : o in
  	          Invariants(TorsionSubgroup(e)) cat
	          Invariants(TorsionSubgroup(QuadraticTwist(e)))]
                  where e is Reduction(Emin,pp[1])
                  where _,Emin is LocalInformation(E,pp[1]) 				      
	    };
      end if;
    end while;
    return bool, p;	
  end function;

  
  //////////////////////////////////////////////////////////////////////////////////////////
  // The Jacobian of a general curve isn't implemented in Magma.
  //////////////////////////////////////////////////////////////////////////////////////////
  // This function computes the order of a divisor D.
  // Optional paramater: N is a known multiple of the order of D (and we assume that N*D = 0)
  // Optional paramater: B, give up after B
  //////////////////////////////////////////////////////////////////////////////////////////
  // Note: usually it is much faster to compute this mod a prime, but often it is convenient
  // to compute this globally.
  //////////////////////////////////////////////////////////////////////////////////////////  

  order := function(D : N := 0, B := N eq 0 select 500 else N)
    currentBound := 0;
    found := false;
    if N eq 0 then
      while not found and currentBound le B do
        currentBound := currentBound + 1;
        found := IsPrincipal(currentBound * D);
      end while;
        return found select currentBound else 0, found;
    else
      for n in [d : d in Divisors(N) | not d eq N] do
        found := IsPrincipal(n*(D));
        if found then
          currentBound := n; break;
        end if;
      end for;
      return found select currentBound else N, true;
     end if;
  end function;

  //testD  := Divisor(Xsm![-1 , -1 , 1 , 1 , 1]) - Divisor(Xsm![-2 , -1 , 1 , 1 , 0]);
  //order(testD : N := 728);


  //////////////////////////////////////////////////////////////////////////////////////////
  // twoTorsionRank
  //////////////////////////////////////////////////////////////////////////////////////////
  // This function computes the 2-rank of the torsion of Jacobian of a hyperelliptic curve
  //////////////////////////////////////////////////////////////////////////////////////////  
  // INPUT: an even degree polynomial
  // OUTPUT: the 2-rank of J(Q)_tors, where J is the Jacobian of the hyperelliptic curve
  //         y^2 = f
  //////////////////////////////////////////////////////////////////////////////////////////

  twoTorsionRank := function(f)
    assert IsSquarefree(f);
    G := GaloisGroup(f);
    n := Degree(G);
    V := VectorSpace(GF(2), n);
    A := Matrix([V.i-V.n : i in [1..n-1]]);
    At := Transpose(A);
    ann := [A * (PermutationMatrix(GF(2), g)-1) * At : g in Generators(G)];
    ann := BlockMatrix(1, #ann, ann);
    return Dimension(Kernel(ann))-1;
  end function;


  /////////////////////////////////////////////////////////////////////////
  // realComponentsX1
  //
  // Input: an integer N
  // 
  // Output: the number of real connected components of X_1(N)
  //
  // (Formula from Real Components Of Modular Curves by Andrew Snowden)
  ///////////////////////////////////////////////////////////////////////// 

  function realComponentsX1(N);
    r := Valuation(N,2);
    M := Integers()!(N/2^r);

    u,m := UnitGroup(ResidueClassRing(M));
    uSub := sub<u | [Inverse(m)(-1), Inverse(m)(2)]>;

    return N eq 4 select 1 
      else (not N eq 4) and r ge 2 select Integers()!((1/4)*EulerPhi(N))
      else Integers()!(#u/#uSub);
  end function;


