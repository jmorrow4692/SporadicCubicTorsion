//////////////////////////////////////////////////////////////////////////////////////////
// This is a complete determination of the cubic points on X1(22)
//////////////////////////////////////////////////////////////////////////////////////////
 
/****************************************************************************** 
 Here is a summary of the argument.
  - X1(22) has genus 6 and gonality at least 4
  - A priori, we know that there are 10 rational cusps
  - Over F_3, we find:
	- 10 places of degree 1 on X1(22)
	- 0  places of degree 2 on X1(22)
	- 0  places of degree 3 on X1(22)
  - Since we have 10 rational points already, we are done
******************************************************************************/

  //////////////////////////////////////////////////////////////////////////////////////////
  // from Drew's page
  //////////////////////////////////////////////////////////////////////////////////////////
  
  N := 22;
    F := FiniteField(3);           
  A2<x,y> := AffineSpace(F,2);
    X:=Curve(A2,y^4 + (x^3 + 2*x^2 + x + 2)*y^3 + (x^5 + x^4 + 2*x^3 + 2*x^2 + 1)*y^2 + (x^5 - x^4 - 2*x^3 - x^2 - x)*y - x^4 - x^3);

    
  #Places(X,1); // 10
  #Places(X,2); // 0
  #Places(X,3); // 0
