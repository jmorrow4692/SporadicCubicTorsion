//////////////////////////////////////////////////////////////////////////////////////////
// This is a complete determination of the cubic points on X1(25)
//////////////////////////////////////////////////////////////////////////////////////////
 
/****************************************************************************** 
 Here is a summary of the argument.
  - X1(25) has genus 12 and gonality at least 4
  - Over F_3, we find:
	- 10 places of degree 1  on X1(25)
	- 0  places of degree 2  on X1(25)
	- 0  places of degree 13 on X1(25)
  - A priori, we know that there are 10 rational cusps, so we are done
    (We can work with a singluar model, so this is very fast)
******************************************************************************/

  //////////////////////////////////////////////////////////////////////////////////////////
  // from Drew's page
  //////////////////////////////////////////////////////////////////////////////////////////
  
  N := 25;
    F := FiniteField(3);           
  A2<u,v> := AffineSpace(F,2);
    X := Curve(A2,u*v^5 + (u^4 - 2*u^3 - u^2 + 2*u + 1)*v^4 - (2*u^6 - 2*u^4 + 4*u^3 + 2*u^2 - 2)*v^3 + (u^8 + u^7 - 2*u^6 + u^5 - u^4 - u^3 - 2*u^2 - u + 1)*v^2 + (u^8 + u^7 + 2*u^6 + u^5 - 2*u^4 + u^3 - u^2)*v + u^6);
          
  #Places(X,1); // 10
  #Places(X,2); // 0
  #Places(X,3); // 0
