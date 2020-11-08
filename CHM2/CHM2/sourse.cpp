#include "Header.h"

void main() 
{
	slau<double> A;
	//A.Jacoby();
	//A.Gaus_Zeidel();
	//A.JSMethod(A.X, A.X);
	A.BMethod(A.X, 2, 0.1);
}
