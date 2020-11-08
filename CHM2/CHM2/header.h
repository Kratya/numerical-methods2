#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

template<class mytype>
class slau
{
private:
	vector<mytype> al;
	vector<mytype> ai;
	vector<mytype> B;
	mytype pogr = 1e-14;
	mytype nevas;
	mytype pogrx;
	mytype va;
	mytype w;
	mytype normaB;
	int N, m;

public:
	vector<mytype> X;
	vector<mytype> X1;
	void Readfile();
	//void Jacoby();
	//void Gaus_Zeidel();
	void Printfile(int n);
	mytype otn_nevas();
	mytype otn_nevas(mytype** XTN, int bsize);
	mytype otn_pogr();
	mytype norma(vector <mytype> &f);
	void JSMethod(vector<mytype>& X_0, vector<mytype> &X_1);
	void BMethod(vector<mytype>& X_0, int bsize, mytype w);
};

template<typename mytype>
void slau<mytype>::Readfile() 
{
	int i;
	mytype temp;
	ifstream filein;
	filein.open("in.txt");
	filein >> N;
	filein >> m;
	filein >> nevas;
	filein >> w;

	ai.resize(9);
	al.resize(N * 9);
	B.resize(N);
	X.resize(N);
	X1.resize(N);

	ai[0] = 0 - (m + 1) - 3;
	ai[1] = 0 - (m + 1) - 2;
	ai[2] = 0 - (m + 1) - 1;
	ai[3] = - 1;
	ai[4] = 0;
	ai[5] = 1;
	ai[6] = 1 + m + 1;
	ai[7] = 1 + m + 2;
	ai[8] = 1 + m + 3;

	for (i = 0; i < 9; i++)
	{
		for (int j = 0; j < N; j++) 
		{
			filein >> temp;
			al[i * N + j] = temp;
		}
	}

	for (i = 0; i < N; i++)
	{
		filein >> temp;
		B[i] = (temp);
	}

	for (i = 0; i < N; i++)
	{
		X[i] = 0;
		X1[i] = 0;
	}

	normaB = norma(B);
}

/*template<typename mytype>
void slau<mytype>::Gaus_Zeidel()
{
	Readfile();
	mytype i, temp = 0, k, ks = 4, ke, di, li, x0 = 0, x1 = N;
	int flag = 1;
	di = ks * N;
	cout.precision(15);
	while (nevas > pogr && flag < 30000)
	{
		ks = 4;
		ke = 9;
		for (i = 0; i < N; i++) 
		{
			temp = B[i];
			for (int j = ks; j < ke; j++)
			{
				k = j * N + i;
				li = i + ai[j];
				if (li >= 0) 
				{
					if (li >= N) ke = j;
					else
						temp -= al[k] * X[li];
				}
			}
			X[i] = X[i] + w * temp / al[di + i];
			if (ks > 0) ks--;
		}

		nevas = otn_nevas();
		cout << flag << endl;
		for (int j = 0; j < N; j++)
			cout << X[j] << " ";
		cout << endl;
		flag++;
	}
	Printfile(flag);
	cin.get();
}*/

/*template<typename mytype>
void slau<mytype>::Jacoby() 
{
	Readfile();
	mytype i, temp = 0, k, ks = 4, ke, di, li, flag = 1, x0 = 0, x1 = N;
	di = ks * N;
	cout.precision(15);
	while (nevas > pogr && flag < 30000) 
	{
		ks = 4;
		ke = 9;
		for (i = 0; i < N; i++) 
		{
			temp = B[i];
			for (int j = ks; j < ke; j++) 
			{
				k = j * N + i;
				li = i + ai[j];
				if (li >= 0)
				{
					if (li >= N) ke = j;
					else
						temp -= al[k] * X[li];
				}
			}
			X1[i] = X[i] + w * temp / al[di + i];
			if (ks > 0) ks--;
		}

		X = X1;
		nevas = otn_nevas();
		cout << flag << endl;
		for (int j = 0; j < N; j++)
			cout << X[j] << " ";
		cout << endl;
		flag++;
	}
	Printfile(flag);
	cin.get();
}*/
/*
template<typename mytype>
void slau<mytype>::JSMethod(vector<mytype>& X)
{
	Readfile();
	mytype i, temp = 0, k, ks = 4, ke, di, li, flag = 1, x0 = 0, x1 = N;
	di = ks * N;
	cout.precision(15);
	while (nevas > pogr && flag < 30000)
	{
		ks = 4;
		ke = 9;
		for (i = 0; i < N; i++)
		{
			temp = B[i];
			for (int j = ks; j < ke; j++)
			{
				k = j * N + i;
				li = i + ai[j];
				if (li >= 0)
				{
					if (li >= N) ke = j;
					else
						temp -= al[k] * X[li];
				}
			}
			X[i] = X[i] + w * temp / al[di + i];
			if (ks > 0) ks--;
		}

		nevas = otn_nevas();
		cout << flag << endl;
		for (int j = 0; j < N; j++)
			cout << X[j] << " ";
		cout << endl;
		flag++;
	}
	Printfile(flag);
	cin.get();
}
*/
template<typename mytype>
void slau<mytype>::JSMethod(vector<mytype>& X_0, vector<mytype>& X_1)
{
	Readfile();
	mytype temp = 0;
	int i, k, ks = 4, ke, di, li, flag = 1;
	di = ks * N;
	cout.precision(15);
	while (nevas > pogr && flag < 30000)
	{
		ks = 4;
		ke = 9;
		for (i = 0; i < N; i++)
		{
			temp = B[i];
			for (int j = ks; j < ke; j++)
			{
				k = j * N + i;
				li = i + ai[j];
				if (li >= 0)
				{
					if (li >= N) ke = j;
					else
						temp -= al[k] * X_0[li];
				}
			}
			X_1[i] = X_0[i] + w * temp / al[di + i];
			if (ks > 0) ks--;
		}

		X_0 = X_1;

		nevas = otn_nevas();
		cout << flag << endl;
		for (int j = 0; j < N; j++)
			cout << X_0[j] << " ";
		cout << endl;
		flag++;
	}
	Printfile(flag);
	cin.get();
}
// БЛОЧНАЯ РЕЛАКСАЦИЯ 
template<typename mytype>
void slau<mytype>::BMethod(vector<mytype>& X_0, int bsize, mytype w)
{
	Readfile();
	cout << "N = " << N << endl;
	if (bsize > N)
	{
		cout << "err1. Размер блока больше длины диагонали" << endl;
		return;
	}
	if (N % bsize != 0)
	{
		cout << "err2. Невозможно разделить матрицу на равные блоки" << endl;
		return;
	}
	int i, k, li, ks, ke, numb, ibord, jbord;
	numb = N / bsize;
	mytype** tmat = new mytype*[N];

	for (int count = 0; count < N; count++)
		tmat[count] = new mytype[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			tmat[i][j] = 0;

	i = 0;

	ks = 4;
	ke = 9;
	for (i = 0; i < N; i++)
	{
		for (int j = ks; j < ke; j++)
		{
			k = j * N + i;
			li = i + ai[j];
			if (li >= 0)
			{
				if (li >= N) ke = j;
				else
					tmat[i][li] = al[k];
			}
		}
		if (ks > 0) ks--;
	}
	
	mytype**** blocks = new mytype***[numb];
	for (int count = 0; count < numb; count++)
	{
		blocks[count] = new mytype**[numb];
		for (int c = 0; c < numb; c++)
		{
			blocks[count][c] = new mytype *[bsize];
			for (int ci = 0; ci < bsize; ci++)
				blocks[count][c][ci] = new mytype[bsize];
		}
	}

	ibord = 0;
	jbord = 0;
	for (int bi = 0; bi < numb; bi++)
	{
		for (int bj = 0; bj < numb; bj++)
		{
			for (int ci = 0; ci < bsize; ci++)
			{
				for (int cj = 0; cj < bsize; cj++)
					blocks[bi][bj][ci][cj] = tmat[ci + ibord][cj + jbord];
			}
			jbord += bsize;
		}
		jbord = 0;
		ibord += bsize;
	}

	for (int count = 0; count < N; count++)
		delete[] tmat[count];
	delete[] tmat;

	ibord = 0;
	mytype** FT = new mytype * [numb];
	for (int count = 0; count < numb; count++)
	{
		FT[count] = new mytype[bsize];
		for (int ci = 0; ci < bsize; ci++)
		{
			FT[count][ci] = B[ci + ibord];
		}
		ibord += bsize;
	}

	mytype** XT = new mytype * [numb];
	mytype** XTN = new mytype * [numb];
	mytype** R = new mytype * [numb];
	mytype*** L = new mytype ** [numb];
	mytype*** U = new mytype ** [numb];
	mytype** Z = new mytype * [numb];
	mytype** Y = new mytype * [numb];

	for (int count = 0; count < numb; count++)
	{
		L[count] = new mytype * [bsize];
		U[count] = new mytype * [bsize];
		for (int ci = 0; ci < bsize; ci++)
		{
			L[count][ci] = new mytype[bsize];
			U[count][ci] = new mytype[bsize];
			for (int cj = 0; cj < bsize; cj++)
			{
				L[count][ci][cj] = 0;
				U[count][ci][cj] = 0;
			}
		}
	}

	for (int count = 0; count < numb; count++)
	{
		XT[count] = new mytype[bsize];
		XTN[count] = new mytype[bsize];
		R[count] = new mytype[bsize];
		Z[count] = new mytype[bsize];
		Y[count] = new mytype[bsize];
		for (int ci = 0; ci < bsize; ci++)
		{
			XT[count][ci] = 0;
			XTN[count][ci] = 0;
			R[count][ci] = 0;
			Z[count][ci] = 0;
			Y[count][ci] = 0;
			
		}
	}

	int flag = 1;
	do
	{
		for (int icount = 0; icount < numb; icount++)
		{
			for (int ci = 0; ci < bsize; ci++)
			{
				R[icount][ci] = FT[icount][ci];
				for (int jcount = 0; jcount < icount; jcount++)
					for (int cj = 0; cj < bsize; cj++)
						R[icount][ci] -= blocks[icount][jcount][ci][cj] * XTN[jcount][cj];
				for (int jcount = icount; jcount < numb; jcount++)
					for (int cj = 0; cj < bsize; cj++)
						R[icount][ci] -= blocks[icount][jcount][ci][cj] * XT[jcount][cj];
				R[icount][ci] *= w;
			}
		}
		
		mytype sum = 0;
		for (int icount = 0; icount < numb; icount++)
		{
			for (int k = 0; k < bsize; k++)
			{
				U[icount][k][k] = 1;
				for (int ci = k; ci < bsize; ci++)
				{
					sum = 0;
					for (int p = 0; p < k - 1; p++)
						sum += L[icount][ci][p] * U[icount][p][k];
					L[icount][ci][k] = blocks[icount][icount][ci][k] - sum;
				}

				for (int cj = k + 1; cj < bsize; cj++)
				{
					sum = 0;
					for (int p = 0; p < k - 1; p++)
						sum += L[icount][k][p] * U[icount][p][cj];
					U[icount][k][cj] = (blocks[icount][icount][k][cj] - sum) / L[icount][k][k];
				}
			}
		}

		for (int icount = 0; icount < numb; icount++)
		{
			for (int ci = 0; ci < bsize; ci++) 
			{
				sum = 0;
				for (int p = 0; p < ci; p++)
					sum += L[icount][ci][p] * Z[icount][p];
				Z[icount][ci] = (R[icount][ci] - sum) / L[icount][ci][ci];
			}
		}

		for (int icount = 0; icount < numb; icount++)
		{
			for (int ci = bsize - 1; ci >= 0; ci--)
			{
				sum = 0;
				for (int p = bsize - 1; p > ci; p--) 
					sum += U[icount][ci][p] * Y[icount][p];
				Y[icount][ci] = (Z[icount][ci] - sum) / U[icount][ci][ci];
			}
		}
		
		for (int icount = 0; icount < numb; icount++)
		{
			for (int ci = 0; ci < bsize; ci++)
			{
				XT[icount][ci] = XTN[icount][ci];
				XTN[icount][ci] += Y[icount][ci];
			}
		}

		nevas = otn_nevas(XTN, bsize);
		flag++;
	} while (nevas > pogr && flag < 30000);

	ibord = 0;
	for (int count = 0; count < numb; count++)
	{
		for (int ci = 0; ci < bsize; ci++)
		{
			X[ci + ibord] = XTN[count][ci];
		}
		ibord += bsize;
	}

	cout.precision(15);
	for (int count = 0; count < numb; count++)
	{
		for (int ci = 0; ci < bsize; ci++)
		{
			cout << "XT[" << count << "][" << ci << "] = " << XT[count][ci] << endl;
		}
		ibord += bsize;
	}
	cout << endl;

	for (int ci = 0; ci < N; ci++)
		cout << "X[" << ci << "] = " << X[ci] << endl;
	cout << "flag = " << flag << endl;
	cout << "nevas = " << nevas << endl;

	Printfile(flag);

	for (int icount = 0; icount < numb; icount++)
	{
		for (int jcount = 0; jcount < numb; jcount++)
		{
			for (int ci = 0; ci < bsize; ci++)
				delete[] blocks[icount][jcount][ci];
			delete[] blocks[icount][jcount];
		}
		delete[] blocks[icount];
	}
	delete[] blocks;

	for (int icount = 0; icount < numb; icount++)
	{
		for (int ci = 0; ci < bsize; ci++)
		{
			delete[] L[icount][ci];
			delete[] U[icount][ci];
		}
		delete[] L[icount];
		delete[] U[icount];
	}
	delete[] L;
	delete[] U;

	for (int icount = 0; icount < numb; icount++)
	{
		delete[] FT[icount];
		delete[] XT[icount];
		delete[] Y[icount];
		delete[] XTN[icount];
		delete[] R[icount];
		delete[] Z[icount];
	}

	delete[] FT;
	delete[] XT;
	delete[] Y;
	delete[] XTN;
	delete[] R;
	delete[] Z;
}

template<typename mytype>
mytype slau<mytype>::otn_nevas() 
{
	mytype temp = 0;
	int ks = 4, ke = 9, li;
	vector<mytype> f;
	f.resize(N);
	for (int i = 0; i < N; i++) 
	{
		for (int j = ks; j < ke; j++) 
		{
			int k = j * N + i;
			li = i + ai[j];
			if (li >= 0)
			{
				if (li >= N) ke = j;
				else
					temp += al[k] * X[li];
			}
		}
		f[i] = temp;
		temp = 0;
		if (ks > 0) ks--;
	}

	mytype fb = 0;
	for (int i = 0; i < N; i++) 
	{
		f[i] = B[i] - f[i];
		fb += f[i] * f[i];
	}
	//fb = norma(f);
	//b = norma(B);
	return fb / normaB;
}

template<typename mytype>
mytype slau<mytype>::otn_nevas(mytype** XTN, int bsize)
{
	mytype temp = 0;
	int ks = 0, ke = 9,  li;
	vector<mytype> f;
	f.resize(N);
	for (int i = 0; i < N; i++)
	{
		for (int j = ks; j < ke; j++)
		{
			int k = j * N + i;
			li = i + ai[j];
			if (li >= 0)
			{
				if (li >= N) ke = j;
				else
					temp += al[k] * XTN[li / bsize][li % bsize];
			}
		}
		f[i] = temp;
		temp = 0;
		if (ks > 0) ks--;
	}

	mytype fb = 0;
	for (int i = 0; i < N; i++)
	{
		f[i] = B[i] - f[i];
		fb += f[i] * f[i];
	}
	return fb / normaB;
}

template<typename mytype>
mytype slau<mytype>::otn_pogr()
{
	mytype xt, x;
	mytype al = 0;
	for (int i = 0; i < N; i++)
	{
		X1[i] = i + 1;
		al += X1[i] * X1[i];
	}
	xt = sqrt(al);

	al = 0;
	for (int i = 0; i < N; i++) 
	{ 
		X1[i] = X[i] - X1[i];
		al += X1[i] * X1[i];
	}
	x = sqrt(al);
	return x / xt;
}

template<typename mytype>
mytype slau<mytype>::norma(vector<mytype> &f)
{
	mytype al = 0;
	for (int i = 0; i < N; i++)
	{
		al += f[i] * f[i];
	}
	return sqrt(al);
}

template<typename mytype>
void slau<mytype>::Printfile(int n)
{
	ofstream file;
	file.open("out.txt");
	file.precision(15);
	pogrx = otn_pogr();
	va = pogrx / nevas;

	file << n << " " << va << endl;
	for (int i = 0; i < N; i++) file << X[i] << endl;
	file.close();
}