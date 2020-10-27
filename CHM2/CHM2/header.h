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
	vector<mytype> X;
	vector<mytype> X1;
	mytype pogr = 1e-14;
	mytype nevas;
	mytype pogrx;
	mytype va;
	mytype w;
	mytype normaB;
	int N, m;

public:
	void Readfile();
	void Jacoby();
	void Gaus_Zeidel();
	void Printfile(int n);
	mytype otn_nevas();
	mytype otn_pogr();
	mytype norma(vector <mytype> f);
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

template<typename mytype>
void slau<mytype>::Gaus_Zeidel()
{
	Readfile();
	mytype i, temp = 0, temp1 = 0, k, ks = 4, ke, di, li, x0 = 0, x1 = N;
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
}

template<typename mytype>
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
}

template<typename mytype>
mytype slau<mytype>::otn_nevas() 
{
	mytype ks = 3, ke = 9, temp = 0, li;
	vector<mytype> f;
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
		f.push_back(temp);
		temp = 0;
		if (ks > 0) ks--;
	}

	mytype fb = 0;
	for (int i = 0; i < N; i++) 
	{
		f[i] = B[i] - f[i];
		fb += f[i] * f[i];
	}
	fb = norma(f);
	//b = norma(B);
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
mytype slau<mytype>::norma(vector<mytype> f)
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