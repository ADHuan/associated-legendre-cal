// legendre.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>

#define PI acos(-1)

using namespace std;

void belikov(int degree, double theta) {
	stringstream fn;
	fn << "./legendre_belikov_" << degree << '_' << theta << ".txt";


	degree++;
	long double** pp = new long double* [degree];
	long double** p = new long double* [degree];//ans
	theta = theta * PI / 180;
	long double u = sin(theta);
	long double t = cos(theta);

	pp[0] = new long double[1]{ 1.0 };
	//pp[0][1]-->0
	pp[1] = new long double[2]{ t , u };


	for (int i = 2; i < degree; i++) {
		pp[i] = new long double[i + 1];
		for (int j = 0; j <= i; j++) {
			if (j == 0) {
				pp[i][0] = t * pp[i - 1][0] - u / 2.0 * pp[i - 1][1];
			}
			else { 
				long double a1 = pp[i - 1][j];
				long double a2 = pp[i - 1][j + 1];
				long double a3 = pp[i - 1][j - 1];
				if (j == i - 1)
					a2 = 0.0;
				if (j == i) {
					a1 = 0.0;
					a2 = 0.0;
				}
					
				pp[i][j] = t * a1 - u * (1.0 / 4.0 * a2 - a3);
			}
		}
	}

	long double** nnm = new long double* [degree];
	nnm[0] = new long double[1]{ 1.0 };
	nnm[1] = new long double[2]{ 1.0,1.0 };

	for (int i = 2; i < degree; i++) {
		nnm[i] = new long double[i + 1];
		for (int j = 0; j <= i; j++) {
			if (j >= 0 && j <= i - 1) {
				nnm[i][j] = sqrt(1.0 -(double(j) * j) / (double(i) * i)) * nnm[i - 1][j];
			}
			else {
				nnm[i][j] = sqrt(1.0 - (1.0 / (2.0 * i))) * nnm[i - 1][j - 1];
			}		
		}
	}

	//p[0] = new long double[1]{ 1 };
	//p[1] = new long double[2]{ sqrt(3) * t ,  sqrt(3) * u };




	ofstream outfile;
	outfile.open(fn.str(), ios::out | ios::trunc);
	outfile << "   n  " << "   m  " << "Pnm" << endl;

	for (int i = 0; i < degree; i++) {
		p[i] = new long double[i + 1];
		for (int j = 0; j <= i; j++) {
			p[i][j] = sqrt(2.0 * i + 1.0) * nnm[i][j] * pp[i][j];
			outfile.setf(ios::dec);
			outfile.width(4);
			outfile << i << "  ";
			outfile.width(4);
			outfile << j << "  ";
			outfile.setf(ios::scientific);
			outfile.width(18);
			outfile.precision(10);
			outfile << p[i][j] << endl;
		}
	}
	outfile.close();
}

void swarztrauber(int degree, double theta) {
	stringstream fn;
	fn << "./legendre_swarztrauber_" << degree << '_' << theta << ".txt";

	degree++;
	long double** p = new long double* [degree];//ans
	theta = theta * PI / 180;
	long double u = sin(theta);
	long double t = cos(theta);

	p[0] = new long double[1]{ 1.0 };
	//pp[0][1]-->0
	p[1] = new long double[2]{ sqrt(3.0)*t , sqrt(3.0) * u };
	long double k = 2;

	ofstream outfile;
	outfile.open(fn.str(), ios::out | ios::trunc);
	outfile << "   n  " << "   m  " << "Pnm" << endl;
	for(int i = 0; i < 2; i++){
		for (int j = 0; j <= i; j++) {
			outfile.setf(ios::dec);
			outfile.width(4);
			outfile << i << "  ";
			outfile.width(4);
			outfile << j << "  ";
			outfile.setf(ios::scientific);
			outfile.width(18);
			outfile.precision(10);
			outfile << p[i][j] << endl;
		}
	}



	for (int i = 2; i < degree; i++) {
		p[i] = new long double[i + 1];
		for (int j = 0; j <= i; j++) {
			if (j == 0 || j == 1) {
				long double anm = sqrt(((2.0 * i - 1.0) * (2.0 * i + 1.0)) / ((i - j) * (i + j)));
				long double bnm = sqrt(((2.0 * i + 1.0) * (i + j - 1.0) * (i - j - 1.0)) / ((i - j) * (i + j) * (2.0 * i - 3.0)));
				p[i][j] = anm * t * p[i - 1][j] - bnm * p[i - 2][j];
			}
			else {
				if (j == 2)
					k = 2.0;
				else
					k = 1.0;
				long double anm = sqrt(((2.0 * i + 1.0) * (i - j) * (i - j - 1.0)) / ((2.0 * i - 3.0) * (i + j) * (i + j - 1.0)));
				long double bnm = sqrt((k * (2.0 * i + 1.0) * (i + j - 2.0) * (i + j - 3.0)) / ((2.0 * i - 3.0) * (i + j) * (i + j - 1.0)));
				long double cnm = sqrt((k * (i - j + 1.0) * (i - j + 2.0)) / ((i + j) * (i + j - 1.0)));
				long double a1 = p[i - 2][j];
				if (j == i - 1 || j == i)
					a1 = 0;
				p[i][j] = anm * a1 + bnm * p[i - 2][j - 2] - cnm * p[i][j - 2];
			}
			outfile.setf(ios::dec);
			outfile.width(4);
			outfile << i << "  ";
			outfile.width(4);
			outfile << j << "  ";
			outfile.setf(ios::scientific);
			outfile.width(18);
			outfile.precision(10);
			outfile << p[i][j] << endl;
		}
	}
	outfile.close();
}

int main()
{
	int degree;
	double theta;	
	int opt;
	cout << "Please input degree:" << endl;
	cin >> degree;
	cout << "Please input theta:" << endl;
	cin >> theta;
	cout << "Choose method number:" << endl << "1.Belikov method" << endl << "2.Swarztrauber method" << endl;
	cin >> opt;

	if (opt == 2)
		swarztrauber(degree, theta);
	else
		belikov(degree, theta);

	cout << "FINISHED!" << endl << "Output:./legendre_{method}_{degree}_{theta}.txt" << endl;
	system("pause");
}

