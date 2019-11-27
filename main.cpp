#include<iostream>
#include<fstream>
#include<cstring>
#include<math.h>
using namespace std;
const int N = 100;
double q[5][5] = { 0 }, a[5] = { 0 }, d[5] = { 0 };
void ALU(int n)
{
	for (int i = 0; i < 4; i++)
		a[i] = 0;
	if (n == 3)
	{
		double y[3] = { 0 };
		double L[3][3] = { 0 }, U[3][3] = { 0 }, A[3][3] = { 0 };
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A[i][j] = q[i][j];
		for (int i = 0; i < n; i++)		  //计算U、L
		{
			U[0][i] = A[0][i];
			L[i][0] = A[i][0] / A[0][0];
			U[1][i] = A[1][i] - L[1][0] * U[0][i];//u2j = a2j-l21*u1j
			L[i][1] = (A[i][1] - L[i][0] * U[0][1]) / U[1][1]; // li2= (ai2 - li1*u12)/u22
			U[2][i] = A[2][i] - L[2][0] * U[0][i] - L[2][1] * U[1][i]; //u3j = a3j-l31*u1j-l32*u2j
			L[i][2] = (A[i][2] - L[i][0] * U[0][2] - L[i][1] * U[1][2]) / U[2][2]; //li3 = (ai3-li1*u13-li2*u23)/u33
		}
		//A=LU 
		y[0] = d[0];
		y[1] = d[1] - L[1][0] * y[0];
		y[2] = d[2] - L[2][0] * y[0] - L[2][1] * y[1];
		a[2] = y[2] / U[2][2]; //cout << "a[2] = " << a[2] << endl;
		//cout << y[2] << " " << U[2][2] << endl;//xn = yn/unn  y2=0?
		a[1] = (y[1] - U[1][2] * a[2]) / U[1][1];// cout << "a[1] = " << a[1] << endl;
		a[0] = (y[0] - U[0][1] * a[1] - U[0][2] * a[2]) / U[0][0];//x1=(y1-u12*x2-u13*x3)/u11  //为什么a0会变化

	}
	else if (n == 4)
	{
		double y[4] = { 0 };
		double L[4][4] = { 0 }, U[4][4] = { 0 }, A[4][4] = { 0 };
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A[i][j] = q[i][j];
		for (int i = 0; i < n; i++)		  //计算U、L
		{
			U[0][i] = A[0][i];
			L[i][0] = A[i][0] / A[0][0];
			U[1][i] = A[1][i] - L[1][0] * U[0][i];//u2j = a2j-l21*u1j
			L[i][1] = (A[i][1] - L[i][0] * U[0][1]) / U[1][1]; // li2= (ai2 - li1*u12)/u22
			U[2][i] = A[2][i] - L[2][0] * U[0][i] - L[2][1] * U[1][i]; //u3j = a3j-l31*u1j-l32*u2j
			L[i][2] = (A[i][2] - L[i][0] * U[0][2] - L[i][1] * U[1][2]) / U[2][2]; //li3 = (ai3-li1*u13-li2*u23)/u33
			U[3][i] = A[3][i] - L[3][0] * U[0][i] - L[3][1] * U[1][i] - L[3][2] * U[2][i];
			L[i][3] = (A[i][3] - L[i][0] * U[0][3] - L[i][1] * U[1][3] - L[i][2] * U[2][3]) / U[3][3];
		}
		//A=LU 
		y[0] = d[0];
		y[1] = d[1] - L[1][0] * y[0];
		y[2] = d[2] - L[2][0] * y[0] - L[2][1] * y[1];
		y[3] = d[3] - L[3][0] * y[0] - L[3][1] * y[1] - L[3][2] * y[2];
		a[3] = y[3] / U[3][3];
		a[2] = (y[2] - U[2][3]*a[3]) / U[2][2]; 
		a[1] = (y[1] - U[1][2] * a[2] - U[1][3] *a[3]) / U[1][1];// cout << "a[1] = " << a[1] << endl;
		a[0] = (y[0] - U[0][1] * a[1] - U[0][2] * a[2] - U[0][3] * a[3]) / U[0][0];//x1=(y1-u12*x2-u13*x3)/u11  //为什么a0会变化
	}
	else if (n == 5)
	{
		double y[5] = { 0 };
		double L[5][5] = { 0 }, U[5][5] = { 0 }, A[5][5] = { 0 };
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A[i][j] = q[i][j];
		for (int i = 0; i < n; i++)		  //计算U、L
		{
			U[0][i] = A[0][i];
			L[i][0] = A[i][0] / A[0][0];
			U[1][i] = A[1][i] - L[1][0] * U[0][i];//u2j = a2j-l21*u1j
			L[i][1] = (A[i][1] - L[i][0] * U[0][1]) / U[1][1]; // li2= (ai2 - li1*u12)/u22
			U[2][i] = A[2][i] - L[2][0] * U[0][i] - L[2][1] * U[1][i]; //u3j = a3j-l31*u1j-l32*u2j
			L[i][2] = (A[i][2] - L[i][0] * U[0][2] - L[i][1] * U[1][2]) / U[2][2]; //li3 = (ai3-li1*u13-li2*u23)/u33
			U[3][i] = A[3][i] - L[3][0] * U[0][i] - L[3][1] * U[1][i] - L[3][2] * U[2][i]; 
			L[i][3] = (A[i][3] - L[i][0] * U[0][3] - L[i][1] * U[1][3] - L[i][2] * U[2][3]) / U[3][3];
			U[4][i] = A[4][i] - L[4][0] * U[0][i] - L[4][1] * U[1][i] - L[4][2] * U[2][i] - L[4][3] * U[3][i];
			L[i][4] = (A[i][4] - L[i][0] * U[0][4] - L[i][1] * U[1][4] - L[i][2] * U[2][4] - L[i][3] * U[3][4]) / U[4][4];
		}
		//A=LU 
		y[0] = d[0];
		y[1] = d[1] - L[1][0] * y[0];
		y[2] = d[2] - L[2][0] * y[0] - L[2][1] * y[1];
		y[3] = d[3] - L[3][0] * y[0] - L[3][1] * y[1] - L[3][2] * y[2];
		y[4] = d[4] - L[4][0] * y[0] - L[4][1] * y[1] - L[4][2] * y[2] - L[4][3] * y[3];
		a[4] = y[4] / U[4][4];
		a[3] = (y[3] - U[3][4] * a[4]) / U[3][3];
		a[2] = (y[2] - U[2][3] * a[3]- U[2][4] * a[4]) / U[2][2];
		a[1] = (y[1] - U[1][2] * a[2] - U[1][3] * a[3] - U[1][4] * a[4]) / U[1][1];// cout << "a[1] = " << a[1] << endl;
		a[0] = (y[0] - U[0][1] * a[1] - U[0][2] * a[2] - U[0][3] * a[3] - U[0][4] * a[4]) / U[0][0];//x1=(y1-u12*x2-u13*x3)/u11  //为什么a0会变化

	}
}

double erjie(double *x, double*y, int n, double X)
{
	// y = a0 + a1*x +a2 *x^2
	for (int i = 0; i < 5; i++)
	{
		d[i] = 0;
		for (int j = 0; j < 5; j++)
			q[i][j] = 0;
	}
	q[0][0] = n;
	for (int i = 0; i < n; i++)
	{
		q[0][1] = q[0][1] + *(x + i);
		q[1][1] = q[1][1] + (*(x + i)) * (*(x + i));
		q[1][2] = q[1][2] + (*(x + i)) * (*(x + i))*(*(x + i));
		q[2][2] = q[2][2] + (*(x + i)) * (*(x + i))*(*(x + i)) * (*(x + i));
		d[0] = d[0] + *(y + i);
		d[1] = d[1] + (*(x + i))*(*(y + i));
		d[2] = d[2] + (*(x + i))*(*(x + i))*(*(y + i));
	}
	q[0][2] = q[1][1];
	q[2][0] = q[1][1];
	q[2][1] = q[1][2];
	q[1][0] = q[0][1];
	ALU(3);
	//cout << "拟合多项式为：y = " << a[0] << "+" << a[1] << "x+" << a[2] << "x^2" << endl;
	//cin.get();
	return a[0]+a[1]*X + a[2] *(X*X);
}
double sanjie(double *x, double*y, int n, double X)
{
	//y=a0+a1*x+a2*x^2+a3*x^3
	for (int i = 0; i < 5; i++)
	{
		d[i] = 0;
		for (int j = 0; j < 5; j++)
			q[i][j] = 0;
	}
	q[0][0] = n;
	for (int i = 0; i < n; i++)
	{
		q[0][1] = q[0][1] + *(x + i);
		q[1][1] = q[1][1] + (*(x + i)) * (*(x + i));
		q[1][2] = q[1][2] + (*(x + i)) * (*(x + i))*(*(x + i));//3次方
		q[2][2] = q[2][2] + (*(x + i)) * (*(x + i))*(*(x + i)) * (*(x + i));//4次方
		q[3][3] = q[3][3] + (*(x + i)) * (*(x + i))*(*(x + i))* (*(x + i)) * (*(x + i))*(*(x + i));
		q[3][2] = q[3][2] + (*(x + i)) * (*(x + i))*(*(x + i))*(*(x + i)) * (*(x + i));//5次方
		d[0] = d[0] + *(y + i);
		d[1] = d[1] + (*(x + i))*(*(y + i));
		d[2] = d[2] + (*(x + i))*(*(x + i))*(*(y + i));
		d[3] = d[3] + (*(x + i))*(*(x + i))*(*(x + i))*(*(y + i));
	}
	q[2][0] = q[1][1];
	q[2][1] = q[1][2];
	q[1][0] = q[0][1];
	q[0][3] = q[1][2];
	q[1][3] = q[2][2];
	q[3][1] = q[2][2];
	q[2][3] = q[3][2];
	q[3][0] = q[1][2];
	q[0][2] = q[1][1];
	ALU(4);
	//cout << "y=" << a[0] << "+" << a[1] << "x+" << a[2] << "x^2+" << a[3] << "x^4" << endl;
	return a[0] + a[1] * X + a[2] * (X*X)+a[3]*(X*X*X);
}
double sijie(double *x, double*y, int n, double X)
{
	//y=a0+a1*x+a2*x^2+a3*x^3
	for (int i = 0; i < 5; i++)
	{
		d[i] = 0;
		for (int j = 0; j < 5; j++)
			q[i][j] = 0;
	}
	q[0][0] = n;
	for (int i = 0; i < n; i++)
	{
		q[0][1] = q[0][1] + *(x + i);
		q[1][1] = q[1][1] + (*(x + i)) * (*(x + i));
		q[1][2] = q[1][2] + (*(x + i)) * (*(x + i))*(*(x + i));//3次方
		q[2][2] = q[2][2] + (*(x + i)) * (*(x + i))*(*(x + i)) * (*(x + i));//4次方
		q[3][3] = q[3][3] + (*(x + i)) * (*(x + i))*(*(x + i))* (*(x + i)) * (*(x + i))*(*(x + i));
		q[3][2] = q[3][2] + (*(x + i)) * (*(x + i))*(*(x + i))*(*(x + i)) * (*(x + i));//5次方
		q[4][4] = q[4][4] + (*(x + i)) * (*(x + i))*(*(x + i))*(*(x + i)) * (*(x + i))*(*(x + i)) * (*(x + i))*(*(x + i));//8次方
		q[3][4] = q[3][4] + (*(x + i)) * (*(x + i))*(*(x + i)) * (*(x + i))*(*(x + i)) * (*(x + i))*(*(x + i));//7次方
		d[0] = d[0] + *(y + i);
		d[1] = d[1] + (*(x + i))*(*(y + i));
		d[2] = d[2] + (*(x + i))*(*(x + i))*(*(y + i));
		d[3] = d[3] + (*(x + i))*(*(x + i))*(*(x + i))*(*(y + i));
		d[4] = d[4] + (*(x + i))*(*(x + i))*(*(x + i))*(*(x + i))*(*(y + i));
	}
	q[4][3] = q[3][4];
	q[4][2] = q[3][3];
	q[2][4] = q[3][3];
	q[1][4] = q[2][3];
	q[4][1] = q[2][3];
	q[0][4] = q[2][2];
	q[4][0] = q[2][2];
	q[2][0] = q[1][1];
	q[2][1] = q[1][2];
	q[1][0] = q[0][1];
	q[0][3] = q[1][2];
	q[1][3] = q[2][2];
	q[3][1] = q[2][2];
	q[2][3] = q[3][2];
	q[3][0] = q[1][2];
	q[0][2] = q[1][1];
	ALU(5);
	//cout << "y=" << a[0] << "+" << a[1] << "x+" << a[2] << "x^2+" << a[3] << "x^4" << endl;
	return a[0] + a[1] * X + a[2] * (X*X) + a[3] * (X*X*X)+a[4]*(X*X*X*X);
}



int main()
{
	int n = 0;
	//cout << "输入数据个数: ";  //本次实验一为五个值
	//cin >> n;
	//cout << endl;
	double x[N];
	double y[N];
	double *xx;
	xx = x;  //存入x，f(x)值的数组
	double *yy;
	yy = y;
	ifstream xfile, yfile;
	xfile.open("xdata.txt");
	yfile.open("ydata.txt");//打开文本文件
	if (!xfile) cout << "xfile error" << endl;
	if (!yfile) cout << "yfile error" << endl;
	int nnn = 0;
	for (int i = 0; i < N; i++)
	{
		xfile >> x[i];
		yfile >> y[i];
		if (x[i] < 0)
		{
			n = i;     //xdata的数量
		}
		if (y[i] < 0)
		{
			nnn = i;
		}
		if (n != 0 && nnn != 0)
		{
			if (n == nnn)
				break;
			else
			{
				cout << "error" << endl;   //若xdata和ydata文件中读入的数据不一样多，则报错
				cin.get();
			}
		}
	}
	//cout << "n = " << n << endl;
	/*cout << "test x1=" << *(xx+1) << "   y1=" << *(yy+1) << endl;
	cout << "x0=" << *(xx) << endl;*/
	double aa, bb;
	/*cout << "请输入拉格朗日插值的x的值：";
		cin >> aa;
		cout << endl;
	cout<<"拉格朗日插值法得出的答案为："<<LagrangeModule(xx, yy,n,aa)<<endl; */
	/*cout << "请输入牛顿插值的x的值：";
	cin >> bb;
	cout << endl;
	cout << "牛顿插值法得出的答案为：" << NewtonModule(xx, yy, n, bb) << endl;*/
	ifstream setxfile;
	int nn;
	setxfile.open("setxdata.txt");
	double xdata[N];
	for (int i = 0; i < N; i++)
	{
		setxfile >> xdata[i];
		if (xdata[i] < 0)
		{
			nn = i;
			break;
		}
	}
	//cout << "nn = " << nn << endl;
	//cin.get();
	double ydata_2[N];
	double ydata_3[N];
	double ydata_4[N];
	cout <<"二阶多项式的结果为："<< erjie(xx, yy, n, 0.1) << " " << erjie(xx, yy, n, 0.3) << " " << erjie(xx, yy, n, 0.5) << " " << erjie(xx, yy, n, 0.7) << " " << erjie(xx, yy, n, 0.9) << endl;;
	cout << "三阶多项式的结果为：" << sanjie(xx, yy, n, 0.1) << " " << sanjie(xx, yy, n, 0.3) << " " << sanjie(xx, yy, n, 0.5) << " " << sanjie(xx, yy, n, 0.7)<<" "<< erjie(xx, yy, n, 0.9) <<endl;// << " " <<
	cout << "四阶多项式的结果为：" << sijie(xx, yy, n, 0.1) << " " << sijie(xx, yy, n, 0.3) << " " << sijie(xx, yy, n, 0.5) << " " << sijie(xx, yy, n, 0.7) << " " << sijie(xx, yy, n, 0.9) << endl;
	
	for (int i = 0; i < nn; i++)
	{
		ydata_2[i] = erjie(xx, yy, n, xdata[i]);
		ydata_3[i] = sanjie(xx, yy, n, xdata[i]);
		ydata_4[i] = sijie(xx, yy, n, xdata[i]);
	}
	
	ofstream erjied,sanjied,sijied;
	erjied.open("erjie.txt");
	if (!erjied) return 0;
	for (int i = 0; i < nn; i++)
		erjied << ydata_2[i] << " ";
	erjied.close();
	sanjied.open("sanjie.txt");
	if (!sanjied) return 0;
	for (int i = 0; i < nn; i++)
		sanjied << ydata_3[i] << " ";
	sanjied.close();
	sijied.open("sijie.txt");
	if (!sijied) return 0;
	for (int i = 0; i < nn; i++)
		sijied << ydata_4[i] << " ";
	sijied.close();
	cout << "计算完成，输出文件" << endl;
	cin.get(); cin.get();
	
}