#include<iostream>
#include<vector>
#include<cstring>
#include<fstream>
#include<string>
#include<omp.h>
#include <emmintrin.h>
#include <immintrin.h>
#include <windows.h>
#include<time.h>
using namespace std;
const int NUM_THREADS = 4;
//x��Ҫ�������ֵ
void bitmap(unsigned int* a,int x)
{
	size_t index = x >> 5;
	size_t num = x % 32;
	a[index] |= (1 << num);
}
void printb(unsigned int* a, int n)
{
	for (int i = n - 1; i >= 0; i--)
	{
		for (int j = 31; j >= 0; j--)
		{
			if ((a[i] >> j) & 1)
			{
				cout << 32 * i + j << "  ";
			}
		}
	}
}
/*
**����txt�ļ�
*/
void copyTxt(string srcFilename, string dstFilename)
{
	ifstream infile;
	ofstream outfile;
	string temp;
	infile.open(srcFilename, ios::in);
	outfile.open(dstFilename, ios::trunc | ios::out);
	if (infile.good())
	{
		while (!infile.eof())
		{
			getline(infile, temp, '\n');
			outfile << temp << '\n';
		}
	}
	infile.close();
	outfile.close();

}
/*
**�޸�ָ��������
*/
void ResetLine(string file, int beginline,int endline,unsigned int **a)
{
	string bup = ".\\tmp.txt";//�����ļ�
	copyTxt(file, bup);
	ifstream rfile;
	ofstream wfile;
	rfile.open(bup, ios::in);
	wfile.open(file, ios::out | ios::trunc);

	string str;
	int i = 1;
	while (!rfile.eof())
	{
		if (i == beginline)
		{
			
		}
		else
		{
			//rfile.getline()
			getline(rfile, str, '\n');
			wfile << str << '\n';
		}
		i++;
	}
	rfile.close();
	wfile.close();
}
////n�У�m��
//int** arrtobin(int** a, int n,int m)
//{
//
//}
int firstnum(unsigned int* a, int n)
{
	for (int i = n - 1; i >= 0; i--)
	{
		if (a[i] != 0)
		{
			unsigned int result = a[i];
			int num = 0;
			while (result != 0)
			{
				result >>= 1;
				num++;
			}
			int rs = 32 * i + num - 1;
			return rs;
		}
	}
}
bool isallzero(unsigned int* a, int n)
{
	int num = 0;
	for (int i = 0; i < n; i++)
	{
		if (a[i] == 0)
			continue;
		else
		{
			num = 1;
			break;
		}
	}
	return num;
}
void chuanxing(ifstream& xiaoyuanzi, string xiaoyuanzifile,ifstream& beixiaoyuanhang, int n, int m, int q)
{
	int block = 500;
	int size = (n >> 5) + 1;
	long long head, tail, freq; // timers
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	int Rsize = (m + q > block) ? (2*block) : (m + q);
	int Esize = (q > block) ? block : q;
	unsigned int last;
	beixiaoyuanhang >> last;
	while (!beixiaoyuanhang.eof())
	{
		unsigned int** E=new unsigned int* [Esize];
		for (int i = 0; i < Esize; i++)
		{
			E[i] = new unsigned[size];
			for (int j = 0; j < size; j++)
			{
				E[i][j] = 0;
			}
		}

		int row = 0;
		bitmap(E[row], last);
		while (row < block && !beixiaoyuanhang.eof())
		{
			unsigned int now;
			beixiaoyuanhang >> now;
			if (now > last)
			{
				if (row == block-1)
				{
					last = now;
					break;
				}
				row++;
			}
			last = now;
			bitmap(E[row], last);
		}
		row += 1;
		unsigned int last1;
		xiaoyuanzi >> last1;
		while (!xiaoyuanzi.eof())
		{
			unsigned int** R = new unsigned int* [Rsize];
			for (int i = 0; i < Rsize; i++)
			{
				R[i] = new unsigned[size];
				for (int j = 0; j < size; j++)
				{
					R[i][j] = 0;
				}
			}
			int row1 = 0;
			bitmap(R[row1], last1);
			//�˴ζ�ȡrow1����Ԫ��
			while (row1 < block && !xiaoyuanzi.eof())
			{
				unsigned int now;
				xiaoyuanzi >> now;
				if (now > last1)
				{
					if (row1 == block-1)
					{
						last1 = now;
						break;

					}
					row1++;
				}
				last1 = now;
				bitmap(R[row1], last1);
			}
			row1 += 1;
			for (int i = 0; i < row; i++)
			{
				if (isallzero(E[i],size))
				{
					int flag = 0;
					//row1����Ԫ�ӣ�֮��ǵü�
					for (int j = 0; j < row1; j++)
					{
						int first = firstnum(E[i],size);
						//д��
						if (firstnum(R[j], size) < first)
						{
							for (int a = row1; a > j; a--)
							{
								for (int b = 0; b <size; b++)
									R[a][b] = R[a - 1][b];
							}
							for (int b = 0; b < size; b++)
								R[j][b] = E[i][b];
							row1++;

							break;
						}
						if (firstnum(R[j],size) == first)
						{
							int minsize = size;
							for (int a = 0; a < minsize; a++)
							{
								int num = E[i][a] ^ R[j][a];
								E[i][a] = num;
							}

							flag = 1;
							//break;
						}
					}

				}
			}
			delete[]R;

		}
		xiaoyuanzi.close();
		//��Ԫ������open
		xiaoyuanzi.open(xiaoyuanzifile);
		delete[]E;
	}
	
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << (tail - head) * 1000.0 / freq << "          ";
}
void simd(ifstream& xiaoyuanzi, string xiaoyuanzifile, ifstream& beixiaoyuanhang, int n, int m, int q)
{
	int size = (n >> 5) + 1;
	long long head, tail, freq; // timers
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	int Rsize = (m + q > 500) ? 1000 : (m + q);
	int Esize = (q > 500) ? 500 : q;
	unsigned int last;
	beixiaoyuanhang >> last;
	while (!beixiaoyuanhang.eof())
	{
		unsigned int** E = new unsigned int* [Esize];
		for (int i = 0; i < Esize; i++)
		{
			E[i] = new unsigned[size];
			for (int j = 0; j < size; j++)
			{
				E[i][j] = 0;
			}
		}

		int row = 0;
		bitmap(E[row], last);
		while (row < 500 && !beixiaoyuanhang.eof())
		{
			unsigned int now;
			beixiaoyuanhang >> now;
			if (now > last)
			{
				if (row == 499)
				{
					last = now;
					break;
				}
				row++;
			}
			last = now;
			bitmap(E[row], last);
		}
		row += 1;
		unsigned int last1;
		xiaoyuanzi >> last1;
		while (!xiaoyuanzi.eof())
		{
			unsigned int** R = new unsigned int* [Rsize];
			for (int i = 0; i < Rsize; i++)
			{
				R[i] = new unsigned[size];
				for (int j = 0; j < size; j++)
				{
					R[i][j] = 0;
				}
			}
			int row1 = 0;
			bitmap(R[row1], last1);
			//�˴ζ�ȡrow1����Ԫ��
			while (row1 < 500 && !xiaoyuanzi.eof())
			{
				unsigned int now;
				xiaoyuanzi >> now;
				if (now > last1)
				{
					if (row1 == 499)
					{
						last1 = now;
						break;

					}
					row1++;
				}
				last1 = now;
				bitmap(R[row1], last1);
			}
			row1 += 1;
			for (int i = 0; i < row; i++)
			{
				if (isallzero(E[i], size))
				{
					//int flag = 0;
					//row1����Ԫ�ӣ�֮��ǵü�
					for (int j = 0; j < row1; j++)
					{
						int first = firstnum(E[i], size);
						//д��
						if (firstnum(R[j], size) < first)
						{
							for (int a = row1; a > j; a--)
							{
								for (int b = 0; b < size; b++)
									R[a][b] = R[a - 1][b];
							}
							for (int b = 0; b < size; b++)
								R[j][b] = E[i][b];
							row1++;

							break;
						}
						if (firstnum(R[j], size) == first)
						{
							int minsize = size;
							int a;
							for (a = 0; a + 4 <= minsize; a += 4)
							{
								__m128i Ri = _mm_loadu_si128((__m128i*)(&R[j][a]));
								__m128i Ei = _mm_loadu_si128((__m128i*)(&E[i][a]));
								__m128i Result = _mm_xor_si128(Ri, Ei);
								/*int num = bxyh[i].indexnum(a) ^ xyz[j].indexnum(a);
								bxyh[i].setindex(a, num);*/
								_mm_storeu_si128((__m128i*)(&E[i][a]), Result);
							}
							while (a < minsize)
							{
								int num = E[i][a] ^ R[j][a];
								E[i][a] = num;
								a++;
							}
						}
					}

				}
			}
			delete[]R;

		}
		xiaoyuanzi.close();
		//��Ԫ������open
		xiaoyuanzi.open(xiaoyuanzifile);
		delete[]E;
	}

	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << (tail - head) * 1000.0 / freq << "          ";
}
void omp(ifstream& xiaoyuanzi, string xiaoyuanzifile, ifstream& beixiaoyuanhang, int n, int m, int q)
{
	omp_set_num_threads(NUM_THREADS);
	int size = (n >> 5) + 1;
	long long head, tail, freq; // timers
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	int Rsize = (m + q > 500) ? 1000 : (m + q);
	int Esize = (q > 500) ? 500 : q;
	unsigned int last;
	beixiaoyuanhang >> last;
	while (!beixiaoyuanhang.eof())
	{
		unsigned int** E = new unsigned int* [Esize];
		for (int i = 0; i < Esize; i++)
		{
			E[i] = new unsigned[size];
			for (int j = 0; j < size; j++)
			{
				E[i][j] = 0;
			}
		}

		int row = 0;
		bitmap(E[row], last);
		while (row < 500 && !beixiaoyuanhang.eof())
		{
			unsigned int now;
			beixiaoyuanhang >> now;
			if (now > last)
			{
				if (row == 499)
				{
					last = now;
					break;
				}
				row++;
			}
			last = now;
			bitmap(E[row], last);
		}
		row += 1;
		unsigned int last1;
		xiaoyuanzi >> last1;
		while (!xiaoyuanzi.eof())
		{
			unsigned int** R = new unsigned int* [Rsize];
			for (int i = 0; i < Rsize; i++)
			{
				R[i] = new unsigned[size];
				for (int j = 0; j < size; j++)
				{
					R[i][j] = 0;
				}
			}
			int row1 = 0;
			bitmap(R[row1], last1);
			//�˴ζ�ȡrow1����Ԫ��
			while (row1 < 500 && !xiaoyuanzi.eof())
			{
				unsigned int now;
				xiaoyuanzi >> now;
				if (now > last1)
				{
					if (row1 == 499)
					{
						last1 = now;
						break;

					}
					row1++;
				}
				last1 = now;
				bitmap(R[row1], last1);
			}
			row1 += 1;
			for (int i = 0; i < row; i++)
			{
				if (isallzero(E[i], size))
				{
					int flag = 0;
					//row1����Ԫ�ӣ�֮��ǵü�
					for (int j = 0; j < row1; j++)
					{
						int first = firstnum(E[i], size);
						//д��
						if (firstnum(R[j], size) < first)
						{
#pragma omp parallel for
							for (int a = row1; a > j; a--)
							{
								for (int b = 0; b < size; b++)
									R[a][b] = R[a - 1][b];
							}
							for (int b = 0; b < size; b++)
								R[j][b] = E[i][b];
							row1++;

							break;
						}
						if (firstnum(R[j], size) == first)
						{
							int minsize = size;
							#pragma omp parallel for schedule(static,NUM_THREADS)
							for (int a = 0; a < minsize; a++)
							{
								int num = E[i][a] ^ R[j][a];
								E[i][a] = num;
							}
							flag = 1;
							//break;
						}
					}

				}
			}
			delete[]R;

		}
		xiaoyuanzi.close();
		//��Ԫ������open
		xiaoyuanzi.open(xiaoyuanzifile);
		delete[]E;
	}

	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << (tail - head) * 1000.0 / freq << "          ";
}
int main()
{
	//ifstream infile1("D:\\����\\Groebner\\��������1 ��������130��������Ԫ��22������Ԫ��8\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������1 ��������130��������Ԫ��22������Ԫ��8\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������1 ��������130��������Ԫ��22������Ԫ��8\\��Ԫ��.txt";
	//chuanxing(infile2, str,infile1, 130, 22, 8);
	//infile1.close();
	//infile2.close();
	//ifstream infile21("D:\\����\\Groebner\\��������1 ��������130��������Ԫ��22������Ԫ��8\\����Ԫ��.txt");
	//ifstream infile22("D:\\����\\Groebner\\��������1 ��������130��������Ԫ��22������Ԫ��8\\��Ԫ��.txt");
	//simd(infile22, str, infile21, 130, 22, 8);
	//infile21.close();
	//infile22.close();
	//ifstream infile31("D:\\����\\Groebner\\��������1 ��������130��������Ԫ��22������Ԫ��8\\����Ԫ��.txt");
	//ifstream infile32("D:\\����\\Groebner\\��������1 ��������130��������Ԫ��22������Ԫ��8\\��Ԫ��.txt");
	//omp(infile32, str, infile31, 130, 22, 8);


	//ifstream infile1("D:\\����\\Groebner\\��������2 ��������254��������Ԫ��106������Ԫ��53\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������2 ��������254��������Ԫ��106������Ԫ��53\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������2 ��������254��������Ԫ��106������Ԫ��53\\��Ԫ��.txt";
	//chuanxing(infile2, str,infile1, 254, 106, 53);
	//infile1.close();
	//infile2.close();
	//ifstream infile21("D:\\����\\Groebner\\��������2 ��������254��������Ԫ��106������Ԫ��53\\����Ԫ��.txt");
	//ifstream infile22("D:\\����\\Groebner\\��������2 ��������254��������Ԫ��106������Ԫ��53\\��Ԫ��.txt");
	//simd(infile22, str, infile21, 254, 106, 53);
	//infile21.close();
	//infile22.close();
	//ifstream infile31("D:\\����\\Groebner\\��������2 ��������254��������Ԫ��106������Ԫ��53\\����Ԫ��.txt");
	//ifstream infile32("D:\\����\\Groebner\\��������2 ��������254��������Ԫ��106������Ԫ��53\\��Ԫ��.txt");
	//omp(infile32, str, infile31, 254, 106, 53);


	//ifstream infile1("D:\\����\\Groebner\\��������3 ��������562��������Ԫ��170������Ԫ��53\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������3 ��������562��������Ԫ��170������Ԫ��53\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������3 ��������562��������Ԫ��170������Ԫ��53\\��Ԫ��.txt";
	//chuanxing(infile2,str,infile1, 562, 170, 53);
	//infile1.close();
	//infile2.close();
	//ifstream infile21("D:\\����\\Groebner\\��������3 ��������562��������Ԫ��170������Ԫ��53\\����Ԫ��.txt");
	//ifstream infile22("D:\\����\\Groebner\\��������3 ��������562��������Ԫ��170������Ԫ��53\\��Ԫ��.txt");
	//simd(infile22, str, infile21, 562, 170, 53);
	//infile21.close();
	//infile22.close();
	//ifstream infile31("D:\\����\\Groebner\\��������3 ��������562��������Ԫ��170������Ԫ��53\\����Ԫ��.txt");
	//ifstream infile32("D:\\����\\Groebner\\��������3 ��������562��������Ԫ��170������Ԫ��53\\��Ԫ��.txt");
	//omp(infile32, str, infile31, 562, 170, 53);



	//ifstream infile1("D:\\����\\Groebner\\��������4 ��������1011��������Ԫ��539������Ԫ��263\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������4 ��������1011��������Ԫ��539������Ԫ��263\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������4 ��������1011��������Ԫ��539������Ԫ��263\\��Ԫ��.txt";
	//chuanxing(infile2,str, infile1, 1011, 539, 263);
	//infile1.close();
	//infile2.close();
	//ifstream infile21("D:\\����\\Groebner\\��������4 ��������1011��������Ԫ��539������Ԫ��263\\����Ԫ��.txt");
	//ifstream infile22("D:\\����\\Groebner\\��������4 ��������1011��������Ԫ��539������Ԫ��263\\��Ԫ��.txt");
	//simd(infile22, str,infile21, 1011, 539, 263);
	//infile21.close();
	//infile22.close();
	//ifstream infile31("D:\\����\\Groebner\\��������4 ��������1011��������Ԫ��539������Ԫ��263\\����Ԫ��.txt");
	//ifstream infile32("D:\\����\\Groebner\\��������4 ��������1011��������Ԫ��539������Ԫ��263\\��Ԫ��.txt");
	//omp(infile32, str, infile31, 1011, 539, 263);

	//ifstream infile1("D:\\����\\Groebner\\��������5 ��������2362��������Ԫ��1226������Ԫ��453\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������5 ��������2362��������Ԫ��1226������Ԫ��453\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������5 ��������2362��������Ԫ��1226������Ԫ��453\\��Ԫ��.txt";
	//chuanxing(infile2,str, infile1,2362,1226,453);
	//infile1.close();
	//infile2.close();
	//ifstream infile21("D:\\����\\Groebner\\��������5 ��������2362��������Ԫ��1226������Ԫ��453\\����Ԫ��.txt");
	//ifstream infile22("D:\\����\\Groebner\\��������5 ��������2362��������Ԫ��1226������Ԫ��453\\��Ԫ��.txt");
	//simd(infile22, str,infile21, 2362, 1226, 453);
	//infile21.close();
	//infile22.close();
	//ifstream infile31("D:\\����\\Groebner\\��������5 ��������2362��������Ԫ��1226������Ԫ��453\\����Ԫ��.txt");
	//ifstream infile32("D:\\����\\Groebner\\��������5 ��������2362��������Ԫ��1226������Ԫ��453\\��Ԫ��.txt");
	//omp(infile32, str, infile31, 2362, 1226, 453);

	
	//ifstream infile1("D:\\����\\Groebner\\��������6 ��������3799��������Ԫ��2759������Ԫ��1953\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������6 ��������3799��������Ԫ��2759������Ԫ��1953\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������6 ��������3799��������Ԫ��2759������Ԫ��1953\\��Ԫ��.txt";
	//chuanxing(infile2, str,infile1, 3799, 2759, 1953);
	//infile1.close();
	//infile2.close();
	//ifstream infile21("D:\\����\\Groebner\\��������6 ��������3799��������Ԫ��2759������Ԫ��1953\\����Ԫ��.txt");
	//ifstream infile22("D:\\����\\Groebner\\��������6 ��������3799��������Ԫ��2759������Ԫ��1953\\��Ԫ��.txt");
	//simd(infile22, str,infile21, 3799, 2759, 1953);
	//infile21.close();
	//infile22.close();
	//ifstream infile31("D:\\����\\Groebner\\��������6 ��������3799��������Ԫ��2759������Ԫ��1953\\����Ԫ��.txt");
	//ifstream infile32("D:\\����\\Groebner\\��������6 ��������3799��������Ԫ��2759������Ԫ��1953\\��Ԫ��.txt");
	//omp(infile32, str, infile31, 3799, 2759, 1953);



	//ifstream infile1("D:\\����\\Groebner\\��������7 ��������8399��������Ԫ��6375������Ԫ��4535\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������7 ��������8399��������Ԫ��6375������Ԫ��4535\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������7 ��������8399��������Ԫ��6375������Ԫ��4535\\��Ԫ��.txt";
	//chuanxing(infile2, str,infile1, 8399, 6375, 4535);
	//infile1.close();
	//infile2.close();
	//ifstream infile21("D:\\����\\Groebner\\��������7 ��������8399��������Ԫ��6375������Ԫ��4535\\����Ԫ��.txt");
	//ifstream infile22("D:\\����\\Groebner\\��������7 ��������8399��������Ԫ��6375������Ԫ��4535\\��Ԫ��.txt");
	//simd(infile22, str,infile21, 8399, 6375, 4535);
	//infile21.close();
	//infile22.close();
	//ifstream infile31("D:\\����\\Groebner\\��������7 ��������8399��������Ԫ��6375������Ԫ��4535\\����Ԫ��.txt");
	//ifstream infile32("D:\\����\\Groebner\\��������7 ��������8399��������Ԫ��6375������Ԫ��4535\\��Ԫ��.txt");
	//omp(infile32, str, infile31, 8399, 6375, 4535);
///û�ܳ���
	ifstream infile1("D:\\����\\Groebner\\��������8 ��������23045��������Ԫ��18748������Ԫ��14325\\����Ԫ��.txt");
	ifstream infile2("D:\\����\\Groebner\\��������8 ��������23045��������Ԫ��18748������Ԫ��14325\\��Ԫ��.txt");
	string str = "D:\\����\\Groebner\\��������8 ��������23045��������Ԫ��18748������Ԫ��14325\\��Ԫ��.txt";
	chuanxing(infile2, str, infile1, 23045, 18748, 14325);




	//ifstream infile1("D:\\����\\Groebner\\��������11 ��������85401��������Ԫ��5724������Ԫ��756\\����Ԫ��.txt");
	//ifstream infile2("D:\\����\\Groebner\\��������11 ��������85401��������Ԫ��5724������Ԫ��756\\��Ԫ��.txt");
	//string str = "D:\\����\\Groebner\\��������11 ��������85401��������Ԫ��5724������Ԫ��756\\��Ԫ��.txt";
	//chuanxing(infile2, str, infile1, 85401, 5724, 756);

}