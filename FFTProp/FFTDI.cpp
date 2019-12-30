#include "FFTDI.h"
#include <fstream>
#include "../common/common/include/GEMS_Memory.h"
#include "../common/common/include/constant_var.h"
#include "FFT.h"
#include <vector>
using namespace std;
using namespace Common;

FFTDI::FFTDI(double f, double xp, double yp, double zp, int N, int M)
	:f(f),xp(xp),yp(yp),zp(zp),N(N),M(M)
{
	lamda = C_Speed / f;
	k = 2 * Pi * f / C_Speed;
	ds = lamda / 3.5;
	N2 = 2 * N - 1;
	M2 = 2 * M - 1;

	PowerHExy = 0;
	PowerHEz_x0 = 0;
	PowerHEz_y0 = 0;
	PowerHHx_x0 = 0;
	PowerHHx_y0 = 0;
	PowerHHy_x0 = 0;
	PowerHHy_y0 = 0;
	PowerHHz_x0 = 0;
	PowerHHz_y0 = 0;

	PowerFFTHExy = 0;
	PowerFFTHEz_x0 = 0;
	PowerFFTHEz_y0 = 0;
	PowerFFTHHx_x0 = 0;
	PowerFFTHHx_y0 = 0;
	PowerFFTHHy_x0 = 0;
	PowerFFTHHy_y0 = 0;
	PowerFFTHHz_x0 = 0;
	PowerFFTHHz_y0 = 0;

	threadNum = 4;

	Initialization();
}

FFTDI::~FFTDI()
{
	FreeCal();
}

void FFTDI::FreeCal()
{

	Free_2D(Ex0);
	Free_2D(Ey0);

	Free_2D(Ex1);
	Free_2D(Ey1);
	Free_2D(Ez1);
	Free_2D(Hx1);
	Free_2D(Hy1);
	Free_2D(Hz1);
}

void FFTDI::Initialization()
{
	//��0 ������-���볡�ֲ�
	Ex0 = Allocate_2D(Ex0, N2, M2);
	Ey0 = Allocate_2D(Ey0, N2, M2);
	for (int i = 0; i < N2; i++)
	{
		for (int j = 0; j < M2; j++)
		{
			Ex0[i][j] = (0, 0);			Ey0[i][j] = (0, 0);
		}
	}

	//Ex1 Ey1 Ez1 Hx1 Hy1 Hz1
	Ex1 = Allocate_2D(Ex1, N, M);
	Ey1 = Allocate_2D(Ey1, N, M);
	Ez1 = Allocate_2D(Ez1, N, M);
	Hx1 = Allocate_2D(Hx1, N, M);
	Hy1 = Allocate_2D(Hy1, N, M);
	Hz1 = Allocate_2D(Hz1, N, M);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Ex1[i][j] = (0, 0);			Ey1[i][j] = (0, 0);			Ez1[i][j] = (0, 0);
			Hx1[i][j] = (0, 0);			Hy1[i][j] = (0, 0);			Hz1[i][j] = (0, 0);
		}
	}


	/*
	//���ݺ���
	HExy = Allocate_2D(HExy, N2, M2);
	HEz_x0 = Allocate_2D(HEz_x0, N2, M2);
	HEz_y0 = Allocate_2D(HEz_y0, N2, M2);
	HHx_x0 = Allocate_2D(HHx_x0, N2, M2);
	HHx_y0 = Allocate_2D(HHx_y0, N2, M2);
	HHy_x0 = Allocate_2D(HHy_x0, N2, M2);
	HHy_y0 = Allocate_2D(HHy_y0, N2, M2);
	HHz_x0 = Allocate_2D(HHz_x0, N2, M2);
	HHz_y0 = Allocate_2D(HHz_y0, N2, M2);

	


	for (int i = 0; i < N2; i++)
	{
		for (int j = 0; j < M2; j++)
		{
			HExy[i][j] = (0, 0);
			HEz_x0[i][j] = (0, 0);
			HEz_y0[i][j] = (0, 0);
			HHx_x0[i][j] = (0, 0);
			HHx_y0[i][j] = (0, 0);
			HHy_x0[i][j] = (0, 0);
			HHy_y0[i][j] = (0, 0);
			HHz_x0[i][j] = (0, 0);
			HHz_y0[i][j] = (0, 0);
		}
	}
	*/
}

void FFTDI::Setds(double ds1)
{
	ds = ds1;
}

void FFTDI::SetInput(std::complex<double> ** ExIn, std::complex<double> ** EyIn)
{
	for (int i = 0; i < N2; i++)
		for (int j = 0; j < M2; j++)
		{
			//int tempi, tempj;
			//tempi = i + (N - 1) / 2;
			//tempj = j + (N - 1) / 2;
			if ((i < N) && (j < M))
			{
				Ex0[i][j] = ExIn[i][j];
				Ey0[i][j] = EyIn[i][j];
			}
			else
			{
				Ex0[i][j] = 0;
				Ey0[i][j] = 0;
			}
		/*	double tempA, tempB, tempR, tempI;
			tempA = 2 * Pi / 2 * 2 * (N - 1) / (2 * (N - 1) - 1) * tempi;
			tempB = 2 * Pi / 2 * 2 * (N - 1) / (2 * (N - 1) - 1) * tempj;
			MulComp(cos(tempA), sin(tempA), cos(tempB), sin(tempB), tempR, tempI);
			Ex0[tempi][tempj] = MulComp(Ex0[tempi][tempj], tempR, tempI);
			Ey0[tempi][tempj] = MulComp(Ey0[tempi][tempj], tempR, tempI);
*/
	}

}


void FFTDI::SourceFFT() {
	FFT fft;
	fft.FFT_2D(Ex0, Ex0, N2, M2);
	fft.FFT_2D(Ey0, Ey0, N2, M2);
}

void FFTDI::FFTToExy1() {
//׼����ʱ�����ʹ��ݺ�������
	tempEx1 = Allocate_2D(tempEx1, N2, M2);//���������
	tempEy1 = Allocate_2D(tempEy1, N2, M2);
	HExy = Allocate_2D(HExy, N2, M2);	   //���ݺ���
	for (int i = 0; i < N2; i++)
		for (int j = 0; j < M2; j++) {
			tempEx1[i][j] = (0, 0);	tempEy1[i][j] = (0, 0);
			HExy[i][j] = (0, 0);
		}

	PowerHExy = 0;
//���㴫�ݺ���
	vector<int> starti;
	vector<int> numi;
	vector<int> stopi;
	starti.resize(threadNum);	numi.resize(threadNum);	stopi.resize(threadNum);
	for (int i = 0; i < threadNum; i++) {
		int num = N2 / threadNum;
		int rr = (N2) % threadNum;
		if (i<rr && rr != 0) num += 1;
		numi[i] = num;
		//starti
		if (i == 0) starti[i] = 0;
		else starti[i] = starti[i - 1] + numi[i - 1];
		//stopi
		stopi[i] = starti[i] + numi[i];
	}
	vector<double>PowerHExyList; PowerHExyList.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHExyList[i] = 0; }
	vector<double>PowerFFTHExyList; PowerFFTHExyList.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHExyList[i] = 0; }

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				double x, y, z; //���޸�
				x = (i - N + 1) * ds + xp;
				//x = (i - N + 1) * ds;
				y = (j - N + 1) * ds + yp;
				z = zp;
				//z = z0 ;

				//�м����
				//double k1, k2; //ϵ��
				double Rejkr, Iejkr; //exp(-jkr)
				double tempR1, tempI1;
				double tempR2, tempI2;
				double r = pow((x * x + y * y + z * z), 0.5); //����
				double r5 = pow(r, 5);
				double r3 = pow(r, 3);

				Rejkr = cos(k * r);
				Iejkr = -sin(k * r);

				MulComp(Rejkr, Iejkr, 1, k*r, tempR1, tempI1);
				tempR1 = tempR1 / r3 * z / 2 / Pi;
				tempI1 = tempI1 / r3 * z / 2 / Pi;

				double temp = z / (2 * Pi * r * r * r);
				tempR1 = temp * (cos(k * r) + k * r * sin(k * r));
				tempI1 = temp * (k * r * cos(k * r) - sin(k * r));

				HExy[i][j] = std::complex<double>(tempR1, tempI1);

				PowerHExyList[id] = PowerHExyList[id] + ds * ds * (HExy[i][j].real() * HExy[i][j].real() + HExy[i][j].imag() * HExy[i][j].imag());
			}
		}//endloop

	}//openmp
	for (int i = 0; i < threadNum; i++) {
		PowerHExy = PowerHExyList[i] + PowerHExy;
	}
	//ת��������
	FFT fft;
	fft.FFT_2D(HExy, HExy, N2, M2);

	//����������
	PowerFFTHExy = 0;

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{	
		double temppowH0 = pow((1.0 / (2 * N - 1) / ds), 2);//���ڼ���power����ʱ����
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				//Ϊ����һ��ѭ�� �ڴ˼���fft���Power
				PowerFFTHExy = PowerFFTHExy + temppowH0 * (HExy[i][j].real() * HExy[i][j].real() + HExy[i][j].imag() * HExy[i][j].imag());
			}
		}
	}
	for (int i = 0; i < threadNum; i++) {
		PowerFFTHExy = PowerFFTHExyList[i] + PowerFFTHExy;
	}

	//��������
	double temppowHExy = pow((PowerHExy / PowerFFTHExy), 0.5);

	//�������-����
	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				tempEx1[i][j] = std::complex<double>((Ex0[i][j].real() * HExy[i][j].real() - Ex0[i][j].imag() * HExy[i][j].imag()) * temppowHExy,
					(Ex0[i][j].real() * HExy[i][j].imag() + Ex0[i][j].imag() * HExy[i][j].real()) * temppowHExy);
				tempEy1[i][j] = std::complex<double>((Ey0[i][j].real() * HExy[i][j].real() - Ey0[i][j].imag() * HExy[i][j].imag()) * temppowHExy,
					(Ey0[i][j].real() * HExy[i][j].imag() + Ey0[i][j].imag() * HExy[i][j].real()) * temppowHExy);
			}
		}
	}
	
	//��任-��
	fft.IFFT_2D(tempEx1, tempEx1, N2, M2);
	fft.IFFT_2D(tempEy1, tempEy1, N2, M2);

	//���
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Ex1[i][j] = tempEx1[i + N - 1][j + M - 1];
			Ey1[i][j] = tempEy1[i + N - 1][j + M - 1];
		}
	}

	Free_2D(tempEx1);		Free_2D(tempEy1);	Free_2D(HExy);
}

void FFTDI::FFTToEz1() {
	//׼����ʱ�����ʹ��ݺ�������
	tempEz1 = Allocate_2D(tempEz1, N2, M2);//���������
	HEz_x0 = Allocate_2D(HEz_x0, N2, M2);	   //���ݺ���
	HEz_y0 = Allocate_2D(HEz_y0, N2, M2);
	for (int i = 0; i < N2; i++)
		for (int j = 0; j < M2; j++) {
			tempEz1[i][j] = (0, 0);	
			HEz_x0[i][j] = (0, 0);
			HEz_y0[i][j] = (0, 0);
		}

	PowerHEz_x0 = 0;	PowerHEz_y0 = 0;
	//���㴫�ݺ���
	vector<int> starti;
	vector<int> numi;
	vector<int> stopi;
	starti.resize(threadNum);	numi.resize(threadNum);	stopi.resize(threadNum);
	for (int i = 0; i < threadNum; i++) {
		int num = N2 / threadNum;
		int rr = (N2) % threadNum;
		if (i<rr && rr != 0) num += 1;
		numi[i] = num;
		//starti
		if (i == 0) starti[i] = 0;
		else starti[i] = starti[i - 1] + numi[i - 1];
		//stopi
		stopi[i] = starti[i] + numi[i];
	}
	vector<double>PowerHEz_x0List; PowerHEz_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHEz_x0List[i] = 0; }
	vector<double>PowerHEz_y0List; PowerHEz_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHEz_y0List[i] = 0; }
	vector<double>PowerFFTHEz_x0List; PowerFFTHEz_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHEz_x0List[i] = 0; }
	vector<double>PowerFFTHEz_y0List; PowerFFTHEz_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHEz_y0List[i] = 0; }

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				double x, y, z; //���޸�
				x = (i - N + 1) * ds + xp;
				//x = (i - N + 1) * ds;
				y = (j - N + 1) * ds + yp;
				z = zp;
				//z = z0 ;

				//�м����
				//double k1, k2; //ϵ��
				double Rejkr, Iejkr; //exp(-jkr)
				double tempR1, tempI1;
				double tempR2, tempI2;
				double r = pow((x * x + y * y + z * z), 0.5); //����
				double r5 = pow(r, 5);
				double r3 = pow(r, 3);

				Rejkr = cos(k * r);
				Iejkr = -sin(k * r);

				MulComp(Rejkr, Iejkr, 1, k*r, tempR1, tempI1);
				tempR1 = tempR1 / r3 * z / 2 / Pi;
				tempI1 = tempI1 / r3 * z / 2 / Pi;

				double temp = z / (2 * Pi * r * r * r);
				tempR1 = temp * (cos(k * r) + k * r * sin(k * r));
				tempI1 = temp * (k * r * cos(k * r) - sin(k * r));

				MulComp(Rejkr, Iejkr, -1, -k*r, tempR1, tempI1);
				tempR2 = tempR1;
				tempI2 = tempI1;
				tempR1 = tempR1 / r3 * x / 2 / Pi;
				tempI1 = tempI1 / r3 * x / 2 / Pi;
				tempR2 = tempR2 / r3 * y / 2 / Pi;
				tempI2 = tempI2 / r3 * y / 2 / Pi;
				HEz_x0[i][j] = std::complex<double>(tempR1, tempI1);
				HEz_y0[i][j] = std::complex<double>(tempR2, tempI2);

				PowerHEz_x0List[id] = PowerHEz_x0List[id] + ds * ds * (HEz_x0[i][j].real() * HEz_x0[i][j].real() + HEz_x0[i][j].imag() * HEz_x0[i][j].imag());
				PowerHEz_y0List[id] = PowerHEz_y0List[id] + ds * ds * (HEz_y0[i][j].real() * HEz_y0[i][j].real() + HEz_y0[i][j].imag() * HEz_y0[i][j].imag());

			}
		}//endloop

	}//openmp
	for (int i = 0; i < threadNum; i++) {
		PowerHEz_x0 = PowerHEz_x0List[i] + PowerHEz_x0;
		PowerHEz_y0 = PowerHEz_y0List[i] + PowerHEz_y0;
	}

	//ת��������
	FFT fft;
	fft.FFT_2D(HEz_x0, HEz_x0, N2, M2);
	fft.FFT_2D(HEz_y0, HEz_y0, N2, M2);

	//����������
	PowerFFTHEz_x0 = 0;
	PowerFFTHEz_y0 = 0;

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		double temppowH0 = pow((1.0 / (2 * N - 1) / ds), 2);//���ڼ���power����ʱ����
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				//Ϊ����һ��ѭ�� �ڴ˼���fft���Power
				PowerFFTHEz_x0 = PowerFFTHEz_x0 + temppowH0 * (HEz_x0[i][j].real() * HEz_x0[i][j].real() + HEz_x0[i][j].imag() * HEz_x0[i][j].imag());
				PowerFFTHEz_y0 = PowerFFTHEz_y0 + temppowH0 * (HEz_y0[i][j].real() * HEz_y0[i][j].real() + HEz_y0[i][j].imag() * HEz_y0[i][j].imag());
			}
		}
	}
	for (int i = 0; i < threadNum; i++) {
		PowerFFTHEz_x0 = PowerFFTHEz_x0List[i] + PowerFFTHEz_x0;
		PowerFFTHEz_y0 = PowerFFTHEz_y0List[i] + PowerFFTHEz_y0;
	}

	//��������
	double temppowHEz_x0 = pow((PowerHEz_x0 / PowerFFTHEz_x0), 0.5);
	double temppowHEz_y0 = pow((PowerHEz_y0 / PowerFFTHEz_y0), 0.5);

	//�������-����
	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{

				double tempR1, tempI1; //��ʽ��һ����
				double tempR2, tempI2; //��ʽ�ڶ�����
									   //Ez1
				MulComp(Ex0[i][j].real(), Ex0[i][j].imag(), HEz_x0[i][j].real(), HEz_x0[i][j].imag(), tempR1, tempI1);
				MulComp(Ey0[i][j].real(), Ey0[i][j].imag(), HEz_y0[i][j].real(), HEz_y0[i][j].imag(), tempR2, tempI2);
				tempR1 = tempR1 * temppowHEz_x0;
				tempI1 = tempI1 * temppowHEz_x0;
				tempR2 = tempR2 * temppowHEz_y0;
				tempI2 = tempI2 * temppowHEz_y0;
				tempEz1[i][j] = std::complex<double>(tempR1 + tempR2, tempI1 + tempI2);

			}
		}
	}

	//��任-��
	fft.IFFT_2D(tempEz1, tempEz1, N2, M2);

	//���
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Ez1[i][j] = tempEz1[i + N - 1][j + M - 1];
		}
	}

	//�����ڴ�
	Free_2D(tempEz1);	Free_2D(HEz_x0);	Free_2D(HEz_y0);
}

void FFTDI::FFTToHx1() {
	//׼����ʱ�����ʹ��ݺ�������
	tempHx1 = Allocate_2D(tempHx1, N2, M2);//���������
	HHx_x0 = Allocate_2D(HHx_x0, N2, M2);	   //���ݺ���
	HHx_y0 = Allocate_2D(HHx_y0, N2, M2);
	for (int i = 0; i < N2; i++)
		for (int j = 0; j < M2; j++) {
			tempHx1[i][j] = (0, 0);
			HHx_x0[i][j] = (0, 0);
			HHx_y0[i][j] = (0, 0);
		}

	PowerHHx_x0 = 0;	PowerHHx_y0 = 0;
	//���㴫�ݺ���
	vector<int> starti;
	vector<int> numi;
	vector<int> stopi;
	starti.resize(threadNum);	numi.resize(threadNum);	stopi.resize(threadNum);
	for (int i = 0; i < threadNum; i++) {
		int num = N2 / threadNum;
		int rr = (N2) % threadNum;
		if (i<rr && rr != 0) num += 1;
		numi[i] = num;
		//starti
		if (i == 0) starti[i] = 0;
		else starti[i] = starti[i - 1] + numi[i - 1];
		//stopi
		stopi[i] = starti[i] + numi[i];
	}
	vector<double>PowerHHx_x0List; PowerHHx_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHHx_x0List[i] = 0; }
	vector<double>PowerHHx_y0List; PowerHHx_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHHx_y0List[i] = 0; }
	vector<double>PowerFFTHHx_x0List; PowerFFTHHx_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHHx_x0List[i] = 0; }
	vector<double>PowerFFTHHx_y0List; PowerFFTHHx_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHHx_y0List[i] = 0; }


	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				double x, y, z; //���޸�
				x = (i - N + 1) * ds + xp;
				//x = (i - N + 1) * ds;
				y = (j - N + 1) * ds + yp;
				z = zp;
				//z = z0 ;

				//�м����
				//double k1, k2; //ϵ��
				double Rejkr, Iejkr; //exp(-jkr)
				double tempR1, tempI1;
				double tempR2, tempI2;
				double r = pow((x * x + y * y + z * z), 0.5); //����
				double r5 = pow(r, 5);
				double r3 = pow(r, 3);

				Rejkr = cos(k * r);
				Iejkr = -sin(k * r);

				MulComp(Rejkr, Iejkr, 1, k*r, tempR1, tempI1);
				tempR1 = tempR1 / r3 * z / 2 / Pi;
				tempI1 = tempI1 / r3 * z / 2 / Pi;

				double temp = z / (2 * Pi * r * r * r);
				tempR1 = temp * (cos(k * r) + k * r * sin(k * r));
				tempI1 = temp * (k * r * cos(k * r) - sin(k * r));

				MulComp(Rejkr, Iejkr, -1, -k*r, tempR1, tempI1);
				tempR2 = tempR1;
				tempI2 = tempI1;
				tempR1 = tempR1 / r3 * x / 2 / Pi;
				tempI1 = tempI1 / r3 * x / 2 / Pi;
				tempR2 = tempR2 / r3 * y / 2 / Pi;
				tempI2 = tempI2 / r3 * y / 2 / Pi;
				double pi2wu0 = 2 * Pi * 2 * Pi * Mu0 * f; // �м�ֵ = 2*pi*w*u0 w = 2*pi*f

				MulComp(Rejkr, Iejkr, -3 * k * r, 3 - k * k * r * r, tempR1, tempI1);
				MulComp(sin(k * r), cos(k * r), 1, k * r, tempR2, tempI2);

				HHx_x0[i][j] = std::complex<double>(tempR1 / r5 * x * y / pi2wu0,
					tempI1 / r5 * x * y / pi2wu0);

				HHx_y0[i][j] = std::complex<double>((tempR1 / r5 * (y * y + z * z) - tempR2 / r3 * 2) / pi2wu0,
					(tempI1 / r5 * (y * y + z * z) - tempI2 / r3 * 2) / pi2wu0);

				PowerHHx_x0List[id] = PowerHHx_x0List[id] + ds * ds * (HHx_x0[i][j].real() * HHx_x0[i][j].real() + HHx_x0[i][j].imag() * HHx_x0[i][j].imag());
				PowerHHx_y0List[id] = PowerHHx_y0List[id] + ds * ds * (HHx_y0[i][j].real() * HHx_y0[i][j].real() + HHx_y0[i][j].imag() * HHx_y0[i][j].imag());

			}
		}//endloop

	}//openmp
	for (int i = 0; i < threadNum; i++) {
		PowerHHx_x0 = PowerHHx_x0List[i] + PowerHHx_x0;
		PowerHHx_y0 = PowerHHx_y0List[i] + PowerHHx_y0;
	}

	//ת��������
	FFT fft;
	fft.FFT_2D(HHx_x0, HHx_x0, N2, M2);
	fft.FFT_2D(HHx_y0, HHx_y0, N2, M2);

	//����������
	PowerFFTHHx_x0 = 0;
	PowerFFTHHx_y0 = 0;

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		double temppowH0 = pow((1.0 / (2 * N - 1) / ds), 2);//���ڼ���power����ʱ����
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				//Ϊ����һ��ѭ�� �ڴ˼���fft���Power
				PowerFFTHHx_x0List[id] = PowerFFTHHx_x0List[id] + temppowH0 * (HHx_x0[i][j].real() * HHx_x0[i][j].real() + HHx_x0[i][j].imag() * HHx_x0[i][j].imag());
				PowerFFTHHx_y0List[id] = PowerFFTHHx_y0List[id] + temppowH0 * (HHx_y0[i][j].real() * HHx_y0[i][j].real() + HHx_y0[i][j].imag() * HHx_y0[i][j].imag());
			}
		}
	}
	for (int i = 0; i < threadNum; i++) {
		PowerFFTHHx_x0 = PowerFFTHHx_x0List[i] + PowerFFTHHx_x0;
		PowerFFTHHx_y0 = PowerFFTHHx_y0List[i] + PowerFFTHHx_y0;
	}

	//��������
	double temppowHHx_x0 = pow((PowerHHx_x0 / PowerFFTHHx_x0), 0.5);
	double temppowHHx_y0 = pow((PowerHHx_y0 / PowerFFTHHx_y0), 0.5);

	//�������-����
	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				double tempR1, tempI1; //��ʽ��һ����
				double tempR2, tempI2; //��ʽ�ڶ�����
				//Hx1
				MulComp(Ex0[i][j].real(), Ex0[i][j].imag(), HHx_x0[i][j].real(), HHx_x0[i][j].imag(), tempR1, tempI1);
				MulComp(Ey0[i][j].real(), Ey0[i][j].imag(), HHx_y0[i][j].real(), HHx_y0[i][j].imag(), tempR2, tempI2);
				tempR1 = tempR1 * temppowHHx_x0;
				tempI1 = tempI1 * temppowHHx_x0;
				tempR2 = tempR2 * temppowHHx_y0;
				tempI2 = tempI2 * temppowHHx_y0;
				tempHx1[i][j] = std::complex<double>(tempR1 + tempR2, tempI1 + tempI2);

			}
		}
	}

	//��任-��
	fft.IFFT_2D(tempHx1, tempHx1, N2, M2);

	//���
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Hx1[i][j] = tempHx1[i + N - 1][j + M - 1];
		}
	}

	//�����ڴ�
	Free_2D(tempHx1);	Free_2D(HHx_x0);	Free_2D(HHx_y0);
}

void FFTDI::FFTToHy1() {
	//׼����ʱ�����ʹ��ݺ�������
	tempHy1 = Allocate_2D(tempHy1, N2, M2);//���������
	HHy_x0 = Allocate_2D(HHy_x0, N2, M2);	   //���ݺ���
	HHy_y0 = Allocate_2D(HHy_y0, N2, M2);
	for (int i = 0; i < N2; i++)
		for (int j = 0; j < M2; j++) {
			tempHy1[i][j] = (0, 0);
			HHy_x0[i][j] = (0, 0);
			HHy_y0[i][j] = (0, 0);
		}

	PowerHHy_x0 = 0;	PowerHHy_y0 = 0;
	//���㴫�ݺ���
	vector<int> starti;
	vector<int> numi;
	vector<int> stopi;
	starti.resize(threadNum);	numi.resize(threadNum);	stopi.resize(threadNum);
	for (int i = 0; i < threadNum; i++) {
		int num = N2 / threadNum;
		int rr = (N2) % threadNum;
		if (i < rr && rr != 0) num += 1;
		numi[i] = num;
		//starti
		if (i == 0) starti[i] = 0;
		else starti[i] = starti[i - 1] + numi[i - 1];
		//stopi
		stopi[i] = starti[i] + numi[i];
	}
	vector<double>PowerHHy_x0List; PowerHHy_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHHy_x0List[i] = 0; }
	vector<double>PowerHHy_y0List; PowerHHy_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHHy_y0List[i] = 0; }
	vector<double>PowerFFTHHy_x0List; PowerFFTHHy_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHHy_x0List[i] = 0; }
	vector<double>PowerFFTHHy_y0List; PowerFFTHHy_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHHy_y0List[i] = 0; }

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				double x, y, z; //���޸�
				x = (i - N + 1) * ds + xp;
				//x = (i - N + 1) * ds;
				y = (j - N + 1) * ds + yp;
				z = zp;
				//z = z0 ;

				//�м����
				//double k1, k2; //ϵ��
				double Rejkr, Iejkr; //exp(-jkr)
				double tempR1, tempI1;
				double tempR2, tempI2;
				double r = pow((x * x + y * y + z * z), 0.5); //����
				double r5 = pow(r, 5);
				double r3 = pow(r, 3);

				Rejkr = cos(k * r);
				Iejkr = -sin(k * r);

				MulComp(Rejkr, Iejkr, 1, k*r, tempR1, tempI1);
				tempR1 = tempR1 / r3 * z / 2 / Pi;
				tempI1 = tempI1 / r3 * z / 2 / Pi;

				double temp = z / (2 * Pi * r * r * r);
				tempR1 = temp * (cos(k * r) + k * r * sin(k * r));
				tempI1 = temp * (k * r * cos(k * r) - sin(k * r));

				MulComp(Rejkr, Iejkr, -1, -k*r, tempR1, tempI1);
				tempR2 = tempR1;
				tempI2 = tempI1;
				tempR1 = tempR1 / r3 * x / 2 / Pi;
				tempI1 = tempI1 / r3 * x / 2 / Pi;
				tempR2 = tempR2 / r3 * y / 2 / Pi;
				tempI2 = tempI2 / r3 * y / 2 / Pi;

				double pi2wu0 = 2 * Pi * 2 * Pi * Mu0 * f; // �м�ֵ = 2*pi*w*u0 w = 2*pi*f

				MulComp(Rejkr, Iejkr, -3 * k * r, 3 - k * k * r * r, tempR1, tempI1);
				MulComp(sin(k * r), cos(k * r), 1, k * r, tempR2, tempI2);

				HHy_x0[i][j] = std::complex<double>((-tempR1 / r5 * (x * x + z * z) + tempR2 / r3 * 2) / pi2wu0,
					(-tempI1 / r5 * (x * x + z * z) + tempI2 / r3 * 2) / pi2wu0);

				HHy_y0[i][j] = std::complex<double>(-tempR1 / r5 * x * y / pi2wu0,
					-tempI1 / r5 * x * y / pi2wu0);

				PowerHHy_x0List[id] = PowerHHy_x0List[id] + ds * ds * (HHy_x0[i][j].real() * HHy_x0[i][j].real() + HHy_x0[i][j].imag() * HHy_x0[i][j].imag());
				PowerHHy_y0List[id] = PowerHHy_y0List[id] + ds * ds * (HHy_y0[i][j].real() * HHy_y0[i][j].real() + HHy_y0[i][j].imag() * HHy_y0[i][j].imag());
			}
		}//endloop

	}//openmp
	for (int i = 0; i < threadNum; i++) {
		PowerHHy_x0 = PowerHHy_x0List[i] + PowerHHy_x0;
		PowerHHy_y0 = PowerHHy_y0List[i] + PowerHHy_y0;
	}

	//ת��������
	FFT fft;
	fft.FFT_2D(HHy_x0, HHy_x0, N2, M2);
	fft.FFT_2D(HHy_y0, HHy_y0, N2, M2);

	//����������
	PowerFFTHHy_x0 = 0;
	PowerFFTHHy_y0 = 0;

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		double temppowH0 = pow((1.0 / (2 * N - 1) / ds), 2);//���ڼ���power����ʱ����
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				//Ϊ����һ��ѭ�� �ڴ˼���fft���Power
				PowerFFTHHy_x0List[id] = PowerFFTHHy_x0List[id] + temppowH0 * (HHy_x0[i][j].real() * HHy_x0[i][j].real() + HHy_x0[i][j].imag() * HHy_x0[i][j].imag());
				PowerFFTHHy_y0List[id] = PowerFFTHHy_y0List[id] + temppowH0 * (HHy_y0[i][j].real() * HHy_y0[i][j].real() + HHy_y0[i][j].imag() * HHy_y0[i][j].imag());
			}
		}
	}//omp
	for (int i = 0; i < threadNum; i++) {
		PowerFFTHHy_x0 = PowerFFTHHy_x0List[i] + PowerFFTHHy_x0;
		PowerFFTHHy_y0 = PowerFFTHHy_y0List[i] + PowerFFTHHy_y0;
	}
	//��������
	double temppowHHy_x0 = pow((PowerHHy_x0 / PowerFFTHHy_x0), 0.5);
	double temppowHHy_y0 = pow((PowerHHy_y0 / PowerFFTHHy_y0), 0.5);

	//�������-����
	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		double temppowH0 = pow((1.0 / (2 * N - 1) / ds), 2);//���ڼ���power����ʱ����
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				double tempR1, tempI1; //��ʽ��һ����
				double tempR2, tempI2; //��ʽ�ڶ�����
				//Hy1
				MulComp(Ex0[i][j].real(), Ex0[i][j].imag(), HHy_x0[i][j].real(), HHy_x0[i][j].imag(), tempR1, tempI1);
				MulComp(Ey0[i][j].real(), Ey0[i][j].imag(), HHy_y0[i][j].real(), HHy_y0[i][j].imag(), tempR2, tempI2);
				tempR1 = tempR1 * temppowHHy_x0;
				tempI1 = tempI1 * temppowHHy_x0;
				tempR2 = tempR2 * temppowHHy_y0;
				tempI2 = tempI2 * temppowHHy_y0;
				tempHy1[i][j] = std::complex<double>(tempR1 + tempR2, tempI1 + tempI2);
			}
		}
	}//omp

	//��任-��
	fft.IFFT_2D(tempHy1, tempHy1, N2, M2);

	//���
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Hy1[i][j] = tempHy1[i + N - 1][j + M - 1];
		}
	}

	//�����ڴ�
	Free_2D(tempHy1);	Free_2D(HHy_x0);	Free_2D(HHy_y0);
}

void FFTDI::FFTToHz1() {
	//׼����ʱ�����ʹ��ݺ�������
	tempHz1 = Allocate_2D(tempHz1, N2, M2);//���������
	HHz_x0 = Allocate_2D(HHz_x0, N2, M2);	   //���ݺ���
	HHz_y0 = Allocate_2D(HHz_y0, N2, M2);
	for (int i = 0; i < N2; i++)
		for (int j = 0; j < M2; j++) {
			tempHz1[i][j] = (0, 0);
			HHz_x0[i][j] = (0, 0);
			HHz_y0[i][j] = (0, 0);
		}

	PowerHHz_x0 = 0;	PowerHHz_y0 = 0;
	//���㴫�ݺ���
	vector<int> starti;
	vector<int> numi;
	vector<int> stopi;
	starti.resize(threadNum);	numi.resize(threadNum);	stopi.resize(threadNum);
	for (int i = 0; i < threadNum; i++) {
		int num = N2 / threadNum;
		int rr = (N2) % threadNum;
		if (i<rr && rr != 0) num += 1;
		numi[i] = num;
		//starti
		if (i == 0) starti[i] = 0;
		else starti[i] = starti[i - 1] + numi[i - 1];
		//stopi
		stopi[i] = starti[i] + numi[i];
	}
	vector<double>PowerHHz_x0List; PowerHHz_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHHz_x0List[i] = 0; }
	vector<double>PowerHHz_y0List; PowerHHz_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerHHz_y0List[i] = 0; }
	vector<double>PowerFFTHHz_x0List; PowerFFTHHz_x0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHHz_x0List[i] = 0; }
	vector<double>PowerFFTHHz_y0List; PowerFFTHHz_y0List.resize(threadNum);
	for (int i = 0; i < threadNum; i++) { PowerFFTHHz_y0List[i] = 0; }


	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				double x, y, z; //���޸�
				x = (i - N + 1) * ds + xp;
				//x = (i - N + 1) * ds;
				y = (j - N + 1) * ds + yp;
				z = zp;
				//z = z0 ;

				//�м����
				//double k1, k2; //ϵ��
				double Rejkr, Iejkr; //exp(-jkr)
				double tempR1, tempI1;
				double tempR2, tempI2;
				double r = pow((x * x + y * y + z * z), 0.5); //����
				double r5 = pow(r, 5);
				double r3 = pow(r, 3);

				Rejkr = cos(k * r);
				Iejkr = -sin(k * r);

				MulComp(Rejkr, Iejkr, 1, k*r, tempR1, tempI1);
				tempR1 = tempR1 / r3 * z / 2 / Pi;
				tempI1 = tempI1 / r3 * z / 2 / Pi;

				double temp = z / (2 * Pi * r * r * r);
				tempR1 = temp * (cos(k * r) + k * r * sin(k * r));
				tempI1 = temp * (k * r * cos(k * r) - sin(k * r));

				MulComp(Rejkr, Iejkr, -1, -k*r, tempR1, tempI1);
				tempR2 = tempR1;
				tempI2 = tempI1;
				tempR1 = tempR1 / r3 * x / 2 / Pi;
				tempI1 = tempI1 / r3 * x / 2 / Pi;
				tempR2 = tempR2 / r3 * y / 2 / Pi;
				tempI2 = tempI2 / r3 * y / 2 / Pi;

				double pi2wu0 = 2 * Pi * 2 * Pi * Mu0 * f; // �м�ֵ = 2*pi*w*u0 w = 2*pi*f

				MulComp(Rejkr, Iejkr, -3 * k * r, 3 - k * k * r * r, tempR1, tempI1);
				MulComp(sin(k * r), cos(k * r), 1, k * r, tempR2, tempI2);

				HHz_x0[i][j] = std::complex<double>(tempR1 / r5 * z * y / pi2wu0,
					tempI1 / r5 * z * y / pi2wu0);

				HHz_y0[i][j] = std::complex<double>(-tempR1 / r5 * x * z / pi2wu0,
					-tempI1 / r5 * x * z / pi2wu0);

				PowerHHz_x0List[id] = PowerHHz_x0List[id] + ds * ds * (HHz_x0[i][j].real() * HHz_x0[i][j].real() + HHz_x0[i][j].imag() * HHz_x0[i][j].imag());
				PowerHHz_y0List[id] = PowerHHz_y0List[id] + ds * ds * (HHz_y0[i][j].real() * HHz_y0[i][j].real() + HHz_y0[i][j].imag() * HHz_y0[i][j].imag());

			}
		}//endloop

	}//openmp
	for (int i = 0; i < threadNum; i++) {
		PowerHHz_x0 = PowerHHz_x0List[i] + PowerHHz_x0;
		PowerHHz_y0 = PowerHHz_y0List[i] + PowerHHz_y0;
	}

	//ת��������
	FFT fft;
	fft.FFT_2D(HHz_x0, HHz_x0, N2, M2);
	fft.FFT_2D(HHz_y0, HHz_y0, N2, M2);

	//����������
	PowerFFTHHz_x0 = 0;
	PowerFFTHHz_y0 = 0;

	omp_set_num_threads(threadNum);
	#pragma omp parallel
	{
		double temppowH0 = pow((1.0 / (2 * N - 1) / ds), 2);//���ڼ���power����ʱ����
		int id = omp_get_thread_num();
		for (int i = starti[id]; i < stopi[id]; i++)
		{
			for (int j = 0; j < M2; j++)
			{
				//Ϊ����һ��ѭ�� �ڴ˼���fft���Power
				PowerFFTHHz_x0List[id] = PowerFFTHHz_x0List[id] + temppowH0 * (HHz_x0[i][j].real() * HHz_x0[i][j].real() + HHz_x0[i][j].imag() * HHz_x0[i][j].imag());
				PowerFFTHHz_y0List[id] = PowerFFTHHz_y0List[id] + temppowH0 * (HHz_y0[i][j].real() * HHz_y0[i][j].real() + HHz_y0[i][j].imag() * HHz_y0[i][j].imag());

			}
		}
	}
	for (int i = 0; i < threadNum; i++) {
		PowerFFTHHz_x0 = PowerFFTHHz_x0List[i] + PowerFFTHHz_x0;
		PowerFFTHHz_y0 = PowerFFTHHz_y0List[i] + PowerFFTHHz_y0;
	}

	//��������
	double temppowHHz_x0 = pow((PowerHHz_x0 / PowerFFTHHz_x0), 0.5);
	double temppowHHz_y0 = pow((PowerHHz_y0 / PowerFFTHHz_y0), 0.5);

	//�������-����
	for (int i = 0; i < N2; i++)
	{
		for (int j = 0; j < M2; j++)
		{
			double tempR1, tempI1; //��ʽ��һ����
			double tempR2, tempI2; //��ʽ�ڶ�����
			//Hz1
			MulComp(Ex0[i][j].real(), Ex0[i][j].imag(), HHz_x0[i][j].real(), HHz_x0[i][j].imag(), tempR1, tempI1);
			MulComp(Ey0[i][j].real(), Ey0[i][j].imag(), HHz_y0[i][j].real(), HHz_y0[i][j].imag(), tempR2, tempI2);
			tempR1 = tempR1 * temppowHHz_x0;
			tempI1 = tempI1 * temppowHHz_x0;
			tempR2 = tempR2 * temppowHHz_y0;
			tempI2 = tempI2 * temppowHHz_y0;
			tempHz1[i][j] = std::complex<double>(tempR1 + tempR2, tempI1 + tempI2);
		}
	}

	//��任-��
	fft.IFFT_2D(tempHz1, tempHz1, N2, M2);

	//���
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Hz1[i][j] = tempHz1[i + N - 1][j + M - 1];
		}
	}

	//�����ڴ�
	Free_2D(tempHz1);	Free_2D(HHz_x0);	Free_2D(HHz_y0);
}

void FFTDI::StartCal()
{
	SourceFFT();
	FFTToExy1();
	FFTToEz1();
	FFTToHx1();
	FFTToHy1();
	FFTToHz1();
}

void FFTDI::output(std::complex<double>** Ex, std::complex<double>** Ey, std::complex<double>** Ez, std::complex<double>** Hx, std::complex<double>** Hy, std::complex<double>** Hz)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Ex[i][j] = Ex1[i][j];
			Ey[i][j] = Ey1[i][j];
			Ez[i][j] = Ez1[i][j];
			Hx[i][j] = Hx1[i][j];
			Hy[i][j] = Hy1[i][j];
			Hz[i][j] = Hz1[i][j];
		}
	}
	FreeCal();
}

void FFTDI::outsource(std::complex<double>** Ex, std::complex<double>** Ey)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Ex[i][j] = Ex1[i][j];
			Ey[i][j] = Ey1[i][j];
		}
	}
	FreeCal();
}

void FFTDI::MulComp(const double r1, const double i1, const double r2, const double i2, double & outR, double & outI)
{
	outR = r1 * r2 - i1 * i2;
	outI = r1 * i2 + r2 * i1;
}

std::complex<double> FFTDI::MulComp(std::complex<double> & In, double r2, double i2)
{
	return std::complex<double>(In.real() * r2 - In.imag() * i2,
		In.real() * i2 + r2 * In.imag());
}

double FFTDI::InserVal(const double x0, const double yUp, const double y0, const double yDown)
{
	if ((yUp == y0) || (yDown == y0))
	{
		return y0;
	}
	else
	{
		//����κ���aX^2+bX+c=Y
		double b1 = int(x0) - 1, a1 = b1 * b1;
		double b2 = int(x0), a2 = b2 * b2;
		double b3 = int(x0) + 1, a3 = b3 * b3;
		double e1 = a2 - a1, f1 = b2 - b1, g1 = y0 - yUp;
		double e2 = a3 - a1, f2 = b3 - b1, g2 = yDown - yUp;
		double y = (g2 * e1 - g1 * e2) / (f2 * e1 - f1 * e2);
		double x = (g1 - f1 * y) / e1;
		double z = yUp - b1 * y - a1 * x;
		return x * x0 * x0 + y * x0 + z;
	}
}



