
/**************************************************************************************
* version 1.0 2016/11/7 modified 20190227
* Designer liyun & MingJin
* E-mail 465120258@qq.com jinmingaps@163.com
* Fuction�� 
* ����ƽ����䴫��
* ���� Ex Ey �õ�Ŀ�����Ex Ey Ez Hx Hy Hz
* ��Ҫ���������Ƶ�ʣ�Ĭ��10.65e9����Ŀ����롢�������N*M��
Ŀ������б�ǣ�x�᷽�򣩡�ds��Ĭ�ϲ���/3.5��
* version 1.2 2019/05/21 ʡ�ڴ��
* E-mail jinmingaps@163.com
* Modification:
* ��������ʡ�ڴ�
***************************************************************************************/
#pragma once
#ifndef FFTDI_H
#define FFTDI_H

#include <cmath>
#include <complex>
#include <omp.h>

namespace Antops
{
	class FFTDI
	{
	public:
		FFTDI(double f = 10.65e9, double xp = 0, double yp = 0, double zp = 1, int N = 361, int M = 361);
		~FFTDI();

		void Initialization();
		void FreeCal();
		void Setds(double ds1); //����ds
		void SetthreadNum(int _threadNum) { threadNum = _threadNum; }

		//�������벢��0
		void SetInput(std::complex<double> ** ExIn, std::complex<double> ** EyIn);

		void StartCal();

		void SourceFFT();
		void FFTToExy1();
		void FFTToEz1();
		void FFTToHx1();
		void FFTToHy1();
		void FFTToHz1();


		//���Ex Ey Ez Hx Hy Hz ������FreeCal
		void output(std::complex<double> ** Ex, std::complex<double> ** Ey, std::complex<double> ** Ez,
			std::complex<double> ** Hx, std::complex<double> ** Hy, std::complex<double> ** Hz);
		void outsource(std::complex<double> ** Ex, std::complex<double> ** Ey);
	private:
		//���帴���ĳ�
		void MulComp(const double r1, const double i1, const double r2, const double i2, double & outR, double & outI);
		std::complex<double> MulComp(std::complex<double> & In, double r2, double i2);
		//��ֵ����
		double InserVal(const double x0, const double yUp, const double y0, const double yDown);

		//�����볡���в�ֵ
		//void InserExEy();

	private:
		double f; // Ƶ�� Ĭ�� 10.65e9
		double lamda; // ����
		double k;

		int threadNum;

		//double theta;
		//double theta_h0;
		double ds;
		double xp;	//����λ��
		double yp;  //����λ��
		double zp; //Ŀ�����

		int N, M; //����ĵ��� һ����Ϊ���� Ĭ��361
		int N2, M2; // N2 = 2 * N -1

		//���ݺ���
		std::complex<double> ** HExy, ** HEz_x0, **HEz_y0;
		std::complex<double> ** HHx_x0, ** HHx_y0,
			** HHy_x0, ** HHy_y0,
			** HHz_x0, ** HHz_y0;

		//��0�������
		std::complex<double> ** Ex0, ** Ey0;

		// ��˹ά������
		double 	PowerHExy, PowerFFTHExy;
		double 	PowerHEz_x0, PowerFFTHEz_x0;
		double 	PowerHEz_y0, PowerFFTHEz_y0;
		double 	PowerHHx_x0, PowerFFTHHx_x0;
		double 	PowerHHx_y0, PowerFFTHHx_y0;
		double 	PowerHHy_x0, PowerFFTHHy_x0;
		double 	PowerHHy_y0, PowerFFTHHy_y0;
		double 	PowerHHz_x0, PowerFFTHHz_x0;
		double 	PowerHHz_y0, PowerFFTHHz_y0;

		//Ex1 Ey1 Ez1 Hx1 Hy1 Hz1
		std::complex<double> ** Ex1, ** Ey1, **Ez1;
		std::complex<double> ** Hx1, ** Hy1, **Hz1;
		std::complex<double> ** tempEx1, ** tempEy1, **tempEz1;
		std::complex<double> ** tempHx1, ** tempHy1, **tempHz1;
	};
}


#endif // !CALCUlATION_H
