/**************************************************************************************
* version 1.0 2019/2/26
* Designer ����
* E-mail jinmingaps@163.com
* Fuction��
* ������ת��
* ���� Eu0, Ev0���õ�������ת���Eut,Evt;
* ��Ҫ����Ĳ�����Ƶ�ʣ�Ŀ����룬�������N*M;
Ŀ������б�ǣ�theta��phi����ds��Ĭ�ϲ���/3.5��
***************************************************************************************/
#pragma once
#ifndef FFTASRotation_H
#define FFTASRotation_H

#include <cmath>
#include <complex>
#include <vector>
#include "../common/common/include/field_base.h"
#include "../common/common/include/coordinate.h"
#include "FFTProp.h"
#include "FFT.h"
using namespace Common;
namespace Antops
{
	class FFTASRotation
	{
	public:
		FFTASRotation();
		~FFTASRotation();
		int Calculate(const FieldBase& in, const FFTPropOption& opt, FieldBase& out);
	
	private:
		double freq; // Ƶ�� Ĭ�� 10.65e9
		double lambda; // ����
		double k;
		double ds;
		double theta_rou, phi_rou;//��ת
		int Nu, Nv; //����ĵ��� һ����Ϊ���� Ĭ��361
		int Nu2, Nv2; // N2 = 2 * N -1

		std::complex<double> ** Eu0, ** Ev0, ** En0;
		std::complex<double> ** Eut, ** Evt, ** Ent;

	};

}

#endif // !CALCUlATION_H
