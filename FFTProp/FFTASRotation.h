/**************************************************************************************
* version 1.0 2019/2/26
* Designer 金铭
* E-mail jinmingaps@163.com
* Fuction：
* 坐标旋转：
* 输入 Eu0, Ev0，得到坐标旋转后的Eut,Evt;
* 需要输入的参数，频率，目标距离，计算点数N*M;
目标面倾斜角（theta和phi）、ds（默认波长/3.5）
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
		double freq; // 频率 默认 10.65e9
		double lambda; // 波长
		double k;
		double ds;
		double theta_rou, phi_rou;//旋转
		int Nu, Nv; //计算的点数 一般设为奇数 默认361
		int Nu2, Nv2; // N2 = 2 * N -1

		std::complex<double> ** Eu0, ** Ev0, ** En0;
		std::complex<double> ** Eut, ** Evt, ** Ent;

	};

}

#endif // !CALCUlATION_H
