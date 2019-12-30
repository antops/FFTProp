/*
*	created by jinming/liyun
*   function 
*/
#pragma once
#include <vector>
#include "../common/common/include/coordinate.h"
#include "../common/common/include/field_base.h"

using Common::Coordinate;
using Common::FieldBase;

namespace Antops
{
	struct FFTPropOption
	{
		bool prop_first_flag = true;  // �Ƿ��ȴ�������ת
		Coordinate coor; //�������λ��
		int omp_num = 1;
	};

	class FFTProp
	{
	public:
		FFTProp(const FFTPropOption& opt);

		virtual ~FFTProp();

		int Calculate(const FieldBase& in, FieldBase& out);

	private:
		FFTPropOption opt_;
	};
}

