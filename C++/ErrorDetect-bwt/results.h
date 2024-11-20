#include <iostream>
#include <functional>
#include <vector>

using namespace std;

class results{
public:
	results(const string &pid, double errval, double errabs) :PointID(pid), errVal(errval), errAbs(errabs){};
	string PointID;
	double errAbs;
	double errVal;
	/* �����ľ���ֵ��С����������������� */
	bool operator > (const results &res)const {
		return errAbs > res.errAbs;
	}
	/* ���������������������� */
	bool operator < (const results &res)const {
		return atoi(PointID.c_str()) < atoi(res.PointID.c_str());
	}

};