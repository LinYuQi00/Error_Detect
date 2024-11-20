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
	/* 按误差的绝对值大小降序排序基本操作符 */
	bool operator > (const results &res)const {
		return errAbs > res.errAbs;
	}
	/* 按点号升序排序基本操作符 */
	bool operator < (const results &res)const {
		return atoi(PointID.c_str()) < atoi(res.PointID.c_str());
	}

};