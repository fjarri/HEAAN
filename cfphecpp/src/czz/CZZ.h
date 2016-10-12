#ifndef UTILS_CZZ_H_
#define UTILS_CZZ_H_

#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class CZZ {
public:
	ZZ r;
	ZZ i;

	CZZ(ZZ r, ZZ i) : r(r), i(i) {};
	CZZ() : r(ZZ::zero()), i(ZZ::zero()) {};

	CZZ operator+(const CZZ& o);
	CZZ operator+(const ZZ& m);
	void operator+=(const CZZ& o);

	CZZ operator-(const CZZ& o);
	CZZ operator-(const ZZ& m);
	void operator-=(const CZZ& o);

	CZZ operator *(const CZZ& o);
	CZZ operator *(const ZZ& o);
	CZZ operator *(const long& o);
	void operator *=(const CZZ& o);

	CZZ operator /(const ZZ& o);
	CZZ operator /(const long& o);

	CZZ operator %(const ZZ& o);
	void operator%=(const ZZ& o);

	CZZ operator <<(const long& s);
	CZZ operator >>(const long& s);

	CZZ sqr();
	void sqrThis();
	ZZ norm();

	string toString();
};

#endif /* UTILS_CZZ_H_ */