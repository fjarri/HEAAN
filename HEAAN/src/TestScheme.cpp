/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "TestScheme.h"

#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>

#include "Ciphertext.h"
#include "EvaluatorUtils.h"
#include "Ring.h"
#include "Scheme.h"
#include "SchemeAlgo.h"
#include "SecretKey.h"
#include "StringUtils.h"
#include "TimeUtils.h"
#include "SerializationUtils.h"

using namespace std;
using namespace NTL;


ZZ mysum2(ZZ* arr, long n)
{
	ZZ a = to_ZZ(0);
	for (int i = 0; i < n; i++)
	{
		a += arr[i];
	}
	return a;
}

uint64_t mysum2(uint64_t* arr, long n)
{
	uint64_t a = 0;
	for (int i = 0; i < n; i++)
	{
		a += arr[i];
	}
	return a;
}


//----------------------------------------------------------------------------------
//   STANDARD TESTS
//----------------------------------------------------------------------------------


void TestScheme::testEncrypt(long logq, long logp, long logn) {
	cout << "!!! START TEST ENCRYPT !!!" << endl;
	srand(time(NULL));
	SetNumThreads(1);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher;

	timeutils.start("Encrypt");
	scheme.encrypt(cipher, mvec, n, logp, logq);
	timeutils.stop("Encrypt");

	timeutils.start("Decrypt");
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
	timeutils.stop("Decrypt");

	//StringUtils::compare(mvec, dvec, n, "val");

	for (int i = 0; i < n; i++)
	{
		cout << mvec[i].real() << " " << dvec[i].real() << endl;
		cout << mvec[i].imag() << " " << dvec[i].imag() << endl;
	}

	cout << "!!! END TEST ENCRYPT !!!" << endl;
}

void TestScheme::testEncryptSingle(long logq, long logp) {
	cout << "!!! START TEST ENCRYPT SINGLE !!!" << endl;
	srand(time(NULL));
	SetNumThreads(1);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	complex<double> mval = EvaluatorUtils::randomComplex();
	Ciphertext cipher;

	timeutils.start("Encrypt Single");
	scheme.encryptSingle(cipher, mval, logp, logq);
	timeutils.stop("Encrypt Single");

	complex<double> dval = scheme.decryptSingle(secretKey, cipher);

	StringUtils::compare(mval, dval, "val");

	cout << "!!! END TEST ENCRYPT SINGLE !!!" << endl;
}

void TestScheme::testAdd(long logq, long logp, long logn) {
	cout << "!!! START TEST ADD !!!" << endl;

	srand(time(NULL));
	SetNumThreads(1);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* madd = new complex<double>[n];

	for(long i = 0; i < n; i++) {
		madd[i] = mvec1[i] + mvec2[i];
	}

	Ciphertext cipher1, cipher2;
	scheme.encrypt(cipher1, mvec1, n, logp, logq);
	scheme.encrypt(cipher2, mvec2, n, logp, logq);

	timeutils.start("Addition");
	scheme.addAndEqual(cipher1, cipher2);
	timeutils.stop("Addition");

	complex<double>* dadd = scheme.decrypt(secretKey, cipher1);

	//StringUtils::compare(madd, dadd, n, "add");
	for (int i = 0; i < n; i++)
	{
		cout << madd[i].real() << " " << dadd[i].real() << endl;
		cout << madd[i].imag() << " " << dadd[i].imag() << endl;
	}

	cout << "!!! END TEST ADD !!!" << endl;
}

void TestScheme::testMult(long logq, long logp, long logn) {
	cout << "!!! START TEST MULT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(1);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmult = new complex<double>[n];
	for(long i = 0; i < n; i++) {
		mmult[i] = mvec1[i] * mvec2[i];
	}

	Ciphertext cipher1, cipher2;
	scheme.encrypt(cipher1, mvec1, n, logp, logq);
	scheme.encrypt(cipher2, mvec2, n, logp, logq);

	timeutils.start("Multiplication");
	scheme.multAndEqual(cipher1, cipher2);
	timeutils.stop("Multiplication");

	complex<double>* dmult = scheme.decrypt(secretKey, cipher1);

	//StringUtils::compare(mmult, dmult, n, "mult");
	for (int i = 0; i < n; i++)
	{
		cout << mmult[i].real() << " " << dmult[i].real() << endl;
		cout << mmult[i].imag() << " " << dmult[i].imag() << endl;
	}

	cout << "!!! END TEST MULT !!!" << endl;
}

void TestScheme::testimult(long logq, long logp, long logn) {
	cout << "!!! START TEST i MULTIPLICATION !!!" << endl;

	srand(time(NULL));
	SetNumThreads(1);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* imvec = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		imvec[i].real(-mvec[i].imag());
		imvec[i].imag(mvec[i].real());
	}

	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Multiplication by i");
	scheme.imultAndEqual(cipher);
	timeutils.stop("Multiplication by i");

	complex<double>* idvec = scheme.decrypt(secretKey, cipher);

	//StringUtils::compare(imvec, idvec, n, "imult");
	for (int i = 0; i < n; i++)
	{
		cout << imvec[i].real() << " " << idvec[i].real() << endl;
		cout << imvec[i].imag() << " " << idvec[i].imag() << endl;
	}

	cout << "!!! END TEST i MULTIPLICATION !!!" << endl;
}


//----------------------------------------------------------------------------------
//   ROTATE & CONJUGATE
//----------------------------------------------------------------------------------


void TestScheme::testRotateFast(long logq, long logp, long logn, long logr) {
	cout << "!!! START TEST ROTATE FAST !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	long r = (1 << logr);
	scheme.addLeftRotKey(secretKey, r);
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Left Rotate Fast");
	scheme.leftRotateFastAndEqual(cipher, r);
	timeutils.stop("Left Rotate Fast");

	complex<double>* dvec = scheme.decrypt(secretKey, cipher);

	EvaluatorUtils::leftRotateAndEqual(mvec, n, r);
	//StringUtils::compare(mvec, dvec, n, "rot");

	for (int i = 0; i < n; i++)
	{
		cout << mvec[i].real() << " " << dvec[i].real() << endl;
		cout << mvec[i].imag() << " " << dvec[i].imag() << endl;
	}

	cout << "!!! END TEST ROTATE BY POWER OF 2 BATCH !!!" << endl;
}

void TestScheme::testConjugate(long logq, long logp, long logn) {
	cout << "!!! START TEST CONJUGATE !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addConjKey(secretKey);

	long n = (1 << logn);

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mvecconj = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mvecconj[i] = conj(mvec[i]);
	}

	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Conjugate");
	scheme.conjugateAndEqual(cipher);
	timeutils.stop("Conjugate");

	complex<double>* dvecconj = scheme.decrypt(secretKey, cipher);
	//StringUtils::compare(mvecconj, dvecconj, n, "conj");

	for (int i = 0; i < n; i++)
	{
		cout << mvecconj[i].real() << " " << dvecconj[i].real() << endl;
		cout << mvecconj[i].imag() << " " << dvecconj[i].imag() << endl;
	}

	cout << "!!! END TEST CONJUGATE !!!" << endl;
}


//----------------------------------------------------------------------------------
//   POWER & PRODUCT TESTS
//----------------------------------------------------------------------------------


void TestScheme::testPowerOf2(long logq, long logp, long logn, long logdeg) {
	cout << "!!! START TEST POWER OF 2 !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	long degree = 1 << logdeg;
	complex<double>* mvec = new complex<double>[n];
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mvec[i] = EvaluatorUtils::randomCircle();

		//mpow[i] = pow(mvec[i], degree);
		mpow[i] = mvec[i];
		for (int j = 0; j < degree-1; j++)
		{
			mpow[i] = mpow[i] * mvec[i];
		}
	}

	Ciphertext cipher, cpow;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Power of 2");
	algo.powerOf2(cpow, cipher, logp, logdeg);
	timeutils.stop("Power of 2");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	//StringUtils::compare(mpow, dpow, n, "pow2");

	for (int i = 0; i < n; i++)
	{
		cout << mpow[i].real() << " " << dpow[i].real() << endl;
		cout << mpow[i].imag() << " " << dpow[i].imag() << endl;
	}

	cout << "!!! END TEST POWER OF 2 !!!" << endl;
}

//-----------------------------------------

void TestScheme::testPower(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST POWER !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomCircleArray(n);
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		//mpow[i] = pow(mvec[i], degree);
		mpow[i] = mvec[i];
		for (int j = 0; j < degree-1; j++)
		{
			mpow[i] = mpow[i] * mvec[i];
		}
	}

	Ciphertext cipher, cpow;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Power");
	algo.power(cpow, cipher, logp, degree);
	timeutils.stop("Power");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	//StringUtils::compare(mpow, dpow, n, "pow");

	for (int i = 0; i < n; i++)
	{
		cout << mpow[i].real() << " " << dpow[i].real() << endl;
		cout << mpow[i].imag() << " " << dpow[i].imag() << endl;
	}

	cout << "!!! END TEST POWER !!!" << endl;
}


//----------------------------------------------------------------------------------
//   FUNCTION TESTS
//----------------------------------------------------------------------------------


void TestScheme::testInverse(long logq, long logp, long logn, long steps) {
	cout << "!!! START TEST INVERSE !!!" << endl;

	srand(time(NULL));
	//SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomCircleArray(n, 0.1);
	complex<double>* minv = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		minv[i] = 1. / mvec[i];
	}

	Ciphertext cipher, cinv;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Inverse");
	algo.inverse(cinv, cipher, logp, steps);
	timeutils.stop("Inverse");

	complex<double>* dinv = scheme.decrypt(secretKey, cinv);
	//StringUtils::compare(minv, dinv, n, "inv");
	for (int i = 0; i < n; i++)
	{
		cout << minv[i].real() << " " << dinv[i].real() << endl;
		cout << minv[i].imag() << " " << dinv[i].imag() << endl;
	}

	cout << "!!! END TEST INVERSE !!!" << endl;
}

void TestScheme::testLogarithm(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST LOGARITHM !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n, 0.1);
	complex<double>* mlog = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mlog[i] = log(mvec[i] + 1.);
	}

	Ciphertext cipher, clog;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(LOGARITHM);
	algo.function(clog, cipher, LOGARITHM, logp, degree);
	timeutils.stop(LOGARITHM);

	complex<double>* dlog = scheme.decrypt(secretKey, clog);
	//StringUtils::compare(mlog, dlog, n, LOGARITHM);
	for (int i = 0; i < n; i++)
	{
		cout << mlog[i].real() << " " << dlog[i].real() << endl;
		cout << mlog[i].imag() << " " << dlog[i].imag() << endl;
	}

	cout << "!!! END TEST LOGARITHM !!!" << endl;
}

void TestScheme::testExponent(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST EXPONENT !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mvec[i]);
	}

	Ciphertext cipher, cexp;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(EXPONENT);
	algo.function(cexp, cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT);

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	//StringUtils::compare(mexp, dexp, n, EXPONENT);
	for (int i = 0; i < n; i++)
	{
		cout << mexp[i].real() << " " << dexp[i].real() << endl;
		cout << mexp[i].imag() << " " << dexp[i].imag() << endl;
	}

	cout << "!!! END TEST EXPONENT !!!" << endl;
}

void TestScheme::testExponentLazy(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST EXPONENT LAZY !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mvec[i]);
	}
	Ciphertext cipher, cexp;
	scheme.encrypt(cipher, mvec, n, logp, logQ);

	timeutils.start(EXPONENT + " lazy");
	algo.functionLazy(cexp, cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT + " lazy");

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	StringUtils::compare(mexp, dexp, n, EXPONENT);

	cout << "!!! END TEST EXPONENT LAZY !!!" << endl;
}

//-----------------------------------------

void TestScheme::testSigmoid(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST SIGMOID !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mvec[i]) / (1. + exp(mvec[i]));
	}

	Ciphertext cipher, csig;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(SIGMOID);
	algo.function(csig, cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID);

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	//StringUtils::compare(msig, dsig, n, SIGMOID);
	for (int i = 0; i < n; i++)
	{
		cout << msig[i].real() << " " << dsig[i].real() << endl;
		cout << msig[i].imag() << " " << dsig[i].imag() << endl;
	}

	cout << "!!! END TEST SIGMOID !!!" << endl;
}

void TestScheme::testSigmoidLazy(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST SIGMOID LAZY !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mvec[i]) / (1. + exp(mvec[i]));
	}

	Ciphertext cipher, csig;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(SIGMOID + " lazy");
	algo.functionLazy(csig, cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID + " lazy");

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	//StringUtils::compare(msig, dsig, n, SIGMOID);
	for (int i = 0; i < n; i++)
	{
		cout << msig[i].real() << " " << dsig[i].real() << endl;
		cout << msig[i].imag() << " " << dsig[i].imag() << endl;
	}

	cout << "!!! END TEST SIGMOID LAZY !!!" << endl;
}


void TestScheme::testWriteAndRead(long logq, long logp, long logSlots) {
	cout << "!!! START TEST WRITE AND READ !!!" << endl;

	cout << "!!! END TEST WRITE AND READ !!!" << endl;
}



int compare(double x, double y)
{
	return floor(-log2(abs(x - y)));
}


void statistics(complex<double>* mvec, complex<double> *dvec, int slots)
{
	int sum_bits = 0;
	int min_bits = 60;
	int max_bits = 0;

	for (int i = 0; i < slots; i++)
	{
		int b = compare(mvec[i].real(), dvec[i].real());
		sum_bits += b;
		if (b < min_bits) min_bits = b;
		if (b > max_bits) max_bits = b;

		b = compare(mvec[i].imag(), dvec[i].imag());
		sum_bits += b;
		if (b < min_bits) min_bits = b;
		if (b > max_bits) max_bits = b;
	}

	double mean = (double)sum_bits / slots / 2;

	double sum_squares = 0;
	for (int i = 0; i < slots; i++)
	{
		sum_squares += pow(compare(mvec[i].real(), dvec[i].real()) - mean, 2);
		sum_squares += pow(compare(mvec[i].imag(), dvec[i].imag()) - mean, 2);
	}

	double std = sqrt(sum_squares / (2 * slots - 1));

	cout << scientific;
	cout << "min=" << min_bits << " max=" << max_bits << " mean=" << mean << " std=" << std << endl;
}


void TestScheme::testBootstrap(long logq, long logp, long logSlots, long logT) {
	cout << "!!! START TEST BOOTSTRAP !!!" << endl;

	srand(time(NULL));
	//SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	timeutils.start("Key generating");
	scheme.addBootKey(secretKey, logSlots, logq + 4);
	timeutils.stop("Key generated");

	long slots = (1 << logSlots);
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(slots);

	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, slots, logp, logq);

	cout << "cipher logq before: " << cipher.logq << endl;

	scheme.modDownToAndEqual(cipher, logq);
	scheme.normalizeAndEqual(cipher);
	cipher.logq = logQ;
	cipher.logp = logq + 4;

	Ciphertext rot;
	timeutils.start("SubSum");
	for (long i = logSlots; i < logNh; ++i) {
		scheme.leftRotateFast(rot, cipher, (1 << i));
		scheme.addAndEqual(cipher, rot);
	}
	scheme.divByPo2AndEqual(cipher, logNh);
	timeutils.stop("SubSum");

	timeutils.start("CoeffToSlot");
	scheme.coeffToSlotAndEqual(cipher);
	timeutils.stop("CoeffToSlot");

	timeutils.start("EvalExp");
	scheme.evalExpAndEqual(cipher, logT);
	timeutils.stop("EvalExp");

	timeutils.start("SlotToCoeff");
	scheme.slotToCoeffAndEqual(cipher);
	timeutils.stop("SlotToCoeff");

	cipher.logp = logp;
	cout << "cipher logq after: " << cipher.logq << endl;
	cout << cipher.logp << " " << cipher.logq << " " << cipher.n << endl;

	complex<double>* dvec = scheme.decrypt(secretKey, cipher);

	//StringUtils::compare(mvec, dvec, slots, "boot");
	for (int i = 0; i < slots; i++)
	{
		cout << mvec[i].real() << " " << dvec[i].real() << endl;
		cout << mvec[i].imag() << " " << dvec[i].imag() << endl;
	}

	statistics(mvec, dvec, slots);

	cout << "!!! END TEST BOOTSRTAP !!!" << endl;
}


void TestScheme::testBootstrapSingleReal(long logq, long logQ, long logp, long logT) {
	cout << "!!! START TEST BOOTSTRAP SINGLE REAL !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	timeutils.start("Key generating");
	scheme.addBootKey(secretKey, 0, logq + 4);
	timeutils.stop("Key generated");

	double mval = EvaluatorUtils::randomReal();

	Ciphertext cipher;
	scheme.encryptSingle(cipher, mval, logp, logq);

	cout << "cipher logq before: " << cipher.logq << endl;
	scheme.modDownToAndEqual(cipher, logq);
	scheme.normalizeAndEqual(cipher);
	cipher.logq = logQ;

	Ciphertext rot, cconj;
	timeutils.start("SubSum");
	for (long i = 0; i < logNh; ++i) {
		scheme.leftRotateFast(rot, cipher, 1 << i);
		scheme.addAndEqual(cipher, rot);
	}
	scheme.conjugate(cconj, cipher);
	scheme.addAndEqual(cipher, cconj);
	scheme.divByPo2AndEqual(cipher, logN);
	timeutils.stop("SubSum");

	timeutils.start("EvalExp");
	scheme.evalExpAndEqual(cipher, logT);
	timeutils.stop("EvalExp");

	cout << "cipher logq after: " << cipher.logq << endl;

	cipher.logp = logp;
	complex<double> dval = scheme.decryptSingle(secretKey, cipher);

	StringUtils::compare(mval, dval.real(), "boot");

	cout << "!!! END TEST BOOTSRTAP SINGLE REAL !!!" << endl;
}

void TestScheme::test() {
	long logq = 1200;
	long logp = 30;
	long logn = 2; // n = 4

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;

	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn); // n = 4
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher;

	timeutils.start("Encrypt");
	scheme.encrypt(cipher, mvec, n, logp, logq);
	timeutils.stop("Encrypt");

	timeutils.start("Decrypt");
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
	timeutils.stop("Decrypt");

	StringUtils::compare(mvec, dvec, n, "val");
}

