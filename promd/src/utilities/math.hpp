#ifndef MATH_HPP
#define MATH_HPP

#include "config.hpp"
// #include "stdout.hpp"


#define DIM 3UL

#define EPSILONLIMIT 1.0e-10

#if defined(SINGLE)
#define prec "float"
typedef float number;
#else
#define prec "double"
typedef double number;
#endif

// #define RVec number[3]
// typedef number RVec[DIM];
// typedef std::vector<number> vector;
// typedef std::vector<RVec> RVector;

#include "blaze/Blaze.h"

// typedef blaze::DynamicVector<number, true> RVec;
typedef blaze::StaticVector<number, DIM> RVec;
typedef blaze::DynamicVector<RVec, blaze::rowVector> RVector;

// typedef blaze::DynamicVector<int, true> IVec;
typedef blaze::StaticVector<int, DIM> IVec;

typedef blaze::StaticMatrix<number,DIM,DIM>  Tensor;

typedef blaze::DynamicVector<Tensor, blaze::rowVector> TVector;

typedef blaze::DynamicVector<number, blaze::rowVector> NVector;
typedef blaze::DynamicVector<int, blaze::rowVector> IVector;


// periodic boudary conditions
#define nint(x) (((x) > 0) ? (int)((x)+0.5) : (int)((x)-0.5))

#define XX    0
#define YY    1
#define ZZ    2

VHR_ALWAYS_INLINE
void pushBack(RVector & rvector, RVec &rvec)
{
	const int size = rvector.size();

	rvector.extend(1);
	rvector[size] = rvec;
}


VHR_ALWAYS_INLINE
void pushBack(NVector & nvector, number value)
{
	const int size = nvector.size();

	nvector.extend(1);
	nvector[size] = value;
}

VHR_ALWAYS_INLINE
void pushBack(IVector & ivector, int value)
{
	const int size = ivector.size();

	ivector.extend(1);
	ivector[size] = value;
}


VHR_ALWAYS_INLINE
void pushBack(RVector & rvector, const number x, const number y, const number z)
{
    RVec rvec;
    rvec[XX] = x;
    rvec[YY] = y;
    rvec[ZZ] = z;

	const int size = rvector.size();

	rvector.extend(1);
	rvector[size] = rvec;
}

VHR_ALWAYS_INLINE
number sum(RVec &rvec)
{
	int i = 0;
	number sum = 0.0;
	for(i = 0; i < rvec.size(); i++)
		sum += rvec[i];

	return sum;
}


VHR_ALWAYS_INLINE
number trace(Tensor &tensor)
{
	int i = 0;
	number sum = 0.0;
	for(i = 0; i < DIM; i++)
		sum += tensor(i, i);

	return sum;
}


VHR_ALWAYS_INLINE
RVec rvecNINT(RVec rvec)
{
	RVec r;
	int i;
	for(i = 0; i < DIM; i++)
		r[i] = nint(rvec[i]);

	return r;
}


VHR_ALWAYS_INLINE
number cosine(RVec &u, RVec &v)
{
	return (trans(u) * v) / (std::sqrt((trans(u) * u)) * std::sqrt((trans(v) * v)));
}

#endif //MATH_HPP
