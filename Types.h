/**
*	Header for basic types for real-valued numbers. By default, this program uses float precision,
*	but users can change data types easily here.
*/

#pragma once
#include <complex>

// Do platform-specific check for x86 or x64 architecture.
#if _WIN32 || _WIN64
#if _WIN64
#define _USE64BIT
#else
#define _USE32BIT
#endif
#endif
#if __GNUC__
#if __x86_64__ || __ppc64__
#define _USE64BIT
#else
#define _USE32BIT
#endif
#endif

namespace MathParser
{
	// Define precision of real/complex numbers, change to 'double' type as necessary.
	typedef float real;
	typedef std::complex<real> complex;

	// Define high precision reals as max possible precision value.
	typedef long double highp;

	// Define highest precision int.
#ifdef _USE64BIT
	typedef int64_t highpInt;
#elif _USE32BIT
	typedef int32_t highpInt;
#endif

	// Define highest precision uint.
#ifdef _USE64BIT
	typedef uint64_t highpUint;
#elif _USE32BIT
	typedef uint32_t highpUint;
#endif

}