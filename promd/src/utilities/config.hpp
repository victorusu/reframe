#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <string.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(USE_MPI)
#include <mpi.h>
#endif

#define VHR_RESTRICT __restrict__
// #define VHR_ALWAYS_INLINE __attribute__((always_inline)) inline
// #define VHR_INLINE inline

#if ( defined(_MSC_VER) || defined(__INTEL_COMPILER) )
#  define VHR_INLINE __forceinline
#else
#  define VHR_INLINE inline
#endif

#if defined(__GNUC__)
#  define VHR_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#  define VHR_ALWAYS_INLINE VHR_INLINE
#endif


#endif //CONFIG_HPP


