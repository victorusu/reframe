//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/unmqr.h
//  \brief Header file for the CLAPACK unmqr wrapper functions
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_UNMQR_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_UNMQR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
extern "C" {

void cunmqr_( char* side, char* trans, int* m, int* n, int* k, float*  A, int* lda, float*  tau, float*  C, int* ldc, float*  work, int* lwork, int* info );
void zunmqr_( char* side, char* trans, int* m, int* n, int* k, double* A, int* lda, double* tau, double* C, int* ldc, double* work, int* lwork, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK FUNCTIONS TO MULTIPLY Q FROM A QR DECOMPOSITION WITH A MATRIX (UNMQR)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK functions to multiply Q from a QR decomposition with a matrix (unmqr) */
//@{
inline void unmqr( char side, char trans, int m, int n, int k, const complex<float>* A, int lda,
                   const complex<float>* tau, complex<float>* C, int ldc, complex<float>* work,
                   int lwork, int* info );

inline void unmqr( char side, char trans, int m, int n, int k, const complex<double>* A, int lda,
                   const complex<double>* tau, complex<double>* C, int ldc, complex<double>* work,
                   int lwork, int* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the multiplication of the single precision Q from a QR decomposition
//        with another matrix.
// \ingroup lapack_decomposition
//
// \param side \c 'L' to apply \f$ Q \f$ or \f$ Q^H \f$ from the left, \c 'R' to apply from the right.
// \param trans \c 'N' for \f$ Q \f$, \c 'C' for \f$ Q^H \f$.
// \param m The number of rows of the matrix \a C \f$[0..\infty)\f$.
// \param n The number of columns of the matrix \a C \f$[0..\infty)\f$.
// \param k The number of elementary reflectors, whose product defines the matrix \a Q.
// \param A Pointer to the first element of the QR decomposed single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \a A \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param C Pointer to the first element of the single precision column-major matrix multiplicator.
// \param ldc The total number of elements between two columns of the matrix \a C \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work.
// \param info Return code of the function call.
// \return void
//
// This function multiplies the \a Q matrix resulting from the QR decomposition of the sgeqrf()
// function with the given general single precision \a m-by-\a n matrix \a C. Depending on the
// settings of \a side and \a trans it overwrites \a C with

   \code
                | side = 'L'    | side = 'R'
   -------------|---------------|---------------
   trans = 'N': | Q * C         | C * Q
   trans = 'C': | ctrans(Q) * C | C * ctrans(Q)
   \endcode

// In case \a side is specified as \c 'L', the \a Q matrix is expected to be of size \a m-by-\a m,
// in case \a side is set to \c 'R', \a Q is expected to be of size \a n-by-\a n.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the cunmqr() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void unmqr( char side, char trans, int m, int n, int k, const complex<float>* A, int lda,
                   const complex<float>* tau, complex<float>* C, int ldc, complex<float>* work,
                   int lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cunmqr_( &side, &trans, &m, &n, &k,
            const_cast<float*>( reinterpret_cast<const float*>( A ) ), &lda,
            const_cast<float*>( reinterpret_cast<const float*>( tau ) ),
            reinterpret_cast<float*>( C ), &ldc, reinterpret_cast<float*>( work ),
            &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the multiplication of the double precision Q from a QR decomposition
//        with another matrix.
// \ingroup lapack_decomposition
//
// \param side \c 'L' to apply \f$ Q \f$ or \f$ Q^H \f$ from the left, \c 'R' to apply from the right.
// \param trans \c 'N' for \f$ Q \f$, \c 'C' for \f$ Q^H \f$.
// \param m The number of rows of the matrix \a C \f$[0..\infty)\f$.
// \param n The number of columns of the matrix \a C \f$[0..\infty)\f$.
// \param k The number of elementary reflectors, whose product defines the matrix \a Q.
// \param A Pointer to the first element of the QR decomposed double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \a A \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param C Pointer to the first element of the double precision column-major matrix multiplicator.
// \param ldc The total number of elements between two columns of the matrix \a C \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work.
// \param info Return code of the function call.
// \return void
//
// This function multiplies the \a Q matrix resulting from the QR decomposition of the sgeqrf()
// function with the given general double precision \a m-by-\a n matrix \a C. Depending on the
// settings of \a side and \a trans it overwrites \a C with

   \code
                | side = 'L'    | side = 'R'
   -------------|---------------|---------------
   trans = 'N': | Q * C         | C * Q
   trans = 'C': | ctrans(Q) * C | C * ctrans(Q)
   \endcode

// In case \a side is specified as \c 'L', the \a Q matrix is expected to be of size \a m-by-\a m,
// in case \a side is set to \c 'R', \a Q is expected to be of size \a n-by-\a n.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the zunmqr() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void unmqr( char side, char trans, int m, int n, int k, const complex<double>* A, int lda,
                   const complex<double>* tau, complex<double>* C, int ldc, complex<double>* work,
                   int lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zunmqr_( &side, &trans, &m, &n, &k,
            const_cast<double*>( reinterpret_cast<const double*>( A ) ), &lda,
            const_cast<double*>( reinterpret_cast<const double*>( tau ) ),
            reinterpret_cast<double*>( C ), &ldc, reinterpret_cast<double*>( work ),
            &lwork, info );
}
//*************************************************************************************************

} // namespace blaze

#endif
