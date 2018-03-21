//=================================================================================================
/*!
//  \file blaze/math/simd/Multiplication.h
//  \brief Header file for the SIMD multiplication functionality
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

#ifndef _BLAZE_MATH_SIMD_MULTIPLICATION_H_
#define _BLAZE_MATH_SIMD_MULTIPLICATION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/BasicTypes.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  SIMD MULTIPLICATION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication of two vectors of 16-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE2 and AVX2.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator*( const simd_i16_t<T>& a, const simd_i16_t<T>& b ) noexcept
#if BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi16( (~a).value, (~b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mullo_epi16( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of 16-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE2 and AVX2.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const simd_uint16_t
   operator*( const simd_i16_t<T1>& a, const simd_i16_t<T2>& b ) noexcept
#if BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi16( (~a).value, (~b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mullo_epi16( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 16-bit signed integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side complex values to be scaled.
// \param b The right-hand side scalars.
// \return The result of the scaling operation.
//
// This operation is only available for SSE2 and AVX2.
*/
BLAZE_ALWAYS_INLINE const simd_cint16_t
   operator*( const simd_cint16_t& a, const simd_int16_t& b ) noexcept
#if BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi16( (~a).value, (~b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mullo_epi16( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 16-bit unsigned integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side complex values to be scaled.
// \param b The right-hand side scalars.
// \return The result of the scaling operation.
//
// This operation is only available for SSE2 and AVX2.
*/
BLAZE_ALWAYS_INLINE const simd_cuint16_t
   operator*( const simd_cuint16_t& a, const simd_uint16_t& b ) noexcept
#if BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi16( (~a).value, (~b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mullo_epi16( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 16-bit signed integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side scalars.
// \param b The right-hand side complex values to be scaled.
// \return The result of the scaling operation.
//
// This operation is only available for SSE2 and AVX2.
*/
BLAZE_ALWAYS_INLINE const simd_cint16_t
   operator*( const simd_int16_t& a, const simd_cint16_t& b ) noexcept
#if BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi16( (~a).value, (~b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mullo_epi16( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 16-bit unsigned integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side scalars.
// \param b The right-hand side complex values to be scaled.
// \return The result of the scaling operation.
//
// This operation is only available for SSE2 and AVX2.
*/
BLAZE_ALWAYS_INLINE const simd_cuint16_t
   operator*( const simd_uint16_t& a, const simd_cuint16_t& b ) noexcept
#if BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi16( (~a).value, (~b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mullo_epi16( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of 16-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE2 and AVX2.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator*( const simd_ci16_t<T>& a, const simd_ci16_t<T>& b ) noexcept
#if BLAZE_AVX2_MODE
{
   __m256i x, y, z;
   const __m256i neg( _mm256_set_epi16( 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1 ) );

   x = _mm256_shufflelo_epi16( (~a).value, 0xA0 );
   x = _mm256_shufflehi_epi16( x, 0xA0 );
   z = _mm256_mullo_epi16( x, (~b).value );
   x = _mm256_shufflelo_epi16( (~a).value, 0xF5 );
   x = _mm256_shufflehi_epi16( x, 0xF5 );
   y = _mm256_shufflelo_epi16( (~b).value, 0xB1 );
   y = _mm256_shufflehi_epi16( y, 0xB1 );
   y = _mm256_mullo_epi16( x, y );
   y = _mm256_mullo_epi16( y, neg );
   return _mm256_add_epi16( z, y );
}
#elif BLAZE_SSE2_MODE
{
   __m128i x, y, z;
   const __m128i neg( _mm_set_epi16( 1, -1, 1, -1, 1, -1, 1, -1 ) );

   x = _mm_shufflelo_epi16( (~a).value, 0xA0 );
   x = _mm_shufflehi_epi16( x, 0xA0 );
   z = _mm_mullo_epi16( x, (~b).value );
   x = _mm_shufflelo_epi16( (~a).value, 0xF5 );
   x = _mm_shufflehi_epi16( x, 0xF5 );
   y = _mm_shufflelo_epi16( (~b).value, 0xB1 );
   y = _mm_shufflehi_epi16( y, 0xB1 );
   y = _mm_mullo_epi16( x, y );
   y = _mm_mullo_epi16( y, neg );
   return _mm_add_epi16( z, y );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of 32-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator*( const simd_i32_t<T>& a, const simd_i32_t<T>& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_mullo_epi32( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of 32-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const simd_uint32_t
   operator*( const simd_i32_t<T1>& a, const simd_i32_t<T2>& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_mullo_epi32( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 32-bit signed integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side complex values to be scaled.
// \param b The right-hand side scalars.
// \return The result of the scaling operation.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cint32_t
   operator*( const simd_cint32_t& a, const simd_int32_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_mullo_epi32( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 32-bit unsigned integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side complex values to be scaled.
// \param b The right-hand side scalars.
// \return The result of the scaling operation.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cuint32_t
   operator*( const simd_cuint32_t& a, const simd_uint32_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_mullo_epi32( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 32-bit signed integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side scalars.
// \param b The right-hand side complex values to be scaled.
// \return The result of the scaling operation.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const simd_cint32_t
   operator*( const simd_int32_t& a, const simd_cint32_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_mullo_epi32( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of 32-bit unsigned integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side scalars.
// \param b The right-hand side complex values to be scaled.
// \return The result of the scaling operation.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const simd_cuint32_t
   operator*( const simd_uint32_t& a, const simd_cuint32_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi32( (~a).value, (~b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_mullo_epi32( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of 32-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator*( const simd_ci32_t<T>& a, const simd_ci32_t<T>& b ) noexcept
#if BLAZE_MIC_MODE
{
   __m512i x, y, z;
   const __m512i neg( _mm256_set_epi32( 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1 ) );

   x = _mm512_shuffle_epi32( (~a).value, 0xA0 );
   z = _mm512_mullo_epi32( x, (~b).value );
   x = _mm512_shuffle_epi32( (~a).value, 0xF5 );
   y = _mm512_shuffle_epi32( (~b).value, 0xB1 );
   y = _mm512_mullo_epi32( x, y );
   y = _mm512_mullo_epi32( y, neg );
   return _mm512_add_epi32( z, y );
}
#elif BLAZE_AVX2_MODE
{
   __m256i x, y, z;
   const __m256i neg( _mm256_set_epi32( 1, -1, 1, -1, 1, -1, 1, -1 ) );

   x = _mm256_shuffle_epi32( (~a).value, 0xA0 );
   z = _mm256_mullo_epi32( x, (~b).value );
   x = _mm256_shuffle_epi32( (~a).value, 0xF5 );
   y = _mm256_shuffle_epi32( (~b).value, 0xB1 );
   y = _mm256_mullo_epi32( x, y );
   y = _mm256_mullo_epi32( y, neg );
   return _mm256_add_epi32( z, y );
}
#elif BLAZE_SSE4_MODE
{
   __m128i x, y, z;
   const __m128i neg( _mm_set_epi32( 1, -1, 1, -1 ) );

   x = _mm_shuffle_epi32( (~a).value, 0xA0 );
   z = _mm_mullo_epi32( x, (~b).value );
   x = _mm_shuffle_epi32( (~a).value, 0xF5 );
   y = _mm_shuffle_epi32( (~b).value, 0xB1 );
   y = _mm_mullo_epi32( x, y );
   y = _mm_mullo_epi32( y, neg );
   return _mm_add_epi32( z, y );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of 64-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for AVX-512.
*/
template< typename T >  // Type of the left-hand side operand
BLAZE_ALWAYS_INLINE const T
   operator*( const simd_i64_t<T>& a, const simd_i64_t<T>& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi64( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of 64-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const simd_uint64_t
   operator*( const simd_i64_t<T1>& a, const simd_i64_t<T2>& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi64( (~a).value, (~b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of single precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_float_t
   operator*( const simd_float_t& a, const simd_float_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
{
   return _mm_mul_ps( a.value, b.value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of single precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side complex values to be scaled.
// \param b The right-hand side scalars.
// \return The result of the scaling operation.
//
// This operation is only available for SSE, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cfloat_t
   operator*( const simd_cfloat_t& a, const simd_float_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
{
   return _mm_mul_ps( a.value, b.value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of single precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side scalars.
// \param b The right-hand side complex values to be scaled.
// \return The result of the scaling operation.
//
// This operation is only available for SSE, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cfloat_t
   operator*( const simd_float_t& a, const simd_cfloat_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
{
   return _mm_mul_ps( a.value, b.value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of single precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE3 and AVX.
*/
BLAZE_ALWAYS_INLINE const simd_cfloat_t
   operator*( const simd_cfloat_t& a, const simd_cfloat_t& b ) noexcept
#if BLAZE_AVX_MODE
{
   __m256 x, y, z;

   x = _mm256_shuffle_ps( a.value, a.value, 0xA0 );
   z = _mm256_mul_ps( x, b.value );
   x = _mm256_shuffle_ps( a.value, a.value, 0xF5 );
   y = _mm256_shuffle_ps( b.value, b.value, 0xB1 );
   y = _mm256_mul_ps( x, y );
   return _mm256_addsub_ps( z, y );
}
#elif BLAZE_SSE3_MODE
{
   __m128 x, y, z;

   x = _mm_shuffle_ps( a.value, a.value, 0xA0 );
   z = _mm_mul_ps( x, b.value );
   x = _mm_shuffle_ps( a.value, a.value, 0xF5 );
   y = _mm_shuffle_ps( b.value, b.value, 0xB1 );
   y = _mm_mul_ps( x, y );
   return _mm_addsub_ps( z, y );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of double precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_double_t
   operator*( const simd_double_t& a, const simd_double_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mul_pd( a.value, b.value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of double precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side complex values to be scaled.
// \param b The right-hand side scalars.
// \return The result of the scaling operation.
//
// This operation is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cdouble_t
   operator*( const simd_cdouble_t& a, const simd_double_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mul_pd( a.value, b.value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of a vector of double precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side scalars.
// \param b The right-hand side complex values to be scaled.
// \return The result of the scaling operation.
//
// This operation is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cdouble_t
   operator*( const simd_double_t& a, const simd_cdouble_t& b ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mul_pd( a.value, b.value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication of two vectors of double precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the multiplication.
//
// This operation is only available for SSE3 and AVX.
*/
BLAZE_ALWAYS_INLINE const simd_cdouble_t
   operator*( const simd_cdouble_t& a, const simd_cdouble_t& b ) noexcept
#if BLAZE_AVX_MODE
{
   __m256d x, y, z;

   x = _mm256_shuffle_pd( a.value, a.value, 0 );
   z = _mm256_mul_pd( x, b.value );
   x = _mm256_shuffle_pd( a.value, a.value, 15 );
   y = _mm256_shuffle_pd( b.value, b.value, 5 );
   y = _mm256_mul_pd( x, y );
   return _mm256_addsub_pd( z, y );
}
#elif BLAZE_SSE3_MODE
{
   __m128d x, y, z;

   x = _mm_shuffle_pd( a.value, a.value, 0 );
   z = _mm_mul_pd( x, b.value );
   x = _mm_shuffle_pd( a.value, a.value, 3 );
   y = _mm_shuffle_pd( b.value, b.value, 1 );
   y = _mm_mul_pd( x, y );
   return _mm_addsub_pd( z, y );
}
#else
= delete;
#endif
//*************************************************************************************************

} // namespace blaze

#endif
