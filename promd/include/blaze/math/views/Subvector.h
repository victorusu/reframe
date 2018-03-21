//=================================================================================================
/*!
//  \file blaze/math/views/Subvector.h
//  \brief Header file for all restructuring subvector functions
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

#ifndef _BLAZE_MATH_VIEWS_SUBVECTOR_H_
#define _BLAZE_MATH_VIEWS_SUBVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/expressions/Vector.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsCrossExpr.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/typetraits/IsVecEvalExpr.h>
#include <blaze/math/typetraits/IsVecForEachExpr.h>
#include <blaze/math/typetraits/IsVecScalarDivExpr.h>
#include <blaze/math/typetraits/IsVecScalarMultExpr.h>
#include <blaze/math/typetraits/IsVecSerialExpr.h>
#include <blaze/math/typetraits/IsVecTransExpr.h>
#include <blaze/math/typetraits/IsVecVecAddExpr.h>
#include <blaze/math/typetraits/IsVecVecDivExpr.h>
#include <blaze/math/typetraits/IsVecVecMultExpr.h>
#include <blaze/math/typetraits/IsVecVecSubExpr.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup views
//
// \param vector The vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given vector.
// The following example demonstrates the creation of a dense and sparse subvector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   typedef blaze::DynamicVector<double,columnVector>  DenseVector;
   typedef blaze::CompressedVector<int,rowVector>     SparseVector;

   DenseVector  d;
   SparseVector s;
   // ... Resizing and initialization

   // Creating a dense subvector of size 8, starting from index 4
   blaze::DenseSubvector<DenseVector> dsv = subvector( d, 4UL, 8UL );

   // Creating a sparse subvector of size 7, starting from index 5
   blaze::SparseSubvector<SparseVector> ssv = subvector( s, 5UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is larger
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// Please note that this function creates an unaligned dense or sparse subvector. For instance,
// the creation of the dense subvector is equivalent to the following three function calls:

   \code
   blaze::DenseSubvector<DenseVector>           dsv = subvector<unaligned>( v, 4UL, 8UL );
   blaze::DenseSubvector<DenseVector,unaligned> dsv = subvector           ( v, 4UL, 8UL );
   blaze::DenseSubvector<DenseVector,unaligned> dsv = subvector<unaligned>( v, 4UL, 8UL );
   \endcode

// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions. However, especially in case of dense subvectors this may
// result in considerable performance improvements. In order to create an aligned subvector the
// following function call has to be used:

   \code
   blaze::DenseSubvector<DenseVector,aligned> = subvector<aligned>( v, 4UL, 8UL );
   \endcode

// Note however that in this case the given \a index and \a size are subject to additional checks
// to guarantee proper alignment.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline SubvectorExprTrait_<VT,unaligned>
   subvector( Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup views
//
// \param vector The vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given vector.
// The following example demonstrates the creation of a dense and sparse subvector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   typedef blaze::DynamicVector<double,columnVector>  DenseVector;
   typedef blaze::CompressedVector<int,rowVector>     SparseVector;

   DenseVector  d;
   SparseVector s;
   // ... Resizing and initialization

   // Creating a dense subvector of size 8, starting from index 4
   blaze::DenseSubvector<DenseVector> dsv = subvector( d, 4UL, 8UL );

   // Creating a sparse subvector of size 7, starting from index 5
   blaze::SparseSubvector<SparseVector> ssv = subvector( s, 5UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is larger
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// Please note that this function creates an unaligned dense or sparse subvector. For instance,
// the creation of the dense subvector is equivalent to the following three function calls:

   \code
   blaze::DenseSubvector<DenseVector>           dsv = subvector<unaligned>( v, 4UL, 8UL );
   blaze::DenseSubvector<DenseVector,unaligned> dsv = subvector           ( v, 4UL, 8UL );
   blaze::DenseSubvector<DenseVector,unaligned> dsv = subvector<unaligned>( v, 4UL, 8UL );
   \endcode

// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions. However, especially in case of dense subvectors this may
// result in considerable performance improvements. In order to create an aligned subvector the
// following function call has to be used:

   \code
   blaze::DenseSubvector<DenseVector,aligned> = subvector<aligned>( v, 4UL, 8UL );
   \endcode

// Note however that in this case the given \a index and \a size are subject to additional checks
// to guarantee proper alignment.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline SubvectorExprTrait_<const VT,unaligned>
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup views
//
// \param vector The vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given dense or sparse vector, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned subvector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   typedef blaze::DynamicVector<double,columnVector>  DenseVector;
   typedef blaze::CompressedVector<int,rowVector>     SparseVector;

   DenseVector  d;
   SparseVector s;
   // ... Resizing and initialization

   // Creating an aligned dense subvector of size 8 starting from index 4
   blaze::DenseSubvector<DenseVector,aligned> dsv = subvector<aligned>( d, 4UL, 8UL );

   // Creating an unaligned subvector of size 7 starting from index 3
   blaze::SparseSubvector<SparseVector,unaligned> ssv = subvector<unaligned>( s, 3UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is larger
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions and the given \a index is subject to additional checks to
// guarantee proper alignment. However, especially in case of dense subvectors this may result
// in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of the subvector must be aligned. The following source code gives some examples
// for a double precision dynamic vector, assuming that AVX is available, which packs 4 \c double
// values into a SIMD vector:

   \code
   using blaze::columnVector;

   typedef blaze::DynamicVector<double,columnVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType,aligned>  SubvectorType;

   VectorType d( 17UL );
   // ... Resizing and initialization

   // OK: Starts at the beginning, i.e. the first element is aligned
   SubvectorType dsv1 = subvector<aligned>( d, 0UL, 13UL );

   // OK: Start index is a multiple of 4, i.e. the first element is aligned
   SubvectorType dsv2 = subvector<aligned>( d, 4UL, 7UL );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   SubvectorType dsv3 = subvector<aligned>( d, 8UL, 9UL );

   // Error: Start index is not a multiple of 4, i.e. the first element is not aligned
   SubvectorType dsv4 = subvector<aligned>( d, 5UL, 8UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DisableIf_< Or< IsComputation<VT>, IsTransExpr<VT> >, SubvectorExprTrait_<VT,AF> >
   subvector( Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return SubvectorExprTrait_<VT,AF>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup views
//
// \param vector The vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given dense or sparse vector, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned subvector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   typedef blaze::DynamicVector<double,columnVector>  DenseVector;
   typedef blaze::CompressedVector<int,rowVector>     SparseVector;

   DenseVector  d;
   SparseVector s;
   // ... Resizing and initialization

   // Creating an aligned dense subvector of size 8 starting from index 4
   blaze::DenseSubvector<DenseVector,aligned> dsv = subvector<aligned>( d, 4UL, 8UL );

   // Creating an unaligned subvector of size 7 starting from index 3
   blaze::SparseSubvector<SparseVector,unaligned> ssv = subvector<unaligned>( s, 3UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is larger
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions and the given \a index is subject to additional checks to
// guarantee proper alignment. However, especially in case of dense subvectors this may result
// in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of the subvector must be aligned. The following source code gives some examples
// for a double precision dynamic vector, assuming that AVX is available, which packs 4 \c double
// values into a SIMD vector:

   \code
   using blaze::columnVector;

   typedef blaze::DynamicVector<double,columnVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType,aligned>  SubvectorType;

   VectorType d( 17UL );
   // ... Resizing and initialization

   // OK: Starts at the beginning, i.e. the first element is aligned
   SubvectorType dsv1 = subvector<aligned>( d, 0UL, 13UL );

   // OK: Start index is a multiple of 4, i.e. the first element is aligned
   SubvectorType dsv2 = subvector<aligned>( d, 4UL, 7UL );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   SubvectorType dsv3 = subvector<aligned>( d, 8UL, 9UL );

   // Error: Start index is not a multiple of 4, i.e. the first element is not aligned
   SubvectorType dsv4 = subvector<aligned>( d, 5UL, 8UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline DisableIf_< Or< IsComputation<VT>, IsTransExpr<VT> >, SubvectorExprTrait_<const VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return SubvectorExprTrait_<const VT,AF>( ~vector, index, size );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector addition.
// \ingroup views
//
// \param vector The constant vector/vector addition.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the addition.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector addition.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecVecAddExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) +
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector subtraction.
// \ingroup views
//
// \param vector The constant vector/vector subtraction.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the subtraction.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector subtraction.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecVecSubExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) -
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector multiplication.
// \ingroup views
//
// \param vector The constant vector/vector multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector multiplication.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecVecMultExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) *
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector division.
// \ingroup views
//
// \param vector The constant vector/vector division.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector division.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecVecDivExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) /
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector cross product.
// \ingroup views
//
// \param dv The constant vector/vector cross product.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the cross product.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector cross product.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsCrossExpr<VT>, SubvectorExprTrait_<VT,unaligned> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return SubvectorExprTrait_<VT,unaligned>( ~vector, index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar multiplication.
// \ingroup views
//
// \param vector The constant vector/scalar multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar multiplication.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecScalarMultExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand(), index, size ) * (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar division.
// \ingroup views
//
// \param vector The constant vector/scalar division.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar division.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecScalarDivExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand(), index, size ) / (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector custom operation.
// \ingroup views
//
// \param vector The constant vector custom operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the custom operation.
//
// This function returns an expression representing the specified subvector of the given vector
// custom operation.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecForEachExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return forEach( subvector<AF>( (~vector).operand(), index, size ), (~vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector evaluation operation.
// \ingroup views
//
// \param vector The constant vector evaluation operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the evaluation operation.
//
// This function returns an expression representing the specified subvector of the given vector
// evaluation operation.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecEvalExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return eval( subvector<AF>( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector serialization operation.
// \ingroup views
//
// \param vector The constant vector serialization operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the serialization operation.
//
// This function returns an expression representing the specified subvector of the given vector
// serialization operation.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecSerialExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return serial( subvector<AF>( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector transpose operation.
// \ingroup views
//
// \param vector The constant vector transpose operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the transpose operation.
//
// This function returns an expression representing the specified subvector of the given vector
// transpose operation.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline EnableIf_< IsVecTransExpr<VT>, SubvectorExprTrait_<VT,AF> >
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return trans( subvector<AF>( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
