#include "CUAPI.h"

#ifdef GPU



extern real    *d_dt_Array_T;
extern real   (*d_Flu_Array_T)[NCOMP_FLUID][ CUBE(PS1) ];
extern double (*d_Corner_Array_T)[3];
#ifdef GRAVITY
extern real   (*d_Pot_Array_T)[ CUBE(GRA_NXT) ];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree_dt
// Description :  Free the GPU and CPU memory previously allocated by CUAPI_MemAllocate_dt()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemFree_dt()
{

// free the device memory
   if ( d_dt_Array_T     != NULL )  CUDA_CHECK_ERROR(  cudaFree( d_dt_Array_T     )  );
   if ( d_Flu_Array_T    != NULL )  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_T    )  );
   if ( d_Corner_Array_T != NULL )  CUDA_CHECK_ERROR(  cudaFree( d_Corner_Array_T )  );
#  ifdef GRAVITY
   if ( d_Pot_Array_T    != NULL )  CUDA_CHECK_ERROR(  cudaFree( d_Pot_Array_T    )  );
#  endif

   d_dt_Array_T     = NULL;
   d_Flu_Array_T    = NULL;
   d_Corner_Array_T = NULL;
#  ifdef GRAVITY
   d_Pot_Array_T    = NULL;
#  endif


// free the host memory allocated by CUDA
   for (int t=0; t<2; t++)
   {
      if ( h_dt_Array_T    [t] != NULL )  CUDA_CHECK_ERROR(  cudaFreeHost( h_dt_Array_T    [t] )  );
      if ( h_Flu_Array_T   [t] != NULL )  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_T   [t] )  );
      if ( h_Corner_Array_T[t] != NULL )  CUDA_CHECK_ERROR(  cudaFreeHost( h_Corner_Array_T[t] )  );
#     ifdef GRAVITY
      if ( h_Pot_Array_T   [t] != NULL )  CUDA_CHECK_ERROR(  cudaFreeHost( h_Pot_Array_T   [t] )  );
#     endif

      h_dt_Array_T    [t] = NULL;
      h_Flu_Array_T   [t] = NULL;
      h_Corner_Array_T[t] = NULL;
#     ifdef GRAVITY
      h_Pot_Array_T   [t] = NULL;
#     endif
   } // for (int t=0; t<2; t++)

} // FUNCTION : CUAPI_MemFree_dt


#endif // #ifdef GPU
