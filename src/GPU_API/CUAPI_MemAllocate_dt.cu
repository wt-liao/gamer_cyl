#include "CUAPI.h"

#ifdef GPU



extern real    *d_dt_Array_T;
extern real   (*d_Flu_Array_T)[NCOMP_FLUID][ CUBE(PS1) ];
extern double (*d_Corner_Array_T)[3];
#ifdef GRAVITY
extern real   (*d_Pot_Array_T)[ CUBE(GRA_NXT) ];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate_dt
// Description :  Allocate GPU and CPU memory for the dt solver
//
// Parameter   :  dt_NPG : Number of patch groups evaluated simultaneously by GPU for the dt solver
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate_dt( const int dt_NPG )
{

// size of the global memory arrays
   const int  dt_NP            = 8*dt_NPG;
   const long dt_MemSize_T     = sizeof(real  )*dt_NP;
   const long Flu_MemSize_T    = sizeof(real  )*dt_NP*NCOMP_FLUID*CUBE(PS1);
   const long Corner_MemSize_T = sizeof(double)*dt_NP*3;
#  ifdef GRAVITY
   const long Pot_MemSize_T    = sizeof(real  )*dt_NP*CUBE(GRA_NXT);
#  endif


// output the total memory requirement
   long TotalSize = dt_MemSize_T + Flu_MemSize_T + Corner_MemSize_T;
#  ifdef GRAVITY
   TotalSize += Pot_MemSize_T;
#  endif

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : total memory requirement in GPU dt solver = %ld MB\n", TotalSize/(1<<20) );


// allocate the device memory
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_dt_Array_T,               dt_MemSize_T            )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flu_Array_T,              Flu_MemSize_T           )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Corner_Array_T,           Corner_MemSize_T        )  );
#  ifdef GRAVITY
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_T,              Pot_MemSize_T           )  );
#  endif


// allocate the host memory by CUDA
   for (int t=0; t<2; t++)
   {
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_dt_Array_T    [t], dt_MemSize_T            )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_T   [t], Flu_MemSize_T           )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Corner_Array_T[t], Corner_MemSize_T        )  );
#     ifdef GRAVITY
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_T   [t], Pot_MemSize_T           )  );
#     endif
   } // for (int t=0; t<2; t++)

} // FUNCTION : CUAPI_MemAllocate_dt



#endif // #ifdef GPU
