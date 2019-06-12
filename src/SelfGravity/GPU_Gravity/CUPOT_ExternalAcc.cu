#ifdef __CUDACC__
#include "Macro.h"
#else
#include "GAMER.h"
#endif
#include "CUPOT.h"

#ifdef GRAVITY


// soften length implementation
//#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT

//extern real Rayleigh_Disk_Slope;


//-----------------------------------------------------------------------------------------
// Function    :  CUPOT_ExternalAcc / CPU_ExternlAcc
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function will be invoked by both CPU and GPU
//                2. "__forceinline__" is required since this device function will be invoked
//                   by more than one kernels (e.g., CUPOT_HydroGravitySolver, CUFLU_ComputeFlux)
//                3. The auxiliary array "UserArray" is set by "Init_ExternalAcc_Ptr", which
//                   points to "Init_ExternalAcc()" by default but may be overwritten by various
//                   test problem initializers
//                4. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                     UserArray[4] = soften_length (<=0.0 --> disable)
//                   --> but one can easily modify this file to change the default behavior
//                5. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                X/Y/Z     : Target spatial coordinates in the adopted coordinate system
//                Time      : Current physical time
//                UserArray : User-provided auxiliary array (set by "Init_ExternalAcc_Ptr")
//
// Return      :  Acc
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__forceinline__ __device__
void CUPOT_ExternalAcc( real Acc[], const double X, const double Y, const double Z, const double Time, const double UserArray[] )
#else
void   CPU_ExternalAcc( real Acc[], const double X, const double Y, const double Z, const double Time, const double UserArray[] )
#endif
{
   //### customize for Rayleigh disk test problem
   if (TESTPROB_ID == 21) {
      const real Rayleigh_Disk_Slope = UserArray[0] ;
      
      Acc[0] = -POW( X, -Rayleigh_Disk_Slope ) ;
      //Acc[0] = -POW( X, -2.0 ) ;
      Acc[1] = (real) 0.0 ;
      Acc[2] = (real) 0.0 ;
      
      return ;
   }
   
   //### customize for Kuzmin disk test problem
   if (TESTPROB_ID == 23) {
      
      Acc[0] = (real) 0.0 ;
      Acc[1] = (real) 0.0 ;
      Acc[2] = (real) 0.0 ;
      
      return ;
   }
   
   
   const double Cen[3] = { UserArray[0], UserArray[1], UserArray[2] };
   const real GM       = (real)UserArray[3];
   const real eps      = (real)UserArray[4];

   const real dx         = (real)(X - Cen[0]);
   const real dy         = (real)(Y - Cen[1]);
   const real dz         = (real)(Z - Cen[2]);
   
#  if ( COORDINATE == CARTESIAN )
   const real r          = SQRT( dx*dx + dy*dy + dz*dz );
#  elif ( COORDINATE == CYLINDRICAL )
   const double cos_dtheta = cos( Y-Cen[1] ) ;
   const double sin_dtheta = sin( Y-Cen[1] );
   const real   r          = SQRT( X*X + Cen[0]*Cen[0] - 2*X*Cen[0]*cos_dtheta + dz*dz );
#  endif

// Plummer
#  if   ( defined SOFTEN_PLUMMER )
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif ( defined SOFTEN_RUFFERT )
   const real tmp = EXP( -SQR(r)/SQR(eps) );
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
   const real _r3 = (real)1.0/CUBE(r);
#  endif

#  if ( COORDINATE == CARTESIAN )
   Acc[0] = -GM*_r3*dx;
   Acc[1] = -GM*_r3*dy;
   Acc[2] = -GM*_r3*dz;
#  elif ( COORDINATE == CYLINDRICAL )
   Acc[0] = -GM*_r3 * ( X-Cen[0]*cos_dtheta );
   Acc[1] = -GM*_r3 * ( Cen[0]*sin_dtheta );
   Acc[2] = -GM*_r3 * dz;
#  endif
   

   
   

} // FUNCTION : CUPOT_ExternalAcc / CPU_ExternalAcc



#endif // #ifdef GRAVITY
