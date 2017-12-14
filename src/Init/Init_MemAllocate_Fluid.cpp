#include "GAMER.h"

#ifndef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Fluid
// Description :  Allocate memory for the fluid solver
//
// Note        :  Work when using CPUs only
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Fluid( const int Flu_NPatchGroup )
{

// determine whether or not to allocate the corner array
   bool AllocateCorner = false;

#  ifdef UNSPLIT_GRAVITY
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
      AllocateCorner = true;
#  endif

#  if ( COORDINATE != CARTESIAN )
      AllocateCorner = true;
#  endif


   for (int t=0; t<2; t++)
   {
      h_Flu_Array_F_In [t] = new real [Flu_NPatchGroup][FLU_NIN ][   FLU_NXT   *FLU_NXT   *FLU_NXT    ];
      h_Flu_Array_F_Out[t] = new real [Flu_NPatchGroup][FLU_NOUT][ 8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ];

      if ( amr->WithFlux )
      h_Flux_Array     [t] = new real [Flu_NPatchGroup][9][NFLUX_TOTAL][ 4*PATCH_SIZE*PATCH_SIZE ];

#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_F[t] = new real [Flu_NPatchGroup][USG_NXT_F][USG_NXT_F][USG_NXT_F];
#     endif

      if ( AllocateCorner )
      h_Corner_Array_F [t] = new double [Flu_NPatchGroup][3];

#     ifdef DUAL_ENERGY
      h_DE_Array_F_Out [t] = new char [Flu_NPatchGroup][ 8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ];
#     endif
   }

} // FUNCTION : Init_MemAllocate_Fluid



#endif // #ifndef GPU
