#include "GAMER.h"

#ifndef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_dt
// Description :  Allocate memory for the dt solver
//
// Note        :  Work when using CPUs only
//
// Parameter   :  dt_NPatchGroup : Number of patch groups calculated at a time
//
// Return      :  h_dt_Array_T, h_Flu_Array_T, h_Corner_Array_T, h_Pot_Array_T
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_dt( const int dt_NPatchGroup )
{

   const int dt_NP = 8*dt_NPatchGroup;

   for (int t=0; t<2; t++)
   {
      h_dt_Array_T    [t] = new real   [dt_NP];
      h_Flu_Array_T   [t] = new real   [dt_NP][NCOMP_FLUID][ CUBE(PS1) ];
      h_Corner_Array_T[t] = new double [dt_NP][3];
#     ifdef GRAVITY
      h_Pot_Array_T   [t] = new real   [dt_NP][ CUBE(GRA_NXT) ];
#     endif
   }

} // FUNCTION : Init_MemAllocate_dt



#endif // #ifndef GPU
