#include "GAMER.h"

#ifndef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_dt
// Description :  Free memory previously allocated by the function "Init_MemAllocate_dt"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_dt()
{

   for (int t=0; t<2; t++)
   {
      if ( h_dt_Array_T    [t] != NULL )  delete [] h_dt_Array_T    [t];
      if ( h_Flu_Array_T   [t] != NULL )  delete [] h_Flu_Array_T   [t];
      if ( h_Corner_Array_T[t] != NULL )  delete [] h_Corner_Array_T[t];
#     ifdef GRAVITY
      if ( h_Pot_Array_T   [t] != NULL )  delete [] h_Pot_Array_T   [t];
#     endif

      h_dt_Array_T    [t] = NULL;
      h_Flu_Array_T   [t] = NULL;
      h_Corner_Array_T[t] = NULL;
#     ifdef GRAVITY
      h_Pot_Array_T   [t] = NULL;
#     endif
   } // for (int t=0; t<2; t++)

} // FUNCTION : End_MemFree_dt



#endif // #ifndef GPU
