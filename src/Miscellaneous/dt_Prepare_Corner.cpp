#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  dt_Prepare_Corner
// Description :  Prepare the input array "h_Corner_Array_T" for the dt solver
//
// Note        :  1. The corner coordinates are mainly used for non-Cartesian coordinates and external force
//                2. Corner coordinates are defined as the central coordinates of the first cell located
//                   at the bottom left corner
//
// Parameter   :  lv               : Target refinement level
//                h_Corner_Array_T : Host array to store the prepared data
//                NPG              : Number of patch groups prepared at a time
//                PID0_List        : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void dt_Prepare_Corner( const int lv, double h_Corner_Array_T[][3], const int NPG, const int *PID0_List )
{

   const double dh_half[3] = { 0.5*amr->dh[lv][0], 0.5*amr->dh[lv][1], 0.5*amr->dh[lv][2] };

   int N, PID, PID0;


#  pragma omp parallel for private( N, PID, PID0 ) schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;
         N   = 8*TID + LocalID;

         for (int d=0; d<3; d++)    h_Corner_Array_T[N][d] = amr->patch[0][lv][PID]->EdgeL[d] + dh_half[d];
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : dt_Prepare_Corner
