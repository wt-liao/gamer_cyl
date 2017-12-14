#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Region
// Description :  Check if the element (i,j,k) of the input patch is within the regions allowed to be refined
//
// Note        :  To use this functionality, please turn on the option "OPT__FLAG_REGION" and then specify the
//                target regions in this file
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[0][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//
// Return      :  "true/false"  if the input cell "is/is not" within the region allowed for refinement
//-------------------------------------------------------------------------------------------------------
bool Flag_Region( const int i, const int j, const int k, const int lv, const int PID )
{

   const double Pos[3] = { Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i ),
                           Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 1, j ),
                           Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 2, k ) };

   bool Within = true;


// put the target region below
// ##########################################################################################################
/*
// Example : sphere
   const double dR[3]     = { Pos[0] - amr->BoxCenter[0],
                              Pos[1] - amr->BoxCenter[1],
                              Pos[2] - amr->BoxCenter[2] };
   const double R         = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );
   const double MaxR      = 1.0;

   Within = R <= MaxR;
*/
// ##########################################################################################################


   return Within;

} // FUNCTION : Flag_Region

