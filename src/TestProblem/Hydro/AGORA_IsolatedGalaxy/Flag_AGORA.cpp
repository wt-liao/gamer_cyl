#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_AGORA
// Description :  Flag cells for refinement for the AGORA isolated galaxy simulation
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by "Init_TestProb_Hydro_AGORA_IsolatedGalaxy()"
//                   to replace "Flag_User()"
//                2. Please turn on the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the targeted element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the targeted patch
//                PID       : ID of the targeted patch
//                Threshold : Useless here
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_AGORA( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{

   const double Pos[3] = { Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i ),
                           Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 1, j ),
                           Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 2, k ) };

// flag cells within the target region [BoxEdgeL+Threshold ... BoxEdgeR-Threshold]
   const double EdgeL[3] = { amr->BoxEdgeL[0] + Threshold,
                             amr->BoxEdgeL[1] + Threshold,
                             amr->BoxEdgeL[2] + Threshold };
   const double EdgeR[3] = { amr->BoxEdgeR[0] - Threshold,
                             amr->BoxEdgeR[1] - Threshold,
                             amr->BoxEdgeR[2] - Threshold };

   bool Flag;

   if (  Pos[0] >= EdgeL[0]  &&  Pos[0] < EdgeR[0]  &&
         Pos[1] >= EdgeL[1]  &&  Pos[1] < EdgeR[1]  &&
         Pos[2] >= EdgeL[2]  &&  Pos[2] < EdgeR[2]     )
      Flag = true;

   else
      Flag = false;


   return Flag;

} // FUNCTION : Flag_AGORA

