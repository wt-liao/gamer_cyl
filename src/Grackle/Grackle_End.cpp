#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_End
// Description :  Free the resources used by Grackle
//
// Note        :  1. Invoked by End_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Grackle_End()
{

// nothing to do if Grackle is disabled
   if ( !GRACKLE_ACTIVATE )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   delete grackle_data;
   grackle_data = NULL;
   
   
#  ifdef GRACKLE_H2_SOBOLEV
   
   if (H2_Op_T_Table != NULL) {
      delete [] H2_Op_T_Table;
      H2_Op_T_Table = NULL;
   }
   
   if (H2_Op_Alpha_Table != NULL) {
      delete [] H2_Op_Alpha_Table;
      H2_Op_Alpha_Table = NULL;
   }
   
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Grackle_End



#endif // #ifdef SUPPORT_GRACKLE
