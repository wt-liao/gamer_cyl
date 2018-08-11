#include "GAMER.h"

#if ( !defined GPU  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_PoissonGravity
// Description :  Free memory previously allocated by the function "Init_MemAllocate_PoissonGravity"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_PoissonGravity()
{

   for (int t=0; t<2; t++)
   {
      if ( h_Rho_Array_P    [t] != NULL )    delete [] h_Rho_Array_P    [t];
      if ( h_Pot_Array_P_In [t] != NULL )    delete [] h_Pot_Array_P_In [t];
      if ( h_Pot_Array_P_Out[t] != NULL )    delete [] h_Pot_Array_P_Out[t];
#     ifdef UNSPLIT_GRAVITY
      if ( h_Pot_Array_USG_G[t] != NULL )    delete [] h_Pot_Array_USG_G[t];
      if ( h_Flu_Array_USG_G[t] != NULL )    delete [] h_Flu_Array_USG_G[t];
#     endif
      if ( h_Flu_Array_G    [t] != NULL )    delete [] h_Flu_Array_G    [t];
      if ( h_Corner_Array_G [t] != NULL )    delete [] h_Corner_Array_G [t];
#     ifdef DUAL_ENERGY
      if ( h_DE_Array_G     [t] != NULL )    delete [] h_DE_Array_G     [t];
#     endif

      h_Rho_Array_P    [t] = NULL;
      h_Pot_Array_P_In [t] = NULL;
      h_Pot_Array_P_Out[t] = NULL;
#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_G[t] = NULL;
      h_Flu_Array_USG_G[t] = NULL;
#     endif
      h_Flu_Array_G    [t] = NULL;
      h_Corner_Array_G [t] = NULL;
#     ifdef DUAL_ENERGY
      h_DE_Array_G     [t] = NULL;
#     endif
   }


#  if ( COORDINATE == CARTESIAN )
   if ( GreenFuncK != NULL )  delete [] GreenFuncK;
   
#  elif ( COORDINATE == CYLINDRICAL )
   // 1.0 free memory
   Aux_DeallocateArray2D(KernelFuncK);
   // 1.1
   Aux_DeallocateArray2D(RhoK);
   Aux_DeallocateArray2D(PhiK);
   // 1.2
   delete [] SendBuf_Rho ; 
   delete [] SendBuf_IDPlanXp ;
   delete [] SendBuf_IDPlanYZ ;
   delete [] RecvBuf_Rho ; 
   delete [] RecvBuf_IDPlanXp ;
   delete [] RecvBuf_IDPlanYZ ;
   // 1.3
   delete [] PhiK_All_re ; 
   delete [] PhiK_All_im ;
   if (RANK_IP == 0) {
      delete [] PhiK_local_re ;
      delete [] PhiK_local_im ;
   }
   // 1.4
   if (RANK_IP == 0) {
      delete [] SendBuf_Phi ; 
      delete [] SendBuf_PID ;
      delete [] SendBuf_I ;
   }
   delete [] RecvBuf_Phi ;
   delete [] RecvBuf_PID ;
   delete [] RecvBuf_I ;
   
   
   // 2.0 free new MPI_comm
   MPI_Comm_free(&rank_i_comm );
   MPI_Comm_free(&rank_ip_comm);
   
#  endif // COORDINATE 

} // FUNCTION : End_MemFree_PoissonGravity

#endif // #if ( !defined GPU  &&  defined GRAVITY )
