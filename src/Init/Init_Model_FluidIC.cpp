#ifdef MODEL_IC_FLUID

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FluidField
// Description :  
//
// Note        :  1. only valid for no0amr scheme; i.e., strictly, lv=0 
//                2. need ghostzone density, potential, pressure 
//                   -> apply Prepare_PatchData
//                
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------

void Init_FluidField(){
   
   // parameters 
   const int lv            = 0 ;
   const int IC_Ghost_Size = 1 ; 
   const int IC_NXT        = 2*(PATCH_SIZE+IC_Ghost_Size) ;
   const int NPG_Max       = FLU_GPU_NPGROUP;
   
   int *PID0_List          = NULL;  // list recording the patch indicies with LocalID==0 to be udpated
   int  NPG;                        // number of patch groups to be updated at a time
   int  NTotal;                     // total number of patch groups to be updated
   int  Disp;                       // index displacement in PID0_List
   
   NTotal       = amr->NPatchComma[lv][1] / 8;
   PID0_List    = new int [NTotal];
   for (int t=0; t<NTotal; t++)  PID0_List[t] = 8*t;
   
   
   //### use new 
   // array to store prepared data; 
   // FLU_NIN->NCOMP_TOTAL, FLU_NXT->2*(PATCH_SIZE+FLU_GHOST_SIZE)
   // allocate memory for store fluid+potential field
   real (*IC_Flu_Array)[FLU_NIN][CUBE(IC_NXT)];
   real (*IC_Pot_Array)[CUBE(IC_NXT)];
   
   IC_Flu_Array = new real [FLU_GPU_NPGROUP][FLU_NIN][CUBE(IC_NXT)];
   IC_Pot_array = new real [FLU_GPU_NPGROUP][CUBE(IC_NXT)];
   
   
   // 2.0
   for (int Disp=0; Disp<NTotal; Disp++) {
      NPG = 1; //###
      IC_Find_Pressure(lv, NPG, PID0_List, IC_Flu_Array, IC_Pot_Array, IC_Ghost_Size);
   }
   
   // 3.0
   MPI_Barrier( MPI_COMM_WORLD );   
   
   // 4.0.
   for (int Disp=0; Disp<NTotal; Disp+=NPG_Max) {
      IC_Find_Velocity();
   }
   
   // 5.0 for dual_energy
   
   
   
   // 6.0 delete allocated IC array
   if (IC_Flu_Array != NULL) {
      delete [] IC_Flu_Array;
   }
   
   if (IC_Pot_Array != NULL) {
      delete [] IC_Pot_Array;
   }

}


//-------------------------------------------------------------------------------------------------------
// Function    :  IC_Find_Pressure
// Description :  
//
// Note        :  1. work for Init_FluidField(); 
//                
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------

void IC_Find_Pressure(const int lv, const int NPG, const int *PID0_List){
   
   const double PrepTime      = Time[lv];
   const bool   IntPhase_No   = false;
   
   //### copy from Flu_Prepare
   
   // 1.0 prepare for patch data
   // 1.1 prepare fluid 
   Prepare_PatchData( lv, PrepTime, IC_Flu_Array, IC_Ghost_Size, NPG, PID0_List, _TOTAL,
                      OPT__FLU_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, BC_POT_NONE, MinDens, MinPres, DE_Consistency );
                      
   // 1.2 prepare potential
   //### why OPT__BC_FLU? check DE_Consistency_No? 
   Prepare_PatchData( lv, PrepTime, IC_Pot_Array, IC_Ghost_Size, NPG, PID0_List, _POTE,
                      OPT__GRA_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, DE_Consistency_No );
   
   // 2.0 Loop over all paches; compute for pressure field
#  pragma omp for schedule( runtime )
   for (int P=0; P<NPG; P++)
   {
      for (int k=IC_Ghost_Size; k<PS2+IC_Ghost_Size; k++)
      for (int j=IC_Ghost_Size; j<PS2+IC_Ghost_Size; j++)
      for (int i=IC_Ghost_Size; i<PS2+IC_Ghost_Size; i++)
      {
         
      } // for (k, j, i)
      
   } // OpenMP parallel region

}


#endif   // MODEL_IC_FLUID