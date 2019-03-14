#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void BC_User( real fluid[], const double X, const double Y, const double Z, const double Time,
                     const int lv, double AuxArray[] );

// this function pointer may be overwritten by various test problem initializers
void (*BC_User_Ptr)( real fluid[], const double X, const double Y, const double Z, const double Time,
                     const int lv, double AuxArray[] ) = BC_User;


#if ( COORDINATE == CYLINDRICAL )
#include "TestProb.h"
                     
// boundary condition using gradient field
static void BC_User_xm( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, 
                        const int ArraySizeY, const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar );
static void BC_User_xp( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, 
                        const int ArraySizeY, const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar );
static void BC_User_ym( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, 
                        const int ArraySizeY, const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar );
static void BC_User_yp( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, 
                        const int ArraySizeY, const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar );
static void BC_User_zm( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, 
                        const int ArraySizeY, const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar );
static void BC_User_zp( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, 
                        const int ArraySizeY, const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                        const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar );

static const double rho_ratio_limit = 0.6;

#endif


//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Always return NCOMP_TOTAL variables
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  fluid    : Fluid field to be set
//                X/Y/Z    : Physical coordinates in the adopted coordinate system
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC_User( real fluid[], const double X, const double Y, const double Z, const double Time,
              const int lv, double AuxArray[] )
{

// please put your B.C. here
// ##########################################################################################################
// Example 1 : set to time-independent values for HYDRO
   /*
   const double *C   = amr->BoxCenter;
   const real Height = 100.0;
   const real Width  =  64.0;
   const real Gamma2 = real( 1.0/GAMMA/(GAMMA-1.0) );
   const real Cs     = 1.0;
   const real Rho0   = 1.0;

   fluid[DENS] = Rho0 + Height*EXP(  -( SQR(X-C[0]) + SQR(Y-C[1]) + SQR(Z-C[2]) ) / SQR(Width)  );
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = Cs*Cs*fluid[DENS]*Gamma2 + (real)0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

// remember to set passive scalars as well
   fluid[EINT] = XXX;
   */


// Example 2 : set to time-dependent values for HYDRO

// ##########################################################################################################

} // FUNCTION : BC_User



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_BoundaryCondition_User
// Description :  Fill up the ghost-zone values by the user-specified boundary condition
//
// Note        :  1. Work for the functions "Prepare_PatchData, InterpolateGhostZone, Refine, LB_Refine_AllocateNewPatch"
//                2. The function pointer "BC_User_Ptr" points to "BC_User()" by default but may be overwritten
//                   by various test problem initializers
//
// Parameter   :  Array          : Array to store the prepared data of one patch group (including the ghost-zone data)
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables are NOT included)
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                Time           : Current physical time
//                dh             : Grid size
//                Corner         : Physcial coordinates at the center of the cell (0,0,0) --> Array[0]
//                TVar           : Target variables to be prepared --> only used for preparing the derived variables
//                lv             : Refinement level
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void Flu_BoundaryCondition_User( real *Array, real *PotArray, const int BC_Face, const int NVar_Flu, const int GhostSize, 
                                 const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ, 
                                 const int Idx_Start[], const int Idx_End[], const int TFluVarIdxList[], const double Time, 
                                 const double dh[], const double *Corner, const int TVar, const int lv )
{

// customized USER_BC
   switch ( BC_Face )
   {
      case 0:  BC_User_xm( Array, PotArray, NVar_Flu, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, 
                           TFluVarIdxList, dh, Corner, TVar );  break;
      case 1:  BC_User_xp( Array, PotArray, NVar_Flu, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, 
                           TFluVarIdxList, dh, Corner, TVar );  break;
      case 2:  BC_User_ym( Array, PotArray, NVar_Flu, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, 
                           TFluVarIdxList, dh, Corner, TVar );  break;
      case 3:  BC_User_yp( Array, PotArray, NVar_Flu, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, 
                           TFluVarIdxList, dh, Corner, TVar );  break;
      case 4:  BC_User_zm( Array, PotArray, NVar_Flu, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, 
                           TFluVarIdxList, dh, Corner, TVar );  break;
      case 5:  BC_User_zp( Array, PotArray, NVar_Flu, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, Idx_Start, Idx_End, 
                           TFluVarIdxList, dh, Corner, TVar );  break;
      default: Aux_Error( ERROR_INFO, "incorrect boundary face (%d) !!\n", BC_Face );
   }
   
// ### previous version of USER_BC
   
   /*
   
// check
   if ( BC_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "BC_User_Ptr == NULL !!\n" );  

// starting coordinates in the adopted coordinate system
   const double X0 = Corner[0] + (double)Idx_Start[0]*dh[0];
   const double Y0 = Corner[1] + (double)Idx_Start[1]*dh[1];
   const double Z0 = Corner[2] + (double)Idx_Start[2]*dh[2];

#  if   ( MODEL == HYDRO )
   const bool CheckMinPres_Yes = true;
   const real Gamma_m1         = GAMMA - (real)1.0;
   const bool PrepVx           = ( TVar & _VELX ) ? true : false;
   const bool PrepVy           = ( TVar & _VELY ) ? true : false;
   const bool PrepVz           = ( TVar & _VELZ ) ? true : false;
   const bool PrepPres         = ( TVar & _PRES ) ? true : false;
   const bool PrepTemp         = ( TVar & _TEMP ) ? true : false;

#  elif ( MODEL == MHD   )
#  warning : WAIT MHD !!

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif


// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   int    i, j, k, v2;
   real   BVal[NCOMP_TOTAL];
   double X, Y, Z;

   for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2]; k++, Z+=dh[2])
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
   {
      BC_User_Ptr( BVal, X, Y, Z, Time, lv, NULL );

      for (int v=0; v<NVar_Flu; v++)   Array3D[v][k][j][i] = BVal[ TFluVarIdxList[v] ];


//    derived variables
      v2 = NVar_Flu;

#     if   ( MODEL == HYDRO )
      if ( PrepVx   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMX] / BVal[DENS];
      if ( PrepVy   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMY] / BVal[DENS];
      if ( PrepVz   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMZ] / BVal[DENS];
      if ( PrepPres )   Array3D[ v2 ++ ][k][j][i] = CPU_GetPressure( BVal[DENS], BVal[MOMX], BVal[MOMY], BVal[MOMZ], BVal[ENGY],
                                                                     Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      if ( PrepTemp )   Array3D[ v2 ++ ][k][j][i] = CPU_GetTemperature( BVal[DENS], BVal[MOMX], BVal[MOMY], BVal[MOMZ], BVal[ENGY],
                                                                        Gamma_m1, CheckMinPres_Yes, MIN_PRES );

#     elif ( MODEL == MHD   )
#     warning : WAIT MHD !!

#     elif ( MODEL == ELBDM )
//    no derived variables yet

#     else
#     error : unsupported MODEL !!
#     endif
   } // k,j,i
   */

} // FUNCTION : Flu_BoundaryCondition_User





#if (COORDINATE == CYLINDRICAL)
//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_xm
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. update Array 
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  Array    : Fluid field to be set
//
//-------------------------------------------------------------------------------------------------------
void BC_User_xm( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar ) 
{
#  ifdef UserPotBC
// 1D array -> 3D array
   const int FACE = 0;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   if (PotArray != NULL) {
      Poi_BoundaryCondition_Extrapolation( PotArray, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
                                           Idx_Start, Idx_End, NULL, NULL );
   }
   real (*PotArray3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )PotArray;
   
   const double X0    = Corner[0] + (double)Idx_End[0]*dh[0];
   const double Y0    = Corner[1] + (double)Idx_Start[1]*dh[1];
   const double Z0    = Corner[2] + (double)Idx_Start[2]*dh[2];
   const double GM    = ExtAcc_AuxArray[3] ;
   const double X_ref = X0 + dh[0];
   const double Y_ref = Y0;
   const double Z_ref = Z0;
   const int    i_ref = Idx_End[0]+1 ;
   
   const bool CheckMinPres_Yes = true;
   const real Gamma_m1         = GAMMA - (real)1.0;
   
   double engy_bc, pres_ref, _rho_ref;
   
   double X, Y, Z;
   int    i, j, k;

   for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2];   k++, Z+=dh[2])
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1];   j++, Y+=dh[1])
   {
      pres_ref  = CPU_GetPressure( Array3D[DENS][k][j][i_ref], Array3D[MOMX][k][j][i_ref], 
                                   Array3D[MOMY][k][j][i_ref], Array3D[MOMZ][k][j][i_ref], 
                                   Array3D[ENGY][k][j][i_ref], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      _rho_ref  = 1.0/Array3D[DENS][k][j][i_ref];
                                   
      for (i=Idx_End[0], X=X0; i>=Idx_Start[0]; i--, X-=dh[0])
      {
         // outflow 
         Array3D[DENS][k][j][i] = Array3D[DENS][k][j][i_ref] ;
         Array3D[MOMX][k][j][i] = Array3D[MOMX][k][j][i_ref] ;
         Array3D[MOMY][k][j][i] = (Array3D[MOMY][k][j][i_ref] < 0)? Array3D[MOMY][k][j][i_ref] : 0 ;
         Array3D[MOMZ][k][j][i] = Array3D[MOMZ][k][j][i_ref] ;
      
         Array3D[ENGY][k][j][i] = pres_ref/Gamma_m1 
                                + 0.5*( SQR(Array3D[MOMX][k][j][i]) + SQR(Array3D[MOMY][k][j][i]) + 
                                        SQR(Array3D[MOMZ][k][j][i]) ) / Array3D[DENS][k][j][i] ;
         
#        ifdef SUPPORT_GRACKLE
         if (GRACKLE_PRIMORDIAL != GRACKLE_PRI_CHE_NSPE9)
            Aux_Message(stderr, "User defined BC curretly only support Grackle with 9 spices! \n");
         
         for (int n=Idx_e; n<=Idx_H2II; n++) {
            Array3D[n][k][j][i] = Array3D[n][k][j][i_ref] * ( Array3D[DENS][k][j][i]*_rho_ref );
         }
#        endif //#SUPPORT_GRACKLE
         
#        if ( DUAL_ENERGY == DE_ENPY )
         Array3D[Idx_Enpy][k][j][i] = 
            CPU_Fluid2Entropy(Array3D[DENS][k][j][i], Array3D[MOMX][k][j][i], Array3D[MOMY][k][j][i],
                              Array3D[MOMZ][k][j][i], Array3D[ENGY][k][j][i], Gamma_m1);
#        endif // DUAL_ENERGY == DE_ENPY
         
      } // i
   } // k, j
   
#  endif
}
                        
                        
//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_xp
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. update Array 
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  Array    : Fluid field to be set
//
//-------------------------------------------------------------------------------------------------------                        
void BC_User_xp( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar )
{
#  ifdef UserPotBC
// 1D array -> 3D array
   const int FACE = 1;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   if (PotArray != NULL) {
      Poi_BoundaryCondition_Extrapolation( PotArray, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
                                           Idx_Start, Idx_End, NULL, NULL );
   }
   real (*PotArray3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )PotArray;
   
   const double X0    = Corner[0] + (double)Idx_Start[0]*dh[0];
   const double Y0    = Corner[1] + (double)Idx_Start[1]*dh[1];
   const double Z0    = Corner[2] + (double)Idx_Start[2]*dh[2];
   const double GM    = ExtAcc_AuxArray[3] ;
   const double X_ref = X0 - dh[0];
   const double Y_ref = Y0;
   const double Z_ref = Z0;
   
   const int    i_ref = Idx_Start[0]-1 ;
   
   double sph_rad, star_g, dens, pot_grad, pres_grad, vtheta_square, vtheta;
   
   double X, Y, Z;
   int    i, j, k;
   
   
   for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2]; k++, Z+=dh[2])
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
   {
      // outflow 
      Array3D[DENS][k][j][i] = Array3D[DENS][k][j][i_ref] ;
      Array3D[MOMX][k][j][i] = Array3D[MOMX][k][j][i_ref] ;
      Array3D[MOMY][k][j][i] = Array3D[MOMY][k][j][i_ref] ;
      Array3D[MOMZ][k][j][i] = Array3D[MOMZ][k][j][i_ref] ;
      Array3D[ENGY][k][j][i] = Array3D[ENGY][k][j][i_ref] ;

   }
#  endif // UserPotBC
}
                        
//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_ym
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. update Array 
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  Array    : Fluid field to be set
//
//-------------------------------------------------------------------------------------------------------              
void BC_User_ym( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar )
{
#  ifdef UserPotBC
// 1D array -> 3D array
   const int FACE = 2;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   const double X0 = Corner[0] + (double)Idx_Start[0]*dh[0];
   const double Y0 = Corner[1] + (double)Idx_End[1]*dh[1]  ;
   const double Z0 = Corner[2] + (double)Idx_Start[2]*dh[2];
   
   const int    j_ref = Idx_End[1]+1 ;
   
   double X, Y, Z;
   int    i, j, k;
   
   for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2];   k++, Z+=dh[2])
   for (j=Idx_End[1],   Y=Y0; j>=Idx_Start[1]; j--, Y-=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0];   i++, X+=dh[0])
   {
      // outflow 
      Array3D[DENS][k][j][i] = Array3D[DENS][k][j_ref][i] ;
      Array3D[MOMX][k][j][i] = Array3D[MOMX][k][j_ref][i] ;
      Array3D[MOMY][k][j][i] = Array3D[MOMY][k][j_ref][i] ;
      Array3D[MOMZ][k][j][i] = Array3D[MOMZ][k][j_ref][i] ;
      Array3D[ENGY][k][j][i] = Array3D[ENGY][k][j_ref][i] ;
   }            
#  endif
}

//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_yp
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. update Array 
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  Array    : Fluid field to be set
//
//-------------------------------------------------------------------------------------------------------        
void BC_User_yp( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar )
{
#  ifdef UserPotBC
// 1D array -> 3D array
   const int FACE = 3;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   const double X0 = Corner[0] + (double)Idx_Start[0]*dh[0];
   const double Y0 = Corner[1] + (double)Idx_Start[1]*dh[1];
   const double Z0 = Corner[2] + (double)Idx_Start[2]*dh[2];
   
   const int    j_ref = Idx_Start[1]-1 ;
   
   double X, Y, Z;
   int    i, j, k;
   
   for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2]; k++, Z+=dh[2])
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
   {
      // outflow 
      Array3D[DENS][k][j][i] = Array3D[DENS][k][j_ref][i] ;
      Array3D[MOMX][k][j][i] = Array3D[MOMX][k][j_ref][i] ;
      Array3D[MOMY][k][j][i] = Array3D[MOMY][k][j_ref][i] ;
      Array3D[MOMZ][k][j][i] = Array3D[MOMZ][k][j_ref][i] ;
      Array3D[ENGY][k][j][i] = Array3D[ENGY][k][j_ref][i] ;  
   }
#  endif
}

//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_zm
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. update Array 
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  Array    : Fluid field to be set
//
//-------------------------------------------------------------------------------------------------------              
void BC_User_zm( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar )
{
#  ifdef UserPotBC
// 1D array -> 3D array
   const int FACE = 4;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   const double X0    = Corner[0] + (double)Idx_Start[0]*dh[0];
   const double Y0    = Corner[1] + (double)Idx_Start[1]*dh[1];
   const double Z0    = Corner[2] + (double)Idx_End[2]*dh[2];
   const double GM    = ExtAcc_AuxArray[3] ;
   const double X_ref = X0;
   const double Y_ref = Y0;
   const double Z_ref = Z0 + dh[2];
   const int    k_ref = Idx_End[2]+1 ;
   
   const bool CheckMinPres_Yes = true;
   const real Gamma_m1         = GAMMA - (real)1.0;
   
   double pres_bc, engy_bc, sph_rad, star_g, dens_p1, pres_p2, pot_grad;
   double pres_ref, rho_guess, rho_ref, _rho_ref, RT_ref, Vx_ref, Vy_ref, Vz_ref;
   double X, Y, Z, Z_p1;
   int    i, j, k;
   
   
   if (PotArray != NULL) {
      Poi_BoundaryCondition_Extrapolation( PotArray, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
                                           Idx_Start, Idx_End, NULL, NULL );
   }
   
   real (*PotArray3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )PotArray;
   
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
   {  
      rho_ref   = Array3D[DENS][k_ref][j][i] ;
      _rho_ref  = 1 / rho_ref;
      Vx_ref    = Array3D[MOMX][k_ref][j][i] * _rho_ref;
      Vy_ref    = Array3D[MOMY][k_ref][j][i] * _rho_ref;
      Vz_ref    = Array3D[MOMZ][k_ref][j][i] * _rho_ref;
      
      pres_ref  = CPU_GetPressure( Array3D[DENS][k_ref][j][i], Array3D[MOMX][k_ref][j][i], 
                                   Array3D[MOMY][k_ref][j][i], Array3D[MOMZ][k_ref][j][i], 
                                   Array3D[ENGY][k_ref][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      RT_ref    = pres_ref * _rho_ref;
      
      // fill in bc
      for (k=Idx_End[2], Z=Z0; k>=Idx_Start[2]; k--, Z-=dh[2]) {
         
         // calculate for BC value
         Z_p1     = Z+dh[2];
         sph_rad  = SQRT( Z_p1*Z_p1 + X*X );
         star_g   = - GM * Z_p1 / CUBE(sph_rad) ;
         
         pres_p2  = CPU_GetPressure( Array3D[DENS][k+2][j][i], Array3D[MOMX][k+2][j][i], 
                                     Array3D[MOMY][k+2][j][i], Array3D[MOMZ][k+2][j][i], 
                                     Array3D[ENGY][k+2][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
                                 
         dens_p1  = Array3D[DENS][k+1][j][i];
         pot_grad = PotArray3D[0][k+2][j][i] - PotArray3D[0][k][j][i];
         
         pres_bc  = pres_p2 + dens_p1*pot_grad - dens_p1*star_g*(2.0*dh[2]) ;
         pres_bc  = FMAX( pres_bc, MIN_PRES ) ;
         
         // use isothermal condition to determine rho
         rho_guess = pres_bc/RT_ref;         
         if ( rho_guess < rho_ratio_limit*rho_ref ) Array3D[DENS][k][j][i] = rho_ratio_limit*rho_ref ;
         else                                       Array3D[DENS][k][j][i] = rho_guess ;
         
         pres_bc = Array3D[DENS][k][j][i] * RT_ref;
         pres_bc  = FMAX( pres_bc, MIN_PRES ) ;
         
         Array3D[MOMX][k][j][i] = Array3D[DENS][k][j][i] * Vx_ref;
         Array3D[MOMY][k][j][i] = Array3D[DENS][k][j][i] * Vy_ref;
         Array3D[MOMZ][k][j][i] = Array3D[DENS][k][j][i] * Vz_ref;
         
         engy_bc  = pres_bc/Gamma_m1 + 0.5*( SQR(Array3D[MOMX][k][j][i]) + SQR(Array3D[MOMY][k][j][i]) + 
                                             SQR(Array3D[MOMZ][k][j][i]) ) / Array3D[DENS][k][j][i] ;
         
         Array3D[ENGY][k][j][i] = engy_bc ;
         
         
#        ifdef SUPPORT_GRACKLE
         if (GRACKLE_PRIMORDIAL != GRACKLE_PRI_CHE_NSPE9)
            Aux_Message(stderr, "User defined BC curretly only support Grackle with 9 spices! \n");
         
         for (int n=Idx_e; n<=Idx_H2II; n++) {
            Array3D[n][k][j][i] = Array3D[n][k_ref][j][i] * ( Array3D[DENS][k][j][i]*_rho_ref );
         }
#        endif //#SUPPORT_GRACKLE
         
#        if ( DUAL_ENERGY == DE_ENPY )
         Array3D[Idx_Enpy][k][j][i] = 
            CPU_Fluid2Entropy(Array3D[DENS][k][j][i], Array3D[MOMX][k][j][i], Array3D[MOMY][k][j][i],
                              Array3D[MOMZ][k][j][i], Array3D[ENGY][k][j][i], Gamma_m1);
#        endif // DUAL_ENERGY == DE_ENPY
         
      }         
   } // for (i, j)
#  endif // UserPotBC
   
}

//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_zp
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. update Array 
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  Array    : Fluid field to be set
//
//-------------------------------------------------------------------------------------------------------              
void BC_User_zp( real *Array, real *PotArray, const int NVar_Flu, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                 const int TFluVarIdxList[], const double dh[], const double *Corner, const int TVar )
{
#  ifdef UserPotBC
// 1D array -> 3D array
   const int FACE = 5;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX]    = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;

   
   const double X0    = Corner[0] + (double)Idx_Start[0]*dh[0];
   const double Y0    = Corner[1] + (double)Idx_Start[1]*dh[1];
   const double Z0    = Corner[2] + (double)Idx_Start[2]*dh[2];
   const double GM    = ExtAcc_AuxArray[3] ;
   const double X_ref = X0;
   const double Y_ref = Y0;
   const double Z_ref = Z0 - dh[2];
   const int    k_ref = Idx_Start[2]-1 ;
   
   const bool CheckMinPres_Yes = true;
   const real Gamma_m1         = GAMMA - (real)1.0;
   
   double pres_bc, engy_bc, sph_rad, star_g, dens_m1, pres_m2, pot_grad;
   double pres_ref, rho_guess, rho_ref, _rho_ref, RT_ref, Vx_ref, Vy_ref, Vz_ref;
   double X, Y, Z, Z_m1;
   int    i, j, k;
   
   if (PotArray != NULL) {
      Poi_BoundaryCondition_Extrapolation( PotArray, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
                                           Idx_Start, Idx_End, NULL, NULL );
   }
   
   real (*PotArray3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )PotArray;
   
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
   {  
      rho_ref   = Array3D[DENS][k_ref][j][i] ;
      _rho_ref  = 1 / rho_ref;
      Vx_ref    = Array3D[MOMX][k_ref][j][i] * _rho_ref;
      Vy_ref    = Array3D[MOMY][k_ref][j][i] * _rho_ref;
      Vz_ref    = Array3D[MOMZ][k_ref][j][i] * _rho_ref;
      pres_ref  = CPU_GetPressure( Array3D[DENS][k_ref][j][i], Array3D[MOMX][k_ref][j][i], 
                                   Array3D[MOMY][k_ref][j][i], Array3D[MOMZ][k_ref][j][i], 
                                   Array3D[ENGY][k_ref][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      RT_ref    = pres_ref * _rho_ref;
      
      if (!Aux_IsFinite(Vx_ref) || !Aux_IsFinite(Vy_ref) || !Aux_IsFinite(Vz_ref)) {
         Aux_Message(stdout, "Unphysical ref velocity in BC_USER ZP. \n");
      }
      
      
      // fill in bc
      for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2]; k++, Z+=dh[2]) {
         
         // calculate for BC value
         Z_m1     = Z-dh[2];
         sph_rad  = SQRT( Z_m1*Z_m1 + X*X );
         star_g   = - GM * Z_m1 / CUBE(sph_rad) ;
         
         pres_m2  = CPU_GetPressure( Array3D[DENS][k-2][j][i], Array3D[MOMX][k-2][j][i], 
                                     Array3D[MOMY][k-2][j][i], Array3D[MOMZ][k-2][j][i], 
                                     Array3D[ENGY][k-2][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
                                 
         dens_m1  = Array3D[DENS][k-1][j][i];
         pot_grad = PotArray3D[0][k][j][i] - PotArray3D[0][k-2][j][i];
         
         pres_bc  = pres_m2 - dens_m1*pot_grad + dens_m1*star_g*(2.0*dh[2]) ;
         pres_bc  = FMAX( pres_bc, MIN_PRES ) ;
         
         // use isothermal condition to determine rho
         rho_guess = pres_bc/RT_ref;         
         if ( rho_guess < rho_ratio_limit*rho_ref ) Array3D[DENS][k][j][i] = rho_ratio_limit*rho_ref ;
         else                                       Array3D[DENS][k][j][i] = rho_guess ;
         
         pres_bc = Array3D[DENS][k][j][i] * RT_ref;
         pres_bc  = FMAX( pres_bc, MIN_PRES ) ;
         
         Array3D[MOMX][k][j][i] = Array3D[DENS][k][j][i] * Vx_ref;
         Array3D[MOMY][k][j][i] = Array3D[DENS][k][j][i] * Vy_ref;
         Array3D[MOMZ][k][j][i] = Array3D[DENS][k][j][i] * Vz_ref;
         
         engy_bc  = pres_bc/Gamma_m1 + 0.5*( SQR(Array3D[MOMX][k][j][i]) + SQR(Array3D[MOMY][k][j][i]) + 
                                             SQR(Array3D[MOMZ][k][j][i]) ) / Array3D[DENS][k][j][i] ;
         
         Array3D[ENGY][k][j][i] = engy_bc ;
         
         
#        ifdef SUPPORT_GRACKLE
         if (GRACKLE_PRIMORDIAL != GRACKLE_PRI_CHE_NSPE9)
            Aux_Message(stderr, "User defined BC curretly only support Grackle with 9 spices! \n");
         
         for (int n=Idx_e; n<=Idx_H2II; n++) {
            Array3D[n][k][j][i] = Array3D[n][k_ref][j][i] * ( Array3D[DENS][k][j][i]*_rho_ref );
         }
#        endif //#SUPPORT_GRACKLE
         
#        if ( DUAL_ENERGY == DE_ENPY )
         Array3D[Idx_Enpy][k][j][i] = 
            CPU_Fluid2Entropy(Array3D[DENS][k][j][i], Array3D[MOMX][k][j][i], Array3D[MOMY][k][j][i],
                              Array3D[MOMZ][k][j][i], Array3D[ENGY][k][j][i], Gamma_m1);
#        endif // DUAL_ENERGY == DE_ENPY
         
      }         
   } // for (i, j)
   
#  endif // ifdef UserPotBC
}

#endif // COORDINATE == CYLINDRICAL

