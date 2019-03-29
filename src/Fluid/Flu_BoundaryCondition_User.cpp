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
// 1D array -> 3D array
   const int FACE = 0;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   // extrapolate density
   //Poi_BoundaryCondition_Extrapolation( Array, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
   //                                     Idx_Start, Idx_End, NULL, NULL );
   
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
   
   double sph_rad, star_g, dens, pot_grad, pres_grad, vtheta_square, vtheta, pres;
   double rho_ref, pres_ref;
   
   double X, Y, Z;
   int    i, j, k;

   // for non-self-grvitating test
   pot_grad  = 0.0;
   pres_grad = 0.0;
   
   const double rho_0=1, T_0=1, R_0=1, const_R=1;
   const double slope_p=0, slope_q=0;

   for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2];   k++, Z+=dh[2])
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1];   j++, Y+=dh[1])
   {
      // fill in bc
      rho_ref   = Array3D[DENS][k][j][i_ref] ;
      pres_ref  = CPU_GetPressure( Array3D[DENS][k][j][i_ref], Array3D[MOMX][k][j][i_ref], 
                                   Array3D[MOMY][k][j][i_ref], Array3D[MOMZ][k][j][i_ref], 
                                   Array3D[ENGY][k][j][i_ref], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
                                   
      for (i=Idx_End[0],   X=X0; i>=Idx_Start[0]; i--, X-=dh[0])
      {
      
         // outflow 
         Array3D[DENS][k][j][i] = Array3D[DENS][k][j][i_ref] ;
         Array3D[MOMX][k][j][i] = 0; // Array3D[MOMX][k][j][i_ref] ;
         //Array3D[MOMY][k][j][i] = Array3D[MOMY][k][j][i_ref] ;
         Array3D[MOMZ][k][j][i] = 0; //Array3D[MOMZ][k][j][i_ref] ;
         //Array3D[ENGY][k][j][i] = Array3D[ENGY][k][j][i_ref] ;
      
         // derived other field
         sph_rad  = SQRT( Z*Z + X*X );
         star_g   = - GM * X / CUBE(sph_rad) ;
         dens     = Array3D[DENS][k][j][i]; 
      
         //### BC3: force balance for vtheta
         //vtheta_square = (pres_grad + dens*pot_grad - dens*star_g)*(X/dens); // <- only useful when self-g is active
         vtheta_square = -star_g*X;
         vtheta = (vtheta_square > 0.0)? SQRT(vtheta_square) : 0.0 ;
         Array3D[MOMY][k][j][i] = dens * vtheta ;
      
         //### BC2: no-shear BC
         //Array3D[MOMY][k][j][i] = Array3D[MOMY][k][j][i_ref] * (X_ref/X) ;
      
         // keep pressure the same, but modify K.E., since pres_grad = 0.0 currently
         pres = pres_ref; 
         Array3D[ENGY][k][j][i] = (pres/Gamma_m1) 
                                + 0.5*(SQR(Array3D[MOMX][k][j][i])+SQR(Array3D[MOMY][k][j][i])+SQR(Array3D[MOMZ][k][j][i])) /dens ; 
      
      
         //### unperturbed BC
         /*
         real R_norm      = X / R_0; 
         real rho_mid     = rho_0 * POW(R_norm, slope_p) ;
         real temperature = T_0 * POW(R_norm, slope_q) ;
         real cs_square   = const_R * temperature ; 
         real _sph_r      = 1.0 / SQRT(X*X + Z*Z);
         real _R_norm     = 1.0 / R_norm ;
         real omega_kep   = SQRT(GM/CUBE(X));
         real H           = SQRT(cs_square)/omega_kep;
   
         const real rho         = rho_mid * EXP( GM/cs_square * (_sph_r - 1/X) );
         const real omega       = omega_kep * SQRT(1.0 + (slope_q+slope_p)*SQR(H/X) + slope_q*(1.0-X*_sph_r) ) ;
         const real pressure    = rho*const_R*temperature ;
   
         Array3D[DENS][k][j][i] = rho;
         Array3D[MOMX][k][j][i] = 0.0;
         Array3D[MOMY][k][j][i] = rho* (X*omega) ;
         Array3D[MOMZ][k][j][i] = 0.0;
         Array3D[ENGY][k][j][i] = 0.5*rho*SQR(X*omega) + pressure/(GAMMA-1.0) ;
         */
      
      } // for i
      
   } // for k, j

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
// 1D array -> 3D array
   const int FACE = 1;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   // extrapolate density
   //Poi_BoundaryCondition_Extrapolation( Array, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
   //                                     Idx_Start, Idx_End, NULL, NULL );
   
   const double X0    = Corner[0] + (double)Idx_Start[0]*dh[0];
   const double Y0    = Corner[1] + (double)Idx_Start[1]*dh[1];
   const double Z0    = Corner[2] + (double)Idx_Start[2]*dh[2];
   const double GM    = ExtAcc_AuxArray[3] ;
   const double X_ref = X0 - dh[0];
   const double Y_ref = Y0;
   const double Z_ref = Z0;
   
   const int    i_ref = Idx_Start[0]-1 ;
   
   const bool CheckMinPres_Yes = true;
   const real Gamma_m1         = GAMMA - (real)1.0;
   
   double sph_rad, star_g, dens, pot_grad, pres_grad, vtheta_square, vtheta, pres;
   double rho_ref, pres_ref ;
   
   double X, Y, Z;
   int    i, j, k;
   
   // for non-self-grvitating test
   pot_grad  = 0.0;
   pres_grad = 0.0;
   
   const double rho_0=1, T_0=1, R_0=1, const_R=1;
   const double slope_p=0, slope_q=0;
   
   for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2]; k++, Z+=dh[2])
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   {
      // fill in bc
      rho_ref   = Array3D[DENS][k][j][i_ref] ;
      pres_ref  = CPU_GetPressure( Array3D[DENS][k][j][i_ref], Array3D[MOMX][k][j][i_ref], 
                                   Array3D[MOMY][k][j][i_ref], Array3D[MOMZ][k][j][i_ref], 
                                   Array3D[ENGY][k][j][i_ref], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      
      for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
      {
      
         // outflow 
         Array3D[DENS][k][j][i] = Array3D[DENS][k][j][i_ref] ;
         Array3D[MOMX][k][j][i] = 0; // Array3D[MOMX][k][j][i_ref] ;
         //Array3D[MOMY][k][j][i] = Array3D[MOMY][k][j][i_ref] ;
         Array3D[MOMZ][k][j][i] = 0; // Array3D[MOMZ][k][j][i_ref] ;
         //Array3D[ENGY][k][j][i] = Array3D[ENGY][k][j][i_ref] ;
      
         // derived other field
         sph_rad  = SQRT( Z*Z + X*X );
         star_g   = - GM * X / CUBE(sph_rad) ;
         dens      = Array3D[DENS][k][j][i];    
      
         //vtheta_square = (pres_grad + dens*pot_grad - dens*star_g)*(X/dens);
         vtheta_square = -star_g*X;
         vtheta = (vtheta_square > 0.0)? SQRT(vtheta_square) : 0.0 ;
         Array3D[MOMY][k][j][i] = dens * vtheta ;
      
         //Array3D[MOMY][k][j][i] = dens * vtheta ;
      
         // no shear BC
         //Array3D[MOMY][k][j][i] = Array3D[MOMY][k][j][i_ref] * (X_ref/X) ;
      
         // keep pressure the same, but modify K.E., since pres_grad = 0.0 currently
         pres = pres_ref; 
         Array3D[ENGY][k][j][i] = (pres/Gamma_m1) 
                                + 0.5*(SQR(Array3D[MOMX][k][j][i])+SQR(Array3D[MOMY][k][j][i])+SQR(Array3D[MOMZ][k][j][i]))/dens ; 
      
         /*
         real R_norm      = X / R_0; 
         real rho_mid     = rho_0 * POW(R_norm, slope_p) ;
         real temperature = T_0 * POW(R_norm, slope_q) ;
         real cs_square   = const_R * temperature ; 
         real _sph_r      = 1.0 / SQRT(X*X + Z*Z);
         real _R_norm     = 1.0 / R_norm ;
         real omega_kep   = SQRT(GM/CUBE(X));
         real H           = SQRT(cs_square)/omega_kep;
   
         const real rho         = rho_mid * EXP( GM/cs_square * (_sph_r - 1/X) );
         const real omega       = omega_kep * SQRT(1.0 + (slope_q+slope_p)*SQR(H/X) + slope_q*(1.0-X*_sph_r) ) ;
         const real pressure    = rho*const_R*temperature ;
   
         Array3D[DENS][k][j][i] = rho;
         Array3D[MOMX][k][j][i] = 0.0;
         Array3D[MOMY][k][j][i] = rho* (X*omega) ;
         Array3D[MOMZ][k][j][i] = 0.0;
         Array3D[ENGY][k][j][i] = 0.5*rho*SQR(X*omega) + pressure/(GAMMA-1.0) ;
         */
      
      } // for i
      
   } // for k, j   

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
// 1D array -> 3D array
   const int FACE = 4;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   // extrapolate density
   //Poi_BoundaryCondition_Extrapolation( Array, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
   //                                     Idx_Start, Idx_End, NULL, NULL );
   
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
   double rho_ref, pres_ref;
   double X, Y, Z, Z_p1;
   int    i, j, k;
   
   // for non-self-grvitating test
   pot_grad = 0.0;
   
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
   {  
      // fill in bc
      /*
      rho_ref   = Array3D[DENS][k_ref][j][i] ;
      pres_ref  = CPU_GetPressure( Array3D[DENS][k_ref][j][i], Array3D[MOMX][k_ref][j][i], 
                                   Array3D[MOMY][k_ref][j][i], Array3D[MOMZ][k_ref][j][i], 
                                   Array3D[ENGY][k_ref][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      */
                                   
      for (k=Idx_End[2], Z=Z0; k>=Idx_Start[2]; k--, Z-=dh[2]) {
         
         // outflow 
         Array3D[DENS][k][j][i] = Array3D[DENS][k_ref][j][i] ;
         Array3D[MOMX][k][j][i] = 0; //Array3D[MOMX][k_ref][j][i] ;
         Array3D[MOMY][k][j][i] = Array3D[MOMY][k_ref][j][i] ;
         Array3D[MOMZ][k][j][i] = 0; //Array3D[MOMZ][k_ref][j][i] ;
         
         // calculate for BC value
         Z_p1     = Z+dh[2];
         sph_rad  = SQRT( Z_p1*Z_p1 + X*X );
         star_g   = - GM * Z_p1 / CUBE(sph_rad) ;
         
         pres_p2  = CPU_GetPressure( Array3D[DENS][k+2][j][i], Array3D[MOMX][k+2][j][i], 
                                     Array3D[MOMY][k+2][j][i], Array3D[MOMZ][k+2][j][i], 
                                     Array3D[ENGY][k+2][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
                                 
         dens_p1  = Array3D[DENS][k+1][j][i];
         
         pres_bc  = pres_p2 - dens_p1*star_g*(2.0*dh[2]) ;
         pres_bc  = FMAX( pres_bc, MIN_PRES ) ;
         
         
         // isothermal bc
         //pres_bc  = pres_ref * (Array3D[DENS][k][j][i]/rho_ref);
         
         engy_bc  = pres_bc/Gamma_m1 + 0.5*( SQR(Array3D[MOMX][k][j][i]) + SQR(Array3D[MOMY][k][j][i]) + 
                                             SQR(Array3D[MOMZ][k][j][i]) ) / Array3D[DENS][k][j][i] ;
         
         Array3D[ENGY][k][j][i] = engy_bc ;
      }         
   } // for (i, j)
   
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
// 1D array -> 3D array
   const int FACE = 5;
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX]    = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;
   
   // extrapolate density
   //Poi_BoundaryCondition_Extrapolation( Array, FACE, 1, GhostSize, ArraySizeX, ArraySizeY, ArraySizeZ, 
   //                                     Idx_Start, Idx_End, NULL, NULL );

   
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
   
   double pres_bc, engy_bc, sph_rad, star_g, dens_m1, pres_m2, pot_grad ;
   double pres_ref, rho_ref;
   double X, Y, Z, Z_m1;
   int    i, j, k;
   
   // for non-self-grvitating test
   pot_grad = 0.0;
   
   for (j=Idx_Start[1], Y=Y0; j<=Idx_End[1]; j++, Y+=dh[1])
   for (i=Idx_Start[0], X=X0; i<=Idx_End[0]; i++, X+=dh[0])
   {  
      // fill in bc
      /*
      rho_ref   = Array3D[DENS][k_ref][j][i] ;
      pres_ref  = CPU_GetPressure( Array3D[DENS][k_ref][j][i], Array3D[MOMX][k_ref][j][i], 
                                   Array3D[MOMY][k_ref][j][i], Array3D[MOMZ][k_ref][j][i], 
                                   Array3D[ENGY][k_ref][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
      */
      
      for (k=Idx_Start[2], Z=Z0; k<=Idx_End[2]; k++, Z+=dh[2]) {
         
         // outflow 
         Array3D[DENS][k][j][i] = Array3D[DENS][k_ref][j][i] ;
         Array3D[MOMX][k][j][i] = 0; // Array3D[MOMX][k_ref][j][i] ;
         Array3D[MOMY][k][j][i] = Array3D[MOMY][k_ref][j][i] ;
         Array3D[MOMZ][k][j][i] = 0; // Array3D[MOMZ][k_ref][j][i] ;  
         
         
         // calculate for force balance BC value
         Z_m1     = Z-dh[2];
         sph_rad  = SQRT( Z_m1*Z_m1 + X*X );
         star_g   = - GM * Z_m1 / CUBE(sph_rad) ;
         
         pres_m2  = CPU_GetPressure( Array3D[DENS][k-2][j][i], Array3D[MOMX][k-2][j][i], 
                                     Array3D[MOMY][k-2][j][i], Array3D[MOMZ][k-2][j][i], 
                                     Array3D[ENGY][k-2][j][i], Gamma_m1, CheckMinPres_Yes, MIN_PRES );
                                 
         dens_m1  = Array3D[DENS][k-1][j][i];
         pres_bc  = pres_m2 + dens_m1*star_g*(2.0*dh[2]) ;
         pres_bc  = FMAX( pres_bc, MIN_PRES ) ;
         
         
         // isothermal bc
         //pres_bc  = pres_ref * (Array3D[DENS][k][j][i]/rho_ref);
         
         engy_bc  = pres_bc/Gamma_m1 + 0.5*( SQR(Array3D[MOMX][k][j][i]) + SQR(Array3D[MOMY][k][j][i]) + 
                                             SQR(Array3D[MOMZ][k][j][i]) ) / Array3D[DENS][k][j][i] ;
         
         Array3D[ENGY][k][j][i] = engy_bc ;
      }         
   } // for (i, j)
   
}

#endif // COORDINATE == CYLINDRICAL

