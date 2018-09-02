#include "GAMER.h"

static void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                       const int ii, const int jj, const int kk );




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Part
// Description :  Output part of data in the ASCII form
//
// Parameter   :  Part     : OUTPUT_XY   : XY plane
//                           OUTPUT_YZ   : YZ plane
//                           OUTPUT_XZ   : XZ plane
//                           OUTPUT_X    : X  line
//                           OUTPUT_Y    : Y  line
//                           OUTPUT_Z    : Z  line
//                           OUTPUT_DIAG : diagonal along (+1,+1,+1)
//
//                BaseOnly : Only output the base-level data
//
//                X        : X coordinate in the adopted coordinate system
//                Y        : Y coordinate in the adopted coordinate system
//                Z        : Z coordinate in the adopted coordinate system
//
//                FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Part( const OptOutputPart_t Part, const bool BaseOnly, const double X, const double Y,
                           const double Z, const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the input parameters
   if ( Part != OUTPUT_XY  &&  Part != OUTPUT_YZ  &&  Part != OUTPUT_XZ  &&
        Part != OUTPUT_X   &&  Part != OUTPUT_Y   &&  Part != OUTPUT_Z   &&  Part != OUTPUT_DIAG )
      Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [0 ~ 6] !!\n", Part );

   if (  ( Part == OUTPUT_YZ  ||  Part == OUTPUT_Y  ||  Part == OUTPUT_Z )  &&
         ( X < amr->BoxEdgeL[0]  ||  X >= amr->BoxEdgeR[0] )  )
      Aux_Error( ERROR_INFO, "incorrect X (out of range [%lf<=X<%lf]) !!\n", amr->BoxEdgeL[0], amr->BoxEdgeR[0] );

   if (  ( Part == OUTPUT_XZ  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Z )  &&
         ( Y < amr->BoxEdgeL[1]  ||  Y >= amr->BoxEdgeR[1] )  )
      Aux_Error( ERROR_INFO, "incorrect Y (out of range [%lf<=Y<%lf]) !!\n", amr->BoxEdgeL[1], amr->BoxEdgeR[1] );

   if (  ( Part == OUTPUT_XY  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Y )  &&
         ( Z < amr->BoxEdgeL[2]  ||  Z >= amr->BoxEdgeR[2] )  )
      Aux_Error( ERROR_INFO, "incorrect Z (out of range [%lf<=Z<%lf]) !!\n", amr->BoxEdgeL[2], amr->BoxEdgeR[2] );

   if ( Part == OUTPUT_DIAG )
   {
#     if ( COORDINATE == CARTESIAN )
      if (  !Mis_CompareRealValue( amr->BoxSize[0], amr->BoxSize[1], NULL, false ) ||
            !Mis_CompareRealValue( amr->BoxSize[0], amr->BoxSize[2], NULL, false )  )
         Aux_Error( ERROR_INFO, "simulation domain must be cubic for \"OUTPUT_DIAG\" !!\n" );

      if (  !Mis_CompareRealValue( amr->dh[0][0], amr->dh[0][1], NULL, false )  ||
            !Mis_CompareRealValue( amr->dh[0][0], amr->dh[0][2], NULL, false )    )
         Aux_Error( ERROR_INFO, "currently the Cartesian coordinates only work with cubic cells --> dh[lv=0] = (%20.14e, %20.14e, %20.14e) !!\n",
                    amr->dh[0][0], amr->dh[0][1], amr->dh[0][2] );
#     else
      Aux_Error( ERROR_INFO, "non-Cartesian coordinates do not support %s() !!\n", __FUNCTION__ );
#     endif
   }


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// check if the file already exists
   if ( MPI_Rank == 0 )
   {
      if ( Aux_CheckFileExist(FileName) )
      {
         Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

         FILE *Temp = fopen( FileName, "w" );
         fclose( Temp );
      }
   }


   const int NLv = ( BaseOnly ) ? 1 : NLEVEL;

   int     ii, jj, kk, scale;
   double *dh, XX, YY, ZZ;    // XX,YY,ZZ => physical coordinates of cell left edge
   int    *Corner  = NULL;    // patch corner in scale
   double *EdgeL   = NULL;    // patch corner in physical coord.
   double *EdgeR   = NULL;
   bool    Check_X = false;
   bool    Check_Y = false;
   bool    Check_Z = false;

   switch ( Part )
   {
      case OUTPUT_XY :                                      Check_Z = true;   break;
      case OUTPUT_YZ :  Check_X = true;                                       break;
      case OUTPUT_XZ :                    Check_Y = true;                     break;
      case OUTPUT_X  :                    Check_Y = true;   Check_Z = true;   break;
      case OUTPUT_Y  :  Check_X = true;                     Check_Z = true;   break;
      case OUTPUT_Z  :  Check_X = true;   Check_Y = true;                     break;

      case OUTPUT_DIAG :
      case OUTPUT_PART_NONE : break; // do nothing
   }


   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         FILE *File = fopen( FileName, "a" );

//       output header
         if ( TargetMPIRank == 0 )
         {
            fprintf( File, "#%10s %10s %10s %20s %20s %20s", "i", "j", "k", "X", "Y", "Z" );

            for (int v=0; v<NCOMP_TOTAL; v++)
            fprintf( File, "%14s", FieldLabel[v] );

#           ifdef GRAVITY
            if ( OPT__OUTPUT_POT )
            fprintf( File, "%14s", PotLabel );
#           endif

//          other derived fields
#           if ( MODEL == HYDRO  ||  MODEL == MHD )
            fprintf( File, "%14s", "Pressure" );
#           endif

            fprintf( File, "\n" );
         } // if ( TargetMPIRank == 0 )


//       output data
         for (int lv=0; lv<NLv; lv++)
         {
            dh    = amr->dh   [lv];
            scale = amr->scale[lv];

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
//             output the patch data only if it has no son (if the option "BaseOnly" is turned off)
               if ( amr->patch[0][lv][PID]->son == -1  ||  BaseOnly )
               {
                  Corner = amr->patch[0][lv][PID]->corner;
                  EdgeL  = amr->patch[0][lv][PID]->EdgeL;
                  EdgeR  = amr->patch[0][lv][PID]->EdgeR;

                  if ( Part == OUTPUT_DIAG ) // (+1,+1,+1) diagonal
                  {
                     if ( Corner[0] == Corner[1]  &&  Corner[0] == Corner[2] )
                     {
                        for (int k=0; k<PS1; k++)
                        {
                           kk = Corner[2] + k*scale;

                           WriteFile( File, lv, PID, k, k, k, kk, kk, kk );
                        }
                     }
                  } // if ( Part == OUTPUT_DIAG )


                  else // x/y/z lines || xy/yz/xz slices
                  {
//                   check whether the patch corner is within the target range
                     if (  !Check_X  ||  ( EdgeL[0]<=X && EdgeR[0]>X )  )
                     if (  !Check_Y  ||  ( EdgeL[1]<=Y && EdgeR[1]>Y )  )
                     if (  !Check_Z  ||  ( EdgeL[2]<=Z && EdgeR[2]>Z )  )
                     {
//                      check whether the cell is within the target range
                        for (int k=0; k<PS1; k++)  {  kk = Corner[2] + k*scale;  ZZ = EdgeL[2] + k*dh[2];
                                                      if ( Check_Z && ( ZZ>Z || ZZ+dh[2]<=Z ) )    continue;

                        for (int j=0; j<PS1; j++)  {  jj = Corner[1] + j*scale;  YY = EdgeL[1] + j*dh[1];
                                                      if ( Check_Y && ( YY>Y || YY+dh[1]<=Y ) )    continue;

                        for (int i=0; i<PS1; i++)  {  ii = Corner[0] + i*scale;  XX = EdgeL[0] + i*dh[0];
                                                      if ( Check_X && ( XX>X || XX+dh[0]<=X ) )    continue;

                           WriteFile( File, lv, PID, i, j, k, ii, jj, kk );

                        }}}
                     } // if patch corner is within the target range

                  } // if ( Part == OUTPUT_DIAG ... else ... )
               } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // for (int lv=0; lv<NLv; lv++)

         fclose( File );

      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Part



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteFile
// Description :  Output data to file
//
// Parameter   :  File     : File pointer
//                lv       : Target refinement level
//                PID      : Patch ID
//                i/j/k    : Cell indices within the patch
//                ii/jj/kk : Cell scale indices in the simulation domain
//-------------------------------------------------------------------------------------------------------
void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                const int ii, const int jj, const int kk )
{

   real u[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)   u[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

// output cell indices and coordinates
   fprintf( File, " %10d %10d %10d %20.14e %20.14e %20.14e",
            ii, jj, kk, Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i ),
                        Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 1, j ),
                        Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 2, k ) );

// output all variables in the fluid array
   for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, " %13.6e", u[v] );

// output potential
#  ifdef GRAVITY
   if ( OPT__OUTPUT_POT )
   fprintf( File, " %13.6e", amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i] );
#  endif

// output other derived fields
#  if   ( MODEL == HYDRO )
   const bool CheckMinPres_Yes = true;
   fprintf( File, " %13.6e", CPU_GetPressure(u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], GAMMA-1.0, CheckMinPres_Yes, MIN_PRES) );
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL

   fprintf( File, "\n" );

} // FUNCTION : WriteFile
