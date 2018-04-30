#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Coord_CellIdx2AdoptedCoord
// Description :  Convert the input AMR level, patch ID, and cell index to the coordinate along the target
//                direction of the adopted coordinate system
//
// Note        :  1. Correspondence between different "directions" and "coordinates" in various coordinate systems
//                      Cartesian  : (0/1/2) <--> (x/y/z)
//                      Cylindrical: (0/1/2) <--> (r/phi/z)
//                      Spherical  : (0/1/2) <--> (r/theta/phi)
//                2. Only work on one target direction at a time
//                3. Always work on double precision
//
// Parameter   :  lv  : Target AMR level
//                PID : Target patch ID
//                dim : Target direction
//                idx : Target cell index along the target direction
//
// Return      :  Cartesian  : return x/y/z
//                Cylindrical: return r/phi/z
//                Spherical  : return r/theta/phi
//-------------------------------------------------------------------------------------------------------
double Aux_Coord_CellIdx2AdoptedCoord( const int lv, const int PID, const int dim, const int idx )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "lv = %d lies outside the accepted range !!\n", lv );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "PID = %d lies outside the accepted range (lv %d, NPatch %d) !!\n", PID, lv, amr->num[lv] );
#  endif


   return amr->patch[0][lv][PID]->EdgeL[dim] + (idx+0.5)*amr->dh[lv][dim];

} // FUNCTION : Aux_Coord_CellIdx2AdoptedCoord



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Coord_CellIdx2CartesianCoord
// Description :  Convert the input AMR level, patch ID, and cell index to the Cartesian coordinates
//
// Note        :  1. Correspondence between different "directions" and "coordinates" in various coordinate systems
//                      Cartesian  : (0/1/2) <--> (x/y/z)
//                      Cylindrical: (0/1/2) <--> (r/phi/z)
//                      Spherical  : (0/1/2) <--> (r/theta/phi)
//                2. Always work on double precision
//
// Parameter   :  lv    : Target AMR level
//                PID   : Target patch ID
//                i/j/k : Target cell index along each direction
//                xyz   : Cartesian coordinates to be returned
//
// Return      :  xyz[] -- Cartesian coordinates of the input cell
//-------------------------------------------------------------------------------------------------------
void Aux_Coord_CellIdx2CartesianCoord( const int lv, const int PID, const int i, const int j, const int k, double xyz[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "lv = %d lies outside the accepted range !!\n", lv );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "PID = %d lies outside the accepted range (lv %d, NPatch %d) !!\n", PID, lv, amr->num[lv] );
#  endif


   double XYZ[3];

   XYZ[0] = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i );
   XYZ[1] = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 1, j );
   XYZ[2] = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 2, k );

   Aux_Coord_Adopted2CartesianCoord( XYZ, xyz );

} // FUNCTION : Aux_Coord_CellIdx2CartesianCoord



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Coord_CellIdx2Volume
// Description :  Convert the input AMR level, patch ID, and cell index to the cell volume in the adopted
//                coordinate system
//
// Note        :  1. Cartesian   : dx*dy*dz
//                   Cylindrical : 0.5*(r_right^2 - r_left^2)*dphi*dz
//                   Spherical   : (1/3)*(r_right^3 - r_left^3)*(cos(theta_left) - cos(theta_right))*dphi
//                2. Always work on double precision
//
// Parameter   :  lv    : Target AMR level
//                PID   : Target patch ID
//                i/j/k : Target cell index along each direction
//
// Return      :  Volume of the target cell
//-------------------------------------------------------------------------------------------------------
double Aux_Coord_CellIdx2Volume( const int lv, const int PID, const int i, const int j, const int k )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "lv = %d lies outside the accepted range !!\n", lv );

#  if ( COORDINATE != CARTESIAN )
   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "PID = %d lies outside the accepted range (lv %d, NPatch %d) !!\n", PID, lv, amr->num[lv] );
#  endif
#  endif


   double dv;

#  if   ( COORDINATE == CARTESIAN )
   dv = amr->dh[lv][0]*amr->dh[lv][1]*amr->dh[lv][2];

#  elif ( COORDINATE == CYLINDRICAL )
   const double dr   = amr->dh[lv][0];
   const double dphi = amr->dh[lv][1];
   const double dz   = amr->dh[lv][2];
   const double r    = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i );
   const double r_R  = r + 0.5*dr;
   const double r_L  = r - 0.5*dr;

   dv = 0.5*( SQR(r_R) - SQR(r_L) )*dphi*dz;

#  elif ( COORDINATE == SPHERICAL )
   const double dr      = amr->dh[lv][0];
   const double dtheta  = amr->dh[lv][1];
   const double dphi    = amr->dh[lv][2];
   const double r       = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i );
   const double r_R     = r + 0.5*dr;
   const double r_L     = r - 0.5*dr;
   const double theta   = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 1, j );
   const double theta_R = theta + 0.5*dtheta;
   const double theta_L = theta - 0.5*dtheta;

   dv = 1.0/3.0*( CUBE(r_R) - CUBE(r_L) )*( cos(theta_L) - cos(theta_R) )*dphi;

#  else
#  error : UNSUPPORTED COORDINATE
#  endif // COORDINATE


   return dv;

} // FUNCTION : Aux_Coord_CellIdx2Volume



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Coord_Adopted2CartesianCoord
// Description :  Convert from the adopted coordinate system to the Cartesian coordinates
//
// Note        :  1. Correspondence between different "directions" and "coordinates" in various coordinate systems
//                      Cartesian  : (0/1/2) <--> (x/y/z)
//                      Cylindrical: (0/1/2) <--> (r/phi/z)
//                      Spherical  : (0/1/2) <--> (r/theta/phi)
//                2. Always work on double precision
//
// Parameter   :  in  : Input coordinates
//                out : Output coordinates
//
// Return      :  out[] -- Cartesian coordinates of the input coordinates
//-------------------------------------------------------------------------------------------------------
void Aux_Coord_Adopted2CartesianCoord( const double in[], double out[] )
{

#  if   ( COORDINATE == CARTESIAN )
   out[0] = in[0];
   out[1] = in[1];
   out[2] = in[2];

#  elif ( COORDINATE == CYLINDRICAL )
   out[0] = in[0]*cos( in[1] );
   out[1] = in[0]*sin( in[1] );
   out[2] = in[2];

#  elif ( COORDINATE == SPHERICAL )
   out[0] = in[0]*sin( in[1] )*cos( in[2] );
   out[1] = in[0]*sin( in[1] )*sin( in[2] );
   out[2] = in[0]*cos( in[1] );

#  else
#  error : UNSUPPORTED COORDINATE
#  endif // COORDINATE

} // FUNCTION : Aux_Coord_Adopted2CartesianCoord



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Coord_Cartesian2AdoptedCoord
// Description :  Convert from the Cartesian coordinates to the adopted coordinate system
//
// Note        :  1. Correspondence between different "directions" and "coordinates" in various coordinate systems
//                      Cartesian  : (0/1/2) <--> (x/y/z)
//                      Cylindrical: (0/1/2) <--> (r/phi/z)
//                      Spherical  : (0/1/2) <--> (r/theta/phi)
//                2. Always work on double precision
//
// Parameter   :  in  : Input coordinates
//                out : Output coordinates
//
// Return      :  out[] -- Adopted coordinate system of the input Cartesian coordinates
//-------------------------------------------------------------------------------------------------------
void Aux_Coord_Cartesian2AdoptedCoord( const double in[], double out[] )
{

#  if   ( COORDINATE == CARTESIAN )
   out[0] = in[0];
   out[1] = in[1];
   out[2] = in[2];

#  elif ( COORDINATE == CYLINDRICAL )
   out[0] = sqrt( SQR(in[0]) + SQR(in[1]) );
   out[1] = atan2( in[1], in[0] );           // note that atan2() returns in the range [-pi, pi]
   out[2] = in[2];

   if ( out[1] < 0.0 )  out[1] += 2.0*M_PI;  // define phi in the range [0, 2*pi]

#  elif ( COORDINATE == SPHERICAL )
   out[0] = sqrt( SQR(in[0]) + SQR(in[1]) + SQR(in[2]) );
   out[1] = acos( in[2]/out[0] );
   out[2] = atan2( in[1], in[0] );           // note that atan2() returns in the range [-pi, pi]

   if ( out[2] < 0.0 )  out[2] += 2.0*M_PI;  // define phi in the range [0, 2*pi]

#  else
#  error : UNSUPPORTED COORDINATE
#  endif // COORDINATE

} // FUNCTION : Aux_Coord_Cartesian2AdoptedCoord
