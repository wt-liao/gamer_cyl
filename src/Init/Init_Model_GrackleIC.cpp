#  if defined (SUPPORT_GRACKLE) && defined (MODEL_IC_GRACKLE)

#include "GAMER.h"
#include "TestProb.h"

// some constants in cgs
static double m_H_cgs     = 1.6733e-24   ;
static double m_e_cgs     = 9.1093897e-28;
static double m_He_cgs    = 6.64647356e-24;
static double k_B         = 1.380658e-16 ;
static double planck_h    = 6.6260755e-27;
static double eV_cgs      = 1.6021772e-12;
static double engy_HII    = 13.54*eV_cgs ; // eV
static double engy_HeII   = 24.48*eV_cgs ; 
static double engy_HeIII  = 54.17*eV_cgs ;

static void Find_GrackleIC(const real* PriVar, double* n_grackle, const double X, const double Y, const double Z);
static void Con2Pri(const real* ConVar, real* PriVar);

static double rate_k4_func(const real T);
static double rate_k5_func(const real T);
static double n_HII_func(const real n_cgs, const real T);
static double n_HeII_func(const real n_cgs, const real T);
static double n_H2_func(const real n_cgs, const real T, const real tau_dyn_cgs);

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GrackleField
// Description :  Work for Init_GAMER
//
// Note        :  1.
//                
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------

void Init_GrackleField(){
   if ( MPI_Rank == 0 ) Aux_Message(stdout, "Generating IC for Grackle... \n");
   if ( MPI_Rank == 0 ) Aux_Message(stdout, "Grackle unit: [L]=%12.8e; [T]=%12.8e; [Dens]=%12.8e . \n", 
                                    Che_Units.length_units, Che_Units.time_units, Che_Units.density_units);
   
   const int lv = 0; // currently only implement for no-amr scheme
   real X, Y, Z;
   real  ConVar[NCOMP_TOTAL], PriVar[NCOMP_TOTAL] ;
   real (*fluid)[PS1][PS1][PS1]=NULL;
   
//#     pragma omp parallel for private( fluid, ConVAr, PriVar, X, Y, Z ) schedule( runtime ) num_threads( OMP_NThread )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++) {
      
      const bool NormPassive_No  = false;
      const bool JeansMinPres_No = false;
      double n_grackle[12]; 
      
      fluid = amr->patch[amr->FluSg[lv]][lv][PID]->fluid;
      
      for (int k=0; k<PS1; k++)  {  Z = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 2, k );
      for (int j=0; j<PS1; j++)  {  Y = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 1, j );
      for (int i=0; i<PS1; i++)  {  X = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i );
         // get fluid field
         for (int v=0; v<NCOMP_TOTAL; v++)   ConVar[v] = fluid[v][k][j][i];
         Con2Pri(ConVar, PriVar);
         
         // calculate the necesary field
         Find_GrackleIC(PriVar, n_grackle, X, Y, Z);
         
         // chemistry field
         if ( GRACKLE_PRIMORDIAL == GRACKLE_PRI_CHE_NSPE6 )
            Aux_Message(stderr, "In <%s>, Grackle IC model currently only support species == 9! \n", __FUNCTION__);
         
         if ( GRACKLE_PRIMORDIAL == GRACKLE_PRI_CHE_NSPE9 ) {
            
            fluid[Idx_e    ][k][j][i] = n_grackle[0];
            fluid[Idx_HI   ][k][j][i] = n_grackle[1];
            fluid[Idx_HII  ][k][j][i] = n_grackle[2];
            fluid[Idx_HeI  ][k][j][i] = n_grackle[3];
            fluid[Idx_HeII ][k][j][i] = n_grackle[4];
            fluid[Idx_HeIII][k][j][i] = n_grackle[5];
                  
            fluid[Idx_HM   ][k][j][i] = n_grackle[6];
            fluid[Idx_H2I  ][k][j][i] = n_grackle[7];
            fluid[Idx_H2II ][k][j][i] = n_grackle[8];
            
         }
         
         if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
            Aux_Message(stderr, "In <%s>, Grackle IC model currently only support species == 9! \n", __FUNCTION__);
            fluid[Idx_DI   ][k][j][i] = 0.0;
            fluid[Idx_DII  ][k][j][i] = 0.0;
            fluid[Idx_HDI  ][k][j][i] = 0.0;
         }
         
      }}} // for loop: k, j, i
   } // openmp black
   
   if ( MPI_Rank == 0 ) Aux_Message(stdout, "Generating IC for Grackle... Done! \n");
}  // Init_GrackleField()



//-------------------------------------------------------------------------------------------------------
// Function    :  Find_GrackleIC
// Description :  
//
// Note        :  1.
//                
// Parameter   :  PriVar   : Primary Variable
//                n_grackle: chemistry density field; grackle uses mass density field 
//                           [0]: n_e, [1]: n_HI, [2]: n_HII, [3]: n_HeI, [4]: n_HeII, [5]: n_HeIII
//                           [6]: n_HM, [7]: n_H2, [8]: n_H2II
//-------------------------------------------------------------------------------------------------------

void Find_GrackleIC(const real* PriVar, double* n_grackle, const double X, const double Y, const double Z){
   
   const double rho_unit    = Che_Units.density_units;
   const double L_unit      = Che_Units.length_units;
   const double time_unit   = Che_Units.time_units ;
   const double mass_unit   = rho_unit * POW(L_unit, 3);
   const double _rho_unit   = 1/rho_unit;
   
   //### check if this fraction is mass fraction; should be! 
   const double y_He        = 0.24; // He fraction
   const double H2_FLOOR    = 1e-8; // probably don't need, since the relaxation step will allow the formation of H_2
   const double m_ave_cgs   = m_H_cgs*(1-y_He) + m_He_cgs*y_He;
   
   // dimensionless field
   const double GAMMA_1     = GAMMA - (real)1.0 ;
   const double m_H         = m_H_cgs / mass_unit;
   const double R           = (k_B/m_ave_cgs) * SQR(time_unit/L_unit) ;
   const double rho         = PriVar[DENS]; 
   const double pres        = PriVar[ENGY];
   const double T           = pres / ( rho * R ) ; // Temperature unit is K 
   
   // calculate for tau_dyn
   const double GM          = ExtAcc_AuxArray[3] ;
   const double sph_rad     = SQRT( X*X + Z*Z );
   const double tau_dyn     = SQRT( CUBE(sph_rad) / GM );
   const double tau_dyn_cgs = tau_dyn * time_unit;
   
   // cgs field
   const double rho_cgs     = rho*rho_unit ; 
   const double rho_H_cgs   = rho_cgs * (1-y_He);
   const double rho_He_cgs  = rho_cgs * y_He;
   
   // y_He is the mass fraction, not particle density fraction
   const double n_H_cgs     = rho_H_cgs /(1.0*m_H_cgs ); // 
   const double n_He_cgs    = rho_He_cgs/(1.0*m_He_cgs);
   const double n_cgs       = n_H_cgs + n_He_cgs ;
   
   // find all field data; f: mass fraction
   if ( GRACKLE_PRIMORDIAL == GRACKLE_PRI_CHE_NSPE9 ) {
      // n: particle number density; f: mass fraction
      const double n_HII   = n_HII_func(n_H_cgs, T);
      const double n_H2I   = n_H2_func(0.5*n_H_cgs, T, tau_dyn_cgs) ;
      const double n_HI    = n_H_cgs - n_HII - 2*n_H2I;
      const double n_HM    = TINY_NUMBER; // in grackle, HM and H2II are both in equlibrium
      const double n_H2II  = TINY_NUMBER;
      
      const double n_HeII  = n_HeII_func(n_He_cgs, T); 
      const double n_HeI   = n_He_cgs - n_HeII ;   
      const double n_HeIII = TINY_NUMBER;
      
      const double n_e     = n_HII + n_HeII;
      
      // [0]: n_e, [1]: n_HI, [2]: n_HII, [3]: n_HeI, [4]: n_HeII, [5]: n_HeIII
      n_grackle[0] = n_e     * m_H_cgs  * _rho_unit ;
      n_grackle[1] = n_HI    * m_H_cgs  * _rho_unit ;
      n_grackle[2] = n_HII   * m_H_cgs  * _rho_unit ;
      n_grackle[3] = n_HeI   * m_He_cgs * _rho_unit ;
      n_grackle[4] = n_HeII  * m_He_cgs * _rho_unit ;
      n_grackle[5] = n_HeIII * m_He_cgs * _rho_unit ;
      
      // [6]: n_HM, [7]: n_H2I, [8]: n_H2II
      n_grackle[6] = n_HM    * m_H_cgs     * _rho_unit ;
      n_grackle[7] = n_H2I   * (2*m_H_cgs) * _rho_unit ;
      n_grackle[8] = n_H2II  * (2*m_H_cgs) * _rho_unit ;
      
   }
   
   else 
      Aux_Message(stderr, "In <%s>, Grackle IC model currently only support species == 9! \n", __FUNCTION__);

}  // Find_GrackleIC()


// POPIII
//-------------------------------------------------------------------------------------------------------
// Function    :  rate_k4 in cgs
//-------------------------------------------------------------------------------------------------------
double rate_k4_func(const real T) {
   double rate_k4 = 5.5e-29 * POW(T, -1) ;
   
   return rate_k4 ;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  rate_k5 in cgs
//-------------------------------------------------------------------------------------------------------
double rate_k5_func(const real T) {
   double rate_k5 = 6.5e-7 * POW(T, -0.5) * EXP(-52000/T) * (1-EXP(-6000/T))  ;   
   
   return rate_k5 ;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  n_HII 
//-------------------------------------------------------------------------------------------------------
double n_HII_func(const real n_cgs, const real T) {
   
   const double consts  = 2*M_PI*m_e_cgs*k_B*T / SQR(planck_h) ;
   const double T_ratio = engy_HII / (k_B*T) ; 
   const double g_ratio = 1.0; 
   
   double c1 = g_ratio*POW(consts, 1.5) / n_cgs * EXP(- T_ratio); 
   double x  = 0.5*( -c1 + SQRT(c1*c1 + 4.0*c1) );
   
   double n_HII = n_cgs * x;
      
   return n_HII ;
}


//-------------------------------------------------------------------------------------------------------
// Function    :  n_HeII 
//-------------------------------------------------------------------------------------------------------
double n_HeII_func(const real n_cgs, const real T) {
   
   const double consts  = 2*M_PI*m_e_cgs*k_B*T / SQR(planck_h) ;
   const double T_ratio = engy_HeII / (k_B*T) ; 
   const double g_ratio = 4.0; 
   
   double c1 = g_ratio*POW(consts, 1.5) / n_cgs * EXP(- T_ratio); 
   double x  = 0.5*( -c1 + SQRT(c1*c1 + 4.0*c1) );
   
   double n_HeII = n_cgs * x;
      
   return n_HeII ;
}


//-------------------------------------------------------------------------------------------------------
// Function    :  f_H2
//-------------------------------------------------------------------------------------------------------
double n_H2_func(const real n_cgs, const real T, const real tau_dyn_cgs) {
   
   const double k4    = rate_k4_func(T);
   const double k5    = rate_k5_func(T); 
   const double coeff = k5/(2*k4);
   
   double n_HI, n_H2_estimate, n_H2, tau_chem;
   
   n_HI          = 0.5*( -coeff + SQRT(coeff*coeff + 4*coeff*n_cgs) ) ;
   n_H2_estimate = (k4/k5)*SQR(n_HI) ;
   
   tau_chem = 1 / ( k4*SQR(n_cgs) ) ;
   n_H2     = n_H2_estimate * EXP(-tau_chem / tau_dyn_cgs);

   //n_H2          = FMIN(n_H2_estimate, n_cgs) ;
   
   return n_H2 ;
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Con2Pri
//-------------------------------------------------------------------------------------------------------
void Con2Pri(const real* ConVar, real* PriVar) {
   real dens, _dens, vel_x, vel_y, vel_z, pres, Eint; 
   
   dens  = ConVar[DENS];
   _dens = 1/dens;
   vel_x = ConVar[MOMX] * _dens;
   vel_y = ConVar[MOMY] * _dens;
   vel_z = ConVar[MOMZ] * _dens;
   Eint  = ConVar[ENGY] - 0.5*dens*( SQR(vel_x)+SQR(vel_y)+SQR(vel_z) );
   pres  = FMAX( Eint/(GAMMA-1), MIN_PRES );
   
   PriVar[DENS] = dens;
   PriVar[MOMX] = vel_x;
   PriVar[MOMY] = vel_y;
   PriVar[MOMZ] = vel_z;
   PriVar[ENGY] = pres; 
}

#  endif // defined (SUPPORT_GRACKLE) && defined (MODEL_IC_GRACKLE)



