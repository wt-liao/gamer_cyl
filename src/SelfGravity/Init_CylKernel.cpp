#include "GAMER.h"

#if (COORDINATE == CYLINDRICAL)
#ifdef GRAVITY

extern rfftwnd_plan     FFTW_Plan, FFTW_Plan_Inv;    


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_CylKernel
// Description :  Calculate Kernel function in k-space
//-------------------------------------------------------------------------------------------------------
void Init_CylKernel(){
   
   const double *dh      = amr->dh[0]; 
   const double dh_cube  = dh[0]*dh[1]*dh[2] ;
   const int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2]*2 }; // FFT_Size[0] is redundunt

   double x, xp, y, z ; // (x, y, z) <-> (r, phi, z)
   int ID_planX, ID_planYZ, kk, ii, iip ;
      
# ifdef SERIAL
   Aux_Message(stderr, "Cylindrical Self-Gravity is not yet ready for serial mode! \n");
   
# else
   //### NEED a scheme to determine RANK_I_TOT, RANK_IP_TOT for general case 
      
   // MPI_RANK = RANK_IP*(RANK_I_TOT) + RANK_I 
   RANK_IP         = int(MPI_Rank/RANK_I_TOT);          // const int RANK_IP = 0;
   RANK_I          = MPI_Rank % RANK_I_TOT ;            // const int RANK_I  = MPI_Rank;    
   global_nx_unit  = ceil(NX0_TOT[0]/RANK_I_TOT );
   global_nxp_unit = ceil(NX0_TOT[0]/RANK_IP_TOT);
   
   if (RANK_I  != RANK_I_TOT-1)  global_nx  = global_nx_unit;   
   else                          global_nx  = NX0_TOT[0] - global_nx_unit *(RANK_I_TOT -1);
   if (RANK_IP != RANK_IP_TOT-1) global_nxp = global_nxp_unit;
   else                          global_nxp = NX0_TOT[0] - global_nxp_unit*(RANK_IP_TOT-1);
   
   //Aux_Message(stdout, "(Rank, Rank_I, Rank_Ip) = (%d, %d, %d). \n", MPI_Rank, RANK_I, RANK_IP );
   
   const int  global_nx_start  = RANK_I *global_nx_unit ;
   const int  global_nxp_start = RANK_IP*global_nxp_unit;
   const int  local_ny         = 2*(FFT_Size[1]/2+1);
   const int  local_nz         = FFT_Size[2]; 
   const long slab_size        = local_ny * local_nz ;
   
   // init memory for MPI
   Init_MemAllocate_CylPoisson(slab_size);
   
# endif // ifdef SERIAL, else...
   
   // new MPI Comm
   MPI_Comm_split(MPI_COMM_WORLD, RANK_I,  MPI_Rank, &rank_i_comm );
   MPI_Comm_split(MPI_COMM_WORLD, RANK_IP, MPI_Rank, &rank_ip_comm);
   
   if (rank_i_comm == MPI_COMM_NULL) Aux_Error(ERROR_INFO, "new MPI_Comm initiation failed... \n") ;

   
   // 1.0 build up kernel for FFT
   // 1.1 allocate Kernel K - 
   //   a 1D array (a flattened 2D array), each component is a pointer to another 1D array (a flattened 2D array)
   Aux_AllocateArray2D(KernelFuncK, global_nx*global_nxp, slab_size) ;
   
   // 1.2 build up Kernel in real space
   for (int i=0;  i <global_nx;  i++)  { ii  = i+ global_nx_start;  x  = amr->BoxEdgeL[0] + (ii +0.5)*dh[0];
   for (int ip=0; ip<global_nxp; ip++) { iip = ip+global_nxp_start; xp = amr->BoxEdgeL[0] + (iip+0.5)*dh[0];
                                         ID_planX = i*global_nxp + ip;
   for (int k=0;  k <local_nz;  k++)   { z = ( k <= NX0_TOT[2] ) ? k*dh[2] : (FFT_Size[2]-k)*dh[2] ;
   for (int j=0;  j <local_ny;  j++)   { y = j*dh[1];
      
      ID_planYZ = k*local_ny + j;
      real denominator = SQRT( SQR(x-xp) + (real)2.0*x*xp*((real)1.0 - COS(y)) + SQR(z) ) ;
      
      if (denominator != 0.0) 
         KernelFuncK[ID_planX][ID_planYZ] = (real) -1.0*dh_cube / denominator; 
      else  
         KernelFuncK[ID_planX][ID_planYZ] = (real) 0.0 ;    // mesh does not see itself
         
   }}}}
   
   
   // 2.0 FFT Kernel for each (i, ip)
   for (int i=0;  i<global_nx;   i++){
   for (int ip=0; ip<global_nxp; ip++){
      
      ID_planX = i*global_nxp + ip;
      rfftwnd_one_real_to_complex( FFTW_Plan, KernelFuncK[ID_planX], NULL );      
   }}
      
} // Init_CylKernel


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_CylPoisson
// Description :  allocate memory needed for CylPoisson (promary design for MPI task)
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_CylPoisson(const long slab_size){
   
   if (MPI_Rank == 0) Aux_Message(stdout, "Init_MemAllocate_CylPoisson... ") ;
   
   const int  NSlab            = amr->NPatchComma[0][1]*PS1;   // number of slabs in each node 
   const int  PSSize           = PS1 * PS1;
   const long slab_size_hf     = slab_size / 2 ;
   const long global_nxp_total = global_nxp * NX0_TOT[1] * NX0_TOT[2];
   const long global_nx_total  = global_nx  * NX0_TOT[1] * NX0_TOT[2];
   const long global_nxp_slab  = global_nxp_total / PSSize;
   const long global_nx_slab   = global_nx_total  / PSSize;
   
   
   // 1.0 memory in CPU_CylPoissonSolver
   Aux_AllocateArray2D(RhoK, global_nxp, slab_size) ;
   Aux_AllocateArray2D(PhiK, global_nx , slab_size) ;    //### only RANK_I==0 needs this
   
   // 2.0 memory in Patch2Slab
   SendBuf_Rho      = new real [ NSlab*PSSize ]; 
   SendBuf_IDPlanXp = new int  [ NSlab  ];
   SendBuf_IDPlanYZ = new long [ NSlab  ];
   RecvBuf_Rho      = new real [ global_nxp_total ]; 
   RecvBuf_IDPlanXp = new int  [ global_nxp_slab  ];
   RecvBuf_IDPlanYZ = new long [ global_nxp_slab  ];
   
   // 3.0 memory in Pot_Isolated
   PhiK_All_re      = new real [ global_nx  * slab_size_hf ] ; 
   PhiK_All_im      = new real [ global_nx  * slab_size_hf ] ;
   if (RANK_IP == 0) {
      PhiK_local_re = new real [ global_nx  * slab_size_hf ] ;
      PhiK_local_im = new real [ global_nx  * slab_size_hf ] ;
   }
   
   // 4.0 in Slab2Patch   
   if (RANK_IP == 0) {
      SendBuf_Phi   = new real [ global_nx_total ]; 
      SendBuf_PID   = new long [ global_nx_slab  ];
      SendBuf_I     = new int  [ global_nx_slab  ];
   }
   // can reuse the one in step 2.0
   RecvBuf_Phi      = new real [ NSlab*PSSize ] ;
   RecvBuf_PID      = new long [ NSlab        ] ;
   RecvBuf_I        = new int  [ NSlab        ] ;
   
   //
   if (MPI_Rank == 0) Aux_Message(stdout, "done \n ") ;
   
}  //FUNCTION: Init_MemAllocate_CylPoisson

#endif // GRAVITY
#endif // COORDINATE == CYLINDRICAL