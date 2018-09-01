#include "GAMER.h"

#include <vector>
#include <iterator>

using std::vector;

#if (COORDINATE == CYLINDRICAL)
#ifdef GRAVITY

static void Pot_Isolated( real ** RhoK, real ** PhiK, const long slab_size ) ;


template <class T>
static void Flatten2DVec( vector< vector<T> > Vec, T Array[] ) ;

extern rfftwnd_plan FFTW_Plan, FFTW_Plan_Inv;


//-------------------------------------------------------------------------------------------------------
// Function    :  Patch2Slab
// Description :  (prepare for density)
//-------------------------------------------------------------------------------------------------------
void Patch2Slab(real **RhoK, int SlabID2Rank[], long SlabID2PID[], const double PrepTime, 
                const int local_ny, const int global_nxp_start ) {
                      
   //Aux_Message(stdout, "Rank = %d: In Function <%s>. \n", MPI_Rank, __FUNCTION__);
   
   const int  PSSize           = PS1*PS1;                                // patch slice size
   const int  Scale0           = amr->scale[0];
   const int  NPatchY          = NX0_TOT[1]/PS1;
   const int  NPatchZ          = NX0_TOT[2]/PS1;
   const int  NSlab            = amr->NPatchComma[0][1]*PS1;
   const long NSlabTotal       = NPatchTotal[0]*PS1;
   const int  global_nxp_slab  = global_nxp * NPatchY * NPatchZ;
   const long global_nxp_total = global_nxp * NX0_TOT[1] * NX0_TOT[2];
   
   int  Cr[3], BPos_Xp;
   int  TRANK_IP, TRANK_I, TRank;
   long SlabID;                                 // need this?
   real radius_p; // radius prime (r')
   
   int  SendCount[MPI_NRank], RecvCount[MPI_NRank], SendCount_Rho[MPI_NRank], RecvCount_Rho[MPI_NRank];
   int  SendDisp [MPI_NRank], RecvDisp [MPI_NRank], SendDisp_Rho [MPI_NRank], RecvDisp_Rho [MPI_NRank];
   int  TempBuf_Rank [ NSlab ] ;
   long TempBuf_PID  [ NSlab ], TempBuf_SlabID [ NSlab ] ;
   int  ListAllRank  [ NSlabTotal ];
   long ListAllPID   [ NSlabTotal ], ListAllSlabID  [ NSlabTotal ];
   vector< vector<int>  > TempBuf_IDPlanXp(MPI_NRank);
   vector< vector<long> > TempBuf_IDPlanYZ(MPI_NRank);
   vector< vector<real> > TempBuf_Rho     (MPI_NRank);
   
   
// 2. prepare the temporary send buffer and record lists
   const OptPotBC_t  PotBC_None        = BC_POT_NONE;
   const IntScheme_t IntScheme         = INT_NONE;
   const NSide_t     NSide_None        = NSIDE_00;
   const bool        IntPhase_No       = false;
   const bool        DE_Consistency_No = false;
   const real        MinDens_No        = -1.0;
   const real        MinPres_No        = -1.0;
   const int         GhostSize         = 0;
   const int         NPG               = 1;

   real (*Dens)[PS1][PS1][PS1] = new real [8*NPG][PS1][PS1][PS1];
   
   int idx = 0 ;     // tracking 1D index - TempBuf_Rank, TempBuf_PID, TempBuf_SlabID
   
   for (int r=0; r<MPI_NRank; r++)  SendCount[r] = 0;  // initialization
   
   for (int PID0=0; PID0<amr->NPatchComma[0][1]; PID0+=8) {
//    even with NSIDE_00 and GhostSize=0, we still need OPT__BC_FLU to determine whether periodic BC is adopted
//    also note that we do not check minimum density here since no ghost zones are required
      Prepare_PatchData( 0, PrepTime, Dens[0][0][0], GhostSize, NPG, &PID0, _DENS,
                         IntScheme, UNIT_PATCH, NSide_None, IntPhase_No, OPT__BC_FLU, PotBC_None,
                         MinDens_No, MinPres_No, DE_Consistency_No );
      
      for (int PID=PID0, LocalID=0; PID<PID0+8; PID++, LocalID++) {
         // determine the x/y/z-index in this patch
         for (int d=0; d<3; d++)    Cr[d] = amr->patch[0][0][PID]->corner[d] / Scale0;
         
         for (int ip=0; ip<PS1; ip++) {
            // BPos_Xp = RANK_I*global_nxp + RANK_IP*local_nxp + residual_nxp_in_that_rank
            BPos_Xp  = Cr[0] + ip;
            SlabID   = (long)( BPos_Xp*NPatchZ + int(Cr[2]/PS1) ) * NPatchY + int(Cr[1]/PS1) ;
            radius_p = Aux_Coord_CellIdx2AdoptedCoord(0, PID, 0, ip);
            
            // find TRANK_IP and TRANK_I
            TRANK_IP = int( BPos_Xp/global_nxp_unit );            
            TRANK_I  = 0;  
            TRank    = TRANK_IP*(RANK_I_TOT) + TRANK_I ; 
                        
            TempBuf_IDPlanXp[TRank].push_back( BPos_Xp );
            TempBuf_IDPlanYZ[TRank].push_back( Cr[2]*local_ny + Cr[1] );
            
            for (int k=0; k<PS1; k++) {  
            for (int j=0; j<PS1; j++) {
               TempBuf_Rho[TRank].push_back( Dens[LocalID][k][j][ip] * radius_p );
            }}
            SendCount[TRank] ++;
            
            TempBuf_Rank  [idx] = MPI_Rank ; 
            TempBuf_PID   [idx] = PID ;
            TempBuf_SlabID[idx] = SlabID ;
            idx ++ ;
            
         } // for (ip=0; ip<PS1; ... )
      } // for (PID=PID0, LocalID=0; PID<PID0+8; ... )  
   } // for (PID0=0; PID0<amr->NPatchComma[0][1]; ...)
   
   delete [] Dens;
   
   
// 3.0 prepare SlabID2Rank, SlabID2PID
   int ListAllNSlab[MPI_NRank], NSlabDisp[MPI_NRank]; 
   
   MPI_Allgather( &NSlab, 1, MPI_INT, ListAllNSlab, 1, MPI_INT, MPI_COMM_WORLD ); 
   
   NSlabDisp[0] = 0;
   for (int r=1; r<MPI_NRank; r++) 
      NSlabDisp[r] = NSlabDisp[r-1] + ListAllNSlab[r-1];
   
   MPI_Allgatherv( TempBuf_Rank,   NSlab, MPI_INT,  ListAllRank,   ListAllNSlab, NSlabDisp, MPI_INT,  MPI_COMM_WORLD );
   MPI_Allgatherv( TempBuf_PID,    NSlab, MPI_LONG, ListAllPID,    ListAllNSlab, NSlabDisp, MPI_LONG, MPI_COMM_WORLD );
   MPI_Allgatherv( TempBuf_SlabID, NSlab, MPI_LONG, ListAllSlabID, ListAllNSlab, NSlabDisp, MPI_LONG, MPI_COMM_WORLD );
   
   
   long SID ;
   for (long t=0; t<NSlabTotal; t++) {
      SID              = ListAllSlabID[t]; 
      SlabID2Rank[SID] = ListAllRank[t]; 
      SlabID2PID [SID] = ListAllPID [t];
   }
   
// 3.1 broadcast the number of elements sending to different ranks
   MPI_Alltoall( SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD );

   
   for (int r=0; r<MPI_NRank; r++)  {
      SendCount_Rho[r] = SendCount[r]*PSSize;
      RecvCount_Rho[r] = RecvCount[r]*PSSize;
   }
   
   //// construct SendDisp, RecvDisp
   SendDisp[0]     = 0;
   RecvDisp[0]     = 0;
   SendDisp_Rho[0] = 0;
   RecvDisp_Rho[0] = 0;
   
   for (int r=1; r<MPI_NRank; r++) {
      SendDisp    [r] = SendDisp    [r-1] + SendCount    [r-1] ;
      SendDisp_Rho[r] = SendDisp_Rho[r-1] + SendCount_Rho[r-1] ;
      RecvDisp    [r] = RecvDisp    [r-1] + RecvCount    [r-1] ;
      RecvDisp_Rho[r] = RecvDisp_Rho[r-1] + RecvCount_Rho[r-1] ;
   }
   
   
// 3.2 prepare the send buffer of Rank   
   Flatten2DVec(TempBuf_Rho,      SendBuf_Rho     );
   Flatten2DVec(TempBuf_IDPlanXp, SendBuf_IDPlanXp);
   Flatten2DVec(TempBuf_IDPlanYZ, SendBuf_IDPlanYZ);
   
   
// 4. exchange data by MPI 
// 4.1 first exange data to rank w/ i=0
   MPI_Alltoallv( SendBuf_IDPlanXp, SendCount, SendDisp, MPI_INT,
                  RecvBuf_IDPlanXp, RecvCount, RecvDisp, MPI_INT, MPI_COMM_WORLD );
   
   MPI_Alltoallv( SendBuf_IDPlanYZ, SendCount, SendDisp, MPI_LONG,
                  RecvBuf_IDPlanYZ, RecvCount, RecvDisp, MPI_LONG, MPI_COMM_WORLD );
                  
   MPI_Alltoallv( SendBuf_Rho, SendCount_Rho, SendDisp_Rho, MPI_DOUBLE,
                  RecvBuf_Rho, RecvCount_Rho, RecvDisp_Rho, MPI_DOUBLE, MPI_COMM_WORLD );
                  
// 4.2 broadcast data along ip=const
   MPI_Bcast( RecvBuf_IDPlanXp, global_nxp_slab,  MPI_INT,    0, rank_ip_comm );
   MPI_Bcast( RecvBuf_IDPlanYZ, global_nxp_slab,  MPI_LONG,   0, rank_ip_comm );
   MPI_Bcast( RecvBuf_Rho,      global_nxp_total, MPI_DOUBLE, 0, rank_ip_comm );
   

// 5. store the received density to the padded array "RhoK" for FFTW
   long count = 0 ;
   int ID_planXp, jj ;
   long ID_planYZ;
   
   for (long t=0; t<global_nxp_slab ; t++)
   {
      ID_planXp = RecvBuf_IDPlanXp[t] - global_nxp_start;
      ID_planYZ = RecvBuf_IDPlanYZ[t];

      for (int j=0; j<PS1; j++) {
      for (int i=0; i<PS1; i++) {
         jj = ID_planYZ + j*local_ny + i;
         RhoK[ID_planXp][ jj ] = RecvBuf_Rho[ count ];
         count ++ ;
      }}
   }
   
} // FUNCTION: Patch2Slab


//-------------------------------------------------------------------------------------------------------
// Function    :  Flatten2DArray
// Description :  traverse a 2D vector into a 1D array
//-------------------------------------------------------------------------------------------------------
template <class T>
void Flatten2DVec( vector< vector<T> > Vec, T Array[] ){
   
   typename vector< vector<T> >::iterator row ;
   typename vector<T>::iterator col ;
   
   int curr = 0 ;     // array index
   
   for ( row = Vec.begin();  row != Vec.end();  row++ ) {
   for ( col = row->begin(); col != row->end(); col++ ) {
      Array[curr] = *col ;
      curr ++ ;
   }}
      
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Slab2Patch
// Description :  displace PhiK back to patch data
//                1. fill in TempBufPhi
//                2. flatten TempBufPhi to SendBufPhi for sending
//                3. distribute SendCount
//                4. construct SendDisp, RecvDisp
//-------------------------------------------------------------------------------------------------------
void Slab2Patch(real **PhiK, const int SaveSg, int SlabID2Rank[], long SlabID2PID[],
                const int local_ny, const int global_nx_start, const real Coeff ) {
                   
   //Aux_Message(stdout, "Rank = %d: In Function <%s>. \n", MPI_Rank, __FUNCTION__);
                         
   const int  PSSize        = PS1 * PS1;
   const real fftw_norm     = (real) 1.0 / (real) ( NX0_TOT[1]*((real)2.0*NX0_TOT[2]) ) ;
   const long NRecvSlab     = (long)amr->NPatchComma[0][1]*PS1;   // total number of received patch slices
   const int  NPatchY       = NX0_TOT[1]/PS1;
   const int  NPatchZ       = NX0_TOT[2]/PS1;
   
   vector< vector<real> > TempBuf_Phi( MPI_NRank );
   vector< vector<long> > TempBuf_PID( MPI_NRank );
   vector< vector<int>  > TempBuf_I  ( MPI_NRank );
   
   int  SendCount[MPI_NRank], RecvCount[MPI_NRank], SendCount_Phi[MPI_NRank], RecvCount_Phi[MPI_NRank];
   int  SendDisp [MPI_NRank], RecvDisp [MPI_NRank], SendDisp_Phi [MPI_NRank], RecvDisp_Phi [MPI_NRank];
   int  Cr0, Cr1, Cr2, ii, jj, kk, TRank;
   long SlabID, PID, ID_planYZ;  
   
   // initialization
   for (int r=0; r<MPI_NRank; r++)  SendCount[r] = 0;
   
   // 1. loop over all slabs and save to a 2D TempBuf vector
   if (RANK_IP == 0) {
      
      for (int i=0; i<global_nx; i++ ) {
         Cr0 = global_nx_start + i; 
      
         // loop over all slabs
         for (Cr2 = 0; Cr2 < NX0_TOT[2]; Cr2 += PS1) 
         for (Cr1 = 0; Cr1 < NX0_TOT[1]; Cr1 += PS1) {
         
            SlabID = (long)( Cr0*NPatchZ + int(Cr2/PS1) ) * NPatchY + int(Cr1/PS1) ;
            TRank  = SlabID2Rank[SlabID] ; 
            PID    = SlabID2PID [SlabID] ;
         
            TempBuf_PID[TRank].push_back(PID);
            TempBuf_I  [TRank].push_back(Cr0);
         
            // save Phi in TempBuf_Phi
            for (int k=0; k<PS1; k++) { kk = Cr2 + k;
            for (int j=0; j<PS1; j++) { jj = Cr1 + j;
               ID_planYZ = kk*local_ny + jj ;
               TempBuf_Phi[TRank].push_back( PhiK[i][ID_planYZ] );
            }} // for j, k
         
            SendCount[TRank] ++ ;
         
         } // for Cr2, Cr1
      } // for (int i=0; ... )
   
   
   // 2. flatten 2D TempBuf vector to a 1D SendBuf array 
      Flatten2DVec(TempBuf_Phi, SendBuf_Phi) ;
      Flatten2DVec(TempBuf_PID, SendBuf_PID) ;
      Flatten2DVec(TempBuf_I  , SendBuf_I  ) ; 
   
   } // if RANK_IP == 0
   

   // 3. distribute SendCount in all processors
   MPI_Alltoall( SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD );
   
   for (int r=0; r<MPI_NRank; r++)  {
      SendCount_Phi[r] = SendCount[r]*PSSize;
      RecvCount_Phi[r] = RecvCount[r]*PSSize;
   }

   
   // 4. Construct SendDisp, RecvDisp
   SendDisp[0]     = 0;
   RecvDisp[0]     = 0;
   SendDisp_Phi[0] = 0;
   RecvDisp_Phi[0] = 0;
   
   for (int r=1; r<MPI_NRank; r++) {
      SendDisp    [r] = SendDisp    [r-1] + SendCount    [r-1] ;
      SendDisp_Phi[r] = SendDisp_Phi[r-1] + SendCount_Phi[r-1] ;
      RecvDisp    [r] = RecvDisp    [r-1] + RecvCount    [r-1] ;
      RecvDisp_Phi[r] = RecvDisp_Phi[r-1] + RecvCount_Phi[r-1] ;
   }
   
   
   // 5. distribute SendBuf_Phi, SendBuf_PID, SendBuf_I
   MPI_Alltoallv( SendBuf_Phi, SendCount_Phi, SendDisp_Phi, MPI_DOUBLE, 
                  RecvBuf_Phi, RecvCount_Phi, RecvDisp_Phi, MPI_DOUBLE, MPI_COMM_WORLD );
                  
   MPI_Alltoallv( SendBuf_PID, SendCount, SendDisp, MPI_LONG, 
                  RecvBuf_PID, RecvCount, RecvDisp, MPI_LONG, MPI_COMM_WORLD );
   
   MPI_Alltoallv( SendBuf_I,   SendCount, SendDisp, MPI_INT, 
                  RecvBuf_I,   RecvCount, RecvDisp, MPI_INT, MPI_COMM_WORLD );
   

   // 6. save the patch data back to sandglass 
   int count = 0 ;
   
   for (int t=0; t<NRecvSlab; t++) {
      PID = RecvBuf_PID[t] ; 
      ii  = RecvBuf_I  [t] ;
      int i = ii % PS1 ;
            
      for (int k=0; k<PS1; k++) { 
      for (int j=0; j<PS1; j++) { 
         amr->patch[SaveSg][0][PID]->pot[k][j][i] = RecvBuf_Phi[count] * fftw_norm ; 
         count++ ;
      }}
   }
   
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Pot_Isolated
// Description :  Evaluate the gravitational potential in cyl coordinate
//                1. FFT RhoK for in each rank with size local_nxp 
//                   -> ###
//                   -> ### use reinterpret_cast in FFTW3
//                2. collect RhoK to rach rank to meet size global_nxp
//                   -> ### use rank_ip_comm, instead of MPI_COMM
//                3. integrate to get PhiK with size global_nx in each rank
//                4. distribute PhiK to each rank with local_nx size, and copy data to meet right type
//                   -> ### use reinterpret_cast in FFTW3
//                5. iFFT PhiK back to real space
//                   
//-------------------------------------------------------------------------------------------------------
void Pot_Isolated( real ** RhoK, real ** PhiK, const long slab_size ){
   
   //Aux_Message(stdout, "Rank = %d: In Function <%s>. \n", MPI_Rank, __FUNCTION__);         
   
   fftw_complex *RhoK_cplx, Temp_cplx;
   fftw_complex *PhiK_cplx, *SubKernel;
   real         *PhiK_re_ptr, *PhiK_im_ptr;
   int          ID_planX, CommCount, target_rank;
   const long   slab_size_hf = slab_size/2; 
   
   
   // 1. collect all RhoK along ip=const direction   
   for (int ip=0; ip<global_nxp; ip++) 
      rfftwnd_one_real_to_complex( FFTW_Plan, RhoK[ip], NULL );
   
   // 3. integrate to get PhiK
   // 3.1 set PhiK to zeros
   for (int t=0; t<global_nx*slab_size_hf; t++ )   PhiK_All_re[t] = PhiK_All_im[t] = (real) 0.0;
   
   // 3.2 integrate locally over r' to get partially integrated PhiK in each rank
   for (int i=0; i<global_nx; i++ ){
      PhiK_re_ptr = & PhiK_All_re[i*slab_size_hf] ; 
      PhiK_im_ptr = & PhiK_All_im[i*slab_size_hf] ;
      
      for (int ip=0; ip<global_nxp; ip++){
         ID_planX  = i*global_nxp + ip ;
         
         RhoK_cplx = (fftw_complex *) (RhoK[ip]);
         SubKernel = (fftw_complex *) KernelFuncK[ID_planX];
         
      
         for (long t=0; t<slab_size_hf; t++) {
            Temp_cplx = RhoK_cplx[t];
            
            PhiK_re_ptr[t] += Temp_cplx.re * SubKernel[t].re - Temp_cplx.im * SubKernel[t].im ;
            PhiK_im_ptr[t] += Temp_cplx.re * SubKernel[t].im + Temp_cplx.im * SubKernel[t].re ;            
         }
      } // for (ip=0; ...)
   } // for (i=0; ...)
   
   
   // 4. add PhiK across different rank for total summation/integration
   CommCount = global_nx * slab_size_hf;
   
   MPI_Reduce( PhiK_All_re, PhiK_local_re, CommCount, MPI_DOUBLE, MPI_SUM, 0, rank_i_comm) ; 
   MPI_Reduce( PhiK_All_im, PhiK_local_im, CommCount, MPI_DOUBLE, MPI_SUM, 0, rank_i_comm) ;
   
   
   if (RANK_IP == 0) {
      // copy PhiK_local to PhiK
      // ### can this possibly be done through memcpy
      PhiK_cplx = (fftw_complex *) PhiK[0] ;
      for (int t=0; t<global_nx*slab_size_hf; t++) {
         PhiK_cplx[t].re = PhiK_local_re[t] ;
         PhiK_cplx[t].im = PhiK_local_im[t];  //###
      }
   
   
   // 5. iFFT PhiK back to real space 
      for ( int i=0; i<global_nx; i++ ) {
         PhiK_cplx = (fftw_complex *) PhiK[i] ;
         rfftwnd_one_complex_to_real( FFTW_Plan_Inv, PhiK_cplx, NULL );
      } 
   
   } // if RANK_IP == 0
   
}


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_CylPoissonSolver_FFT
// Description :  
//-------------------------------------------------------------------------------------------------------
void CPU_CylPoissonSolver( const real Poi_Coeff, const int SaveSg, const double PrepTime ){
   
   //if (MPI_Rank == 0) Aux_Message(stdout, "In Function <%s> ... \n", __FUNCTION__);
   
   // determine the FFT size; FFT_Size[0] is redundunt
   int  FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2]*2 };
   int  SlabID2Rank[ NPatchTotal[0]*PS1 ] ;
   long SlabID2PID [ NPatchTotal[0]*PS1 ] ;
   
# ifdef SERIAL
   Aux_Message(stderr, "Cylindrical Self-Gravity is not yet ready for serial mode! \n");
   
# else
   const int  global_nx_start  = RANK_I *global_nx_unit ;
   const int  global_nxp_start = RANK_IP*global_nxp_unit;
   const int  local_ny         = 2*(FFT_Size[1]/2+1);
   const int  local_nz         = FFT_Size[2]; 
   const long slab_size        = local_ny * local_nz ;
   
# endif // ifdef SERIAL, else...
   
   
   // init RhoK array to zerol
   for (int ip=0; ip<global_nxp; ip++)
   for (int t=0;  t<slab_size;  t++ ) {
      RhoK[ip][t] = (real) 0.0;
   }
   
   // ### need to do it to PhiK as well??
   for (int i=0; i<global_nx;  i++)
   for (int t=0; t<slab_size; t++ ) {
      PhiK[i][t] = (real) 0.0;
   }
   
   
   // 1. get (real)RhoK[ip] - Patch2Slab()
   Patch2Slab( RhoK, SlabID2Rank, SlabID2PID, PrepTime, local_ny, global_nxp_start ) ;
   
   // 2.
   if ( OPT__BC_POT == BC_POT_ISOLATED ) 
      Pot_Isolated( RhoK, PhiK, slab_size ) ;
   else
      Aux_Error( ERROR_INFO, "Cylindrical poisson sovler only support isolated boundary condition. \n");
   
   // 3.
   Slab2Patch( PhiK, SaveSg, SlabID2Rank, SlabID2PID, local_ny, global_nx_start, Poi_Coeff ) ;
   

} // CPU_CylPoissonSolver_FFT

#endif // GRAVITY

#endif // COORDINATE == CYLINDRICAL