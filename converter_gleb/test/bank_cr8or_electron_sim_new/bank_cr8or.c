
// Standard C++ Headers:
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>

// 
#include <getopt.h>

// ROOT Headers:
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

extern "C" 
{

#include <signal.h>
#include <errno.h>
#include <ntypes.h>
#include <bostypes.h>
#include <ec.h>
#include <clas_cern.h>
#include <ctype.h>
#include <kinematics.h>
#include <map_manager.h>
#include <trk.h>
#include <clasmdl.h>
#include <utility.h>
#include <pid.h>
#include <makebanks.h>
#include <call.h>
#include <bosddl.h>
#include <tagtnorm.h>
#include <vertex.h>

  void initbos();
  int  getBOS(BOSbank *bcs, int lun, char *list);
  void cleanBanks(BOSbank *bcs);
  void dropAllBanks(BOSbank *bcs, char *list);
  void *getBank(BOSbank *,const char *);
  void open_fpack_unit(char *filename,char *dataname,int unitnum);
  void close_fpack_unit(char *dataname);
  BOSbank bcs_ ;
  BOSbank wcs_ ;

  float ranf_();
  int ranset_( float* );
  void ranlux_( float* , const int* ) ;
 


//  time_t   time( time_t );
//  int   clock();
}





using namespace std;

//clasPART_t*  CreatePART();
//clasHEAD_t*  CreateHEAD( unsigned );

int main()
{

clasMCTK_t*  CreateMCTK();
clasMCTK_t*  CreateMCVX();
clasMCTK_t*  CreateMCEV();
clasHEAD_t*  CreateHEAD( unsigned );
  unsigned EvtNumb = 10000; 

  char* FileName = "test.bos" ;
  char mess[256];

  unsigned EvtReport = 1000000;

  // Get SEED from the time of the system 
  //  and use it to set CERNLIB seed 
  float seed = (float) time( NULL );
  ranset_( &seed );

//  cout << " Clock gives " << seed << " value " << endl ;


  // Initialize bos and open the file for writing
  //
  initbos();
  //  read_ddl();
  bankList( &bcs_, "E=", "HEADPARTMCTKMCVXMCEV" );
  sprintf( mess, "OPEN BOSOUTPUT UNIT=1 FILE=\"%s\" WRITE STATUS=NEW RECL=3600", FileName );
  fparm_c( mess );


  // Loop over all events 
  //
  for ( unsigned iEvt = 0; iEvt < EvtNumb; iEvt++ )
    {
      if ( iEvt%EvtReport == 0 )
     cout << "Generating Event # " << ( iEvt + 1 ) << endl ; 

      CreateHEAD( iEvt);
      CreateMCEV();
      CreateMCTK();
      CreateMCVX(); 

      // Write banks in the list to the file and 
      // then tidy up 
      putBOS( &bcs_, 1, "E" );
      dropAllBanks( &bcs_, "E");
      cleanBanks( &bcs_ );
    }
  // Close the outupt file 
  close_fpack_unit("BOSOUTPUT");   
}


clasHEAD_t*  CreateHEAD( unsigned iEvt )
{
  clasHEAD_t* HEAD ;
  if ( ( HEAD = (clasHEAD_t*)makeBank( &bcs_, "HEAD", 0, sizeof(head_t)/4, 1 ) ) )
    {
      HEAD->head[0].version = 1; 
      HEAD->head[0].nevent = iEvt+1 ;
      HEAD->head[0].nrun = 10;
      HEAD->head[0].evtclass = 7 ;
      HEAD->head[0].type = -2 ;
      HEAD->head[0].time = 0;
      HEAD->head[0].trigbits = 1 ;
      
      return (clasHEAD_t*) HEAD ;
    }
  else
    {
      return (clasHEAD_t*) NULL;
    }   
}

clasMCEV_t* CreateMCEV()
{
  const int len = 3;
  float rndm_vec[len] ;
  clasMCEV_t* MCEV;
  MCEV = (clasMCEV_t*) makeBank( &bcs_, "MCEV", 0, sizeof(mcev_t)/4, 1 ) ;
 
   ranlux_( &rndm_vec[0], &len );
   
      MCEV->mcev[0].i1 = 100000*rndm_vec[0]; 
      MCEV->mcev[0].i2 =100000*rndm_vec[1];  
      
      return (clasMCEV_t*) MCEV;

};


clasMCVX_t* CreateMCVX()
{
  const int len = 3;
  float rndm_vec[len] ;
  clasMCVX_t* MCVX;
  MCVX = (clasMCVX_t*) makeBank( &bcs_, "MCVX", 0, sizeof(mcvx_t)/4, 1 ) ;
 
   ranlux_( &rndm_vec[0], &len );
   
      MCVX->mcvx[0].x = 0.; 
      MCVX->mcvx[0].y = 0.;
      MCVX->mcvx[0].z = -0.4  - 2./2. + 2.*rndm_vec[0];
      MCVX->mcvx[0].tof = 0.;
      MCVX->mcvx[0].flag = 0;
      
   
      
      return (clasMCVX_t*) MCVX;

};
  
clasMCTK_t* CreateMCTK()
{
  const int len = 3;
  float rndm_vec[len] ;
  clasMCTK_t* MCTK;
  MCTK = (clasMCTK_t*) makeBank( &bcs_, "MCTK", 0, sizeof(mctk_t)/4, 1 ) ;



  // Fill proton MCTK bank 

  ranlux_( &rndm_vec[0], &len );

/*  float p_p = 0.2 + 1.2 * rndm_vec[0] ;
//  float p_p = 0.5 ;
//  float Th_p  = ( 5.   + 25. * rndm_vec[1] ) * acos(0.) / 90. ;
  float Th_p  = ( 10.   + 120. * rndm_vec[1] ) * acos(0.) / 90. ;
  float Phi_p = ( -30. + 360. * rndm_vec[2] ) * acos(0.) / 90. ;
//  float Phi_p = ( 180.) * acos(0.) / 90. ;

  PART->part[0].pid = 14 ;
  PART->part[0].p.t  = sqrt( p_p*p_p + 0.938*0.938 ) ;
  PART->part[0].p.space.x = p_p * sin( Th_p ) * cos( Phi_p );
  PART->part[0].p.space.y = p_p * sin( Th_p ) * sin( Phi_p );
  PART->part[0].p.space.z = p_p * cos( Th_p ) ;
  PART->part[0].vert.x = 0;
  PART->part[0].vert.y = 0;
  PART->part[0].vert.z = 0;
  PART->part[0].q = 1 ;
 */
 
float p_p = 0.2 + 1.8 * rndm_vec[0] ;
//  float p_p = 0.5 ;
//  float Th_p  = ( 5.   + 25. * rndm_vec[1] ) * acos(0.) / 90. ;
  float Th_p  = ( 10.   + 40. * rndm_vec[1] ) * acos(0.) / 90. ;
  float Phi_p = ( -30. + 360. * rndm_vec[2] ) * acos(0.) / 90. ;
//  float Phi_p = ( 180.) * acos(0.) / 90. ;
     
      MCTK->mctk[0].cx =  sin( Th_p ) * cos( Phi_p );
      MCTK->mctk[0].cy =  sin( Th_p ) * sin( Phi_p );
      MCTK->mctk[0].cz = cos( Th_p ) ;
      MCTK->mctk[0].pmom = p_p;
      MCTK->mctk[0].mass = 0.001;
      MCTK->mctk[0].charge = -1.;
      MCTK->mctk[0].id = 11;
      MCTK->mctk[0].flag = 0;
      MCTK->mctk[0].beg_vtx = 1;
      MCTK->mctk[0].end_vtx = 0;
      MCTK->mctk[0].parent = 0;
 
  return (clasMCTK_t*) MCTK;

}


