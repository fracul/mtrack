#include <math.h>
#include "types.h"
#include "def.h"
#include "rw.h"

/* Global variables */
extern ring_t ring;
extern tracking_t track;

int transformSbunch_shortrange(in, kb, ebeam, MPmodel, bunch0, bunch)
  int in;
  unsigned int kb;
  const e_beam_t ebeam;
  const bunch_macroparticle_model_t MPmodel;
  bunch_strong_t * bunch0;
  bunch_strong_t * bunch;
{ 
  int    icell;/* jp;*/
  double xeps_gainj;

/*------------------------------------------------------------------------------
 * The Principal Transformations in Each Plane:
 *----------------------------------------------------------------------------*/
 
  for(icell = 0; icell < 2 * ebeam.distrib.Nslc + 1; icell++)
  {
    xeps_gainj = ring.Vrf0/(ring.E0*FKILO)*sin(ring.wrf*bunch0->xtau[icell]+ring.phai0);
    
    bunch->xeps[icell]  = bunch0->xeps[icell]  + xeps_gainj - ring.Urad;
    bunch->xtau[icell]  = bunch0->xtau[icell]  - bunch0->xeps[icell]*ring.T0*ring.ac;
  }
 /* weakweak equivalent:
 
  for(jp = 0; jp < MPmodel.Np; jp++) { 
     xeps_gainj = ring.Vrf0/(ring.E0*FKILO)*sin(ring.wrf*xtau0[jp]+ring.phai0);    

     xepsN[jp]  = xeps0[jp]  + xeps_gainj - ring.Urad;
     xtauN[jp]  = xtau0[jp]  - xepsN[jp]*ring.T0*ring.ac;
  } 
  */
/*----------------------------------------------------------------------------*/
  /*** if(0 > intrabunch_vertical_wakepotentials(kb));  ***/ 
	
  bunch_strong_updateCM(bunch, 2 * ebeam.distrib.Nslc + 1, LON);
  bunch_strong_updateCM(bunch, 2 * ebeam.distrib.Nslc + 1, HOR);
  bunch_strong_updateCM(bunch, 2 * ebeam.distrib.Nslc + 1, VER);

/*----------------------------------------------------------------------------*/


/* No shaking in weakstrong model
    if(in >= track.NrevTot-NrevShk  && logShaker)  if(0 > apply_shaker(in)) ;
 */


/*------------------------------------------------------------------------------
 * Updating the variables: 
 *----------------------------------------------------------------------------*/
  for(icell = 0; icell < 2 * ebeam.distrib.Nslc + 1; icell++)
  {
    bunch0->xtau[icell] = bunch->xtau[icell];
    bunch0->xeps[icell] = bunch->xeps[icell];
  }   
/*-----------------------------------------------------------------------------*/
  return 1;
}

int
transformSbunch_longrange(int in, unsigned int kb, const plane_t plane,
                         const e_beam_t ebeam,
                         const bunch_macroparticle_model_t MPmodel,
                         const bunch_CM_history_strong_t * bCMhist,
                         bunch_strong_t * bCM)
{
/*------------------------------------------------------------------------------
 * Calculate the influence on kb-th bunch of long-range RW wakefields produced 
 * by different bunches.
 * We assume that the field variation over the bunch length can be neglected
 * both for the source bunch and kb-th bunch. The dipole moment of a bunch is
 * thus obtained by simply via the centre of mass position ignoring the
 * longitudinal difference.
 *
 * Multi-turn effect is considered up to NmultiT turns.
 *----------------------------------------------------------------------------*/
  int icell;
  int jmlt;
  
  const int NDlong = 1;   /* Number of division around the ring for the long range transformation */
  const double fNDlong = ((double) NDlong);
  double Psi0;
  double Psi = 0.0;
  double beff3;
  switch(plane)
  {
    case HOR:
      Psi0 = 2.0*M_PI*ring.QH0/fNDlong;
      beff3       = pow(ring.beffH[0], 3);
    break;
    case VER:
      Psi0 = 2.0*M_PI*ring.QV0/fNDlong;
      beff3       = pow(ring.beffV[0], 3);
    break;
    default:
      return -1;
    break;
  }
  
  const int Ntrigger = 50;
  if(track.triggerRW && in < Ntrigger) {
     beff3 =  pow(2.0e-3, 3);
  }
  
  double uuhistjb[30];
  /*double uphistjb[30]; TODO check*/
  double Wrwt, Dipolejb;
  
  const double Pref0 = ring.T0/ring.E0/MPmodel.fNp/FTERA;
  const double rhorw  = 2.8e-08;
  const double kappa  = 1.0/rhorw;
  const double fNharm = ((double) ring.Nharm);
  const double Dub0 = ring.Lc/fNharm;
  double Dup;
  
  const double prfrw = 2.0*ring.R/(beff3)/sqrt(M_PI*kappa*EPSILON0)/fNDlong;
  
  /* Matrix transforming CMs of all other jb bunches to the is-th observation point: */
  Psi  = Psi + Psi0; 
  double cosPsi  = cos(Psi),   sinPsi = sin(Psi);
  double amv11 = cosPsi + ring.alpha1[plane]*sinPsi;
  double amv12 =        + ring.beta1[plane]*sinPsi;
  
  /*double amv21 =        - ring.gamma1[plane]*sinPsi;
  double amv22 = cosPsi - ring.alpha1[plane]*sinPsi; TODO check */
  
  unsigned int jb;
  for(jb = 0; jb < ring.Nharm; jb++)
  {
    if(ebeam.nfFill[jb])
    {
      int j;
      double * up = (plane == HOR) ? bCM->xp : bCM->zp;
      const double * uuhist = (plane == HOR) ? bCMhist->xxhist[jb] : bCMhist->zzhist[jb];
      const double * uphist = (plane == HOR) ? bCMhist->xphist[jb] : bCMhist->zphist[jb];
    
      /*** Transform CM of jb-th bunch to is-th position: (is=0,...,NDlong-1) ***/
      for(j = 0; j < track.Nmlt + 2; j++)
      {
        uuhistjb[j] = amv11*uuhist[j] + amv12*uphist[j];
        /*uphistjb[j] = amv21*uuhist[j] + amv22*uphist[j]; TODO check*/
      }
      
      Dup = 0.0;
/*----------------------------------------------------------------------------*/
/* Getting the time distance between the kb-th bunch and jb-th bunch:         */  
/*----------------------------------------------------------------------------*/
      double Dujk0  = ((double) (jb - kb))*ring.Lc/((double) ring.Nharm);   /***  NB: jb=0 corresponds to the tail ***/
      if(jb == kb)  Dujk0 = ring.Lc;
      if(Dujk0 < 0) Dujk0 = ring.Lc + Dujk0;
      
      if(jb > kb) Dujk0 = ((double)(jb - kb))*Dub0;
      if(jb < kb) Dujk0 = ((double)(ring.Nharm + jb - kb))*Dub0;
      
      for(j = 0; j < track.Nmlt + 1; j++)
      {
        double Dujk = Dujk0 + ((double) j)*ring.Lc;
        double Dtjk = Dujk/C_LIGHT;
	Wrwt = prfrw/sqrt(Dtjk);
/*----------------------------------------------------------------------------*/
/* Getting the dipole moment of jb-th bunch:                                  */  
/*----------------------------------------------------------------------------*/
        if(jb > kb) jmlt = j; 
	else        jmlt = j+1;
        Dipolejb = MPmodel.fNp*uuhistjb[jmlt];
	Dup      = Dup + Dipolejb*Wrwt;
      }
      Dup = Pref0*ring.Ibunch[jb]*Dup;
        
      if(in < track.Nmlt) Dup = 0.0; /* 10 replaced by track.Nmlt */
/*------------------------------------------------------------------------------
 * Applying the obtained wakeforce of the jb-th bunch on the particles
 * of the kb-th bunch:
 *----------------------------------------------------------------------------*/
      for(icell = 0; icell < 2*ebeam.distrib.Nslc+1; icell++)
      {
        up[icell] += Dup;
      }
    }
      
  } /*** End of loop over jb ***/
  
  return 1;
}

int transformSbunch_optic(in, kb, ebeam, bCM)
  int    in, kb;
  const e_beam_t ebeam;
  bunch_strong_t * bCM;
{ 
  int    icell, is;
  
  double yyCMN, ypCMN;
  double xxCMN, xpCMN;
  
  int    NDlong = 1;   /* Number of division around the ring for the long range transformation */
  double fNDlong = ((double) NDlong); 
  int    NmultiTsave = track.NmultiT;

/*----------------------------------------------------------------------------*/
  
  if(in < track.NmultiT) track.NmultiT = in+1; 
  else                   track.NmultiT = NmultiTsave;
  
  /*** if(in == 0) NmultiT = 1;  ***/
  /*** if(in == 1 && NmultiT > 2) NmultiT = 2; ***/
  
  double PsiV0   = 2.0*M_PI*ring.QV0/fNDlong; 
  double PsiH0   = 2.0*M_PI*ring.QH0/fNDlong; 
/*----------------------------------------------------------------------------*/
    
  for(is=0; is<NDlong; is++)
  {
    /*** Transform the CM of kb-th bunch to the is-th observation point: ***/
    double PsiVj  = PsiV0; 
    double cosVj  = cos(PsiVj),   sinVj = sin(PsiVj);
    double amv11j = cosVj + ring.alpha1[VER]*sinVj,   amv12j =         ring.beta1[VER]*sinVj;
    double amv21j =       - ring.gamma1[VER]*sinVj,   amv22j = cosVj - ring.alpha1[VER]*sinVj;
       
    double PsiHj  = PsiH0; 
    double cosHj  = cos(PsiHj),   sinHj = sin(PsiHj);
    double amh11j = cosHj + ring.alpha1[HOR]*sinHj,   amh12j =         ring.beta1[HOR]*sinHj;
    double amh21j =       - ring.gamma1[HOR]*sinHj,   amh22j = cosHj - ring.alpha1[HOR]*sinHj;
    
    for(icell=0; icell<2*ebeam.distrib.Nslc+1; icell++)
    {
        yyCMN = amv11j*bCM->zz[icell] + amv12j*bCM->zp[icell];
        ypCMN = amv21j*bCM->zz[icell] + amv22j*bCM->zp[icell];
        bCM->zz[icell] = yyCMN; 
        bCM->zp[icell] = ypCMN; 

        xxCMN = amh11j*bCM->xx[icell] + amh12j*bCM->xp[icell];
        xpCMN = amh21j*bCM->xx[icell] + amh22j*bCM->xp[icell];
        bCM->xx[icell] = xxCMN; 
        bCM->xp[icell] = xpCMN; 
        
    }/*** End of loop over icell ***/
         
  } /*** End of loop over is ***/

/*----------------------------------------------------------------------------*/
  track.NmultiT = NmultiTsave;
/*----------------------------------------------------------------------------*/
  return 1;
}
