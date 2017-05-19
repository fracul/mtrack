#include <stdlib.h>
#include <math.h>
#include "def.h"
#include "types.h"
#include "bunch.h"
#include "transform_weak.h"
#include "cyclic_array.h"
#include "confmpi.h"
#include "feedback_rf.h"

/* Global variables */
//extern ring_t ring;
extern tracking_t track;


//*****************************************************//
//*************** Optics section **********************//
//*****************************************************//

int
transform_weak_bunch_optic(weak_bunch_t * bunch, const bunch_macroparticle_model_t * bunchModel, int iseed, int kb, const ring_t * ring)
{
  unsigned int jp;
  double xeps_gainj;  
  int j;
  
  /* Number of division around the ring for the long range transformation.
   If that ever will be changed, the idea is to alter transverse optics trafo
   with longrange RW. Think of the leap-frog scheme! */
//   int    NDlong = 1;
//   double fNDlong = ((double) NDlong); 
  
  double PsiH0, PsiV0;

  /* Longitudinal trafo, ALWAYS */
  const double damp = 1.0 - ring->De;
  const double excite = 2.0 * bunchModel->slope_sgm.v[LON] * sqrt(ring->T0/ring->taue);
  const double quantumV_pos = bunchModel->pos_sgm.v[VER] * sqrt(2.0 * ring->T0 / ring->tauz);
  const double quantumV_slope = bunchModel->slope_sgm.v[VER] * sqrt(2.0 * ring->T0 / ring->tauz);
  const double quantumH_pos = bunchModel->pos_sgm.v[HOR] * sqrt(2.0 * ring->T0 / ring->taux);
  const double quantumH_slope = bunchModel->slope_sgm.v[HOR] * sqrt(2.0 * ring->T0 / ring->taux);
  
  particle_t particle0[bunch->Np];
  particle_t * particle;
  double rffocus[bunch->Np];
  double eps_const = ring->Vrf0 / (ring->E0 * FKILO);
  active_HC_t * aHC;
  
  if(track.EnableIdealHC == 0)
  {
  /* Storing initial position and longitudinal trafo */  
    for(jp = 0; jp < bunch->Np; jp++)
    {
      particle0[jp] = bunch->particles[jp];
      particle = &(bunch->particles[jp]);
      /* xeps_gainj: energy gain by the RF cavity */
      xeps_gainj = eps_const * sin(ring->wrf * particle0[jp].pos.xtau + ring->phai0);   
      
      // adding potential of active HCs
      for(j = 0; j < ring->active_HC_size; j++)
      {
        aHC = &(ring->active_HC[j]);

        xeps_gainj += aHC->Vpeak/(ring->E0 * FGIGA) * sin(ring->wrf * aHC->nHC * particle0[jp].pos.xtau + aHC->phi_aHC + aHC->nHC * kb * 2 * M_PI);

      }      
      
      /* Urad: energy loss by radiation */
      particle->slope.xtau = particle0[jp].slope.xtau  + xeps_gainj - ring->Urad;
      
      if(track.EnableQuantum)
      {     
        particle->slope.xtau = particle->slope.xtau * damp + excite * box_muller(iseed);
        rffocus[jp] = (1.0 + particle0[jp].slope.xtau)/(1.0 + particle0[jp].slope.xtau + xeps_gainj);
      }
      particle->pos.xtau = particle0[jp].pos.xtau - particle->slope.xtau * ring->T0*ring->ac;
    }
  }
  else
  {  
    /* Including an ideal HC */
    for(jp = 0; jp < bunch->Np; jp++)
    {
      particle0[jp] = bunch->particles[jp];
      particle = &(bunch->particles[jp]);
      /* xeps_gainj: energy gain by the RF cavity & active HC */     
      xeps_gainj = eps_const * (sin(ring->wrf * particle0[jp].pos.xtau + ring->phai0) + ring->HC_k * sin(ring->m_aHC * ring->wrf * particle0[jp].pos.xtau + ring->m_aHC * ring->phi_n));    

      particle->slope.xtau = particle0[jp].slope.xtau  + xeps_gainj - ring->Urad;
      if(track.EnableQuantum)
      {     
        particle->slope.xtau = particle->slope.xtau * damp + excite * box_muller(iseed);
        rffocus[jp] = (1.0 + particle0[jp].slope.xtau)/(1.0 + particle0[jp].slope.xtau + xeps_gainj);
      }
      particle->pos.xtau = particle0[jp].pos.xtau - particle->slope.xtau * ring->T0*ring->ac;
    }
  }

      
    //**********************************************************************************//  
    if(track.TrackPlane[VER])     
    {
      double PsiVj, cosVj, sinVj, amv11j, amv21j, amv12j, amv22j;
      double dQVj = 0.0;
      double EVj, EHj;

      PsiV0   = 2.0*M_PI*ring->QV0;
//       PsiV0   = 2.0*M_PI*ring->QV0/fNDlong;
//       for(is=0; is<NDlong; is++)
//       {
        for(jp = 0; jp < bunch->Np; jp++)
        {
          particle = &(bunch->particles[jp]);
      
	  if (ring->AmpDependdQV) {
	    EVj = ring->gamma1[VER]*particle0[jp].pos.z*particle0[jp].pos.z
	      +2*ring->alpha1[VER]*particle0[jp].pos.z*particle0[jp].slope.z
	      +ring->beta1[VER]*particle0[jp].slope.z*particle0[jp].slope.z;
	    EHj = ring->gamma1[VER]*particle0[jp].pos.x*particle0[jp].pos.x
	      +2*ring->alpha1[VER]*particle0[jp].pos.x*particle0[jp].slope.x
	      +ring->beta1[VER]*particle0[jp].slope.x*particle0[jp].slope.x;
	    dQVj = ring->CVH*EHj+ring->CVV*EVj;
	  }

          PsiVj  = PsiV0 * (1.0 + ring->Gziz * particle0[jp].slope.xtau+dQVj);
          cosVj  = cos(PsiVj);
          sinVj  = sin(PsiVj);
          amv11j = cosVj + ring->alpha1[VER]*sinVj,   amv12j =         ring->beta1[VER]*sinVj;
          amv21j =       - ring->gamma1[VER]*sinVj,   amv22j = cosVj - ring->alpha1[VER]*sinVj;
      
          particle->pos.z = amv11j*particle0[jp].pos.z + amv12j*particle0[jp].slope.z;
          particle->slope.z = amv21j*particle0[jp].pos.z + amv22j*particle0[jp].slope.z;
          
          if(track.EnableQuantum)
          {
            particle->pos.z = particle->pos.z + quantumV_pos * box_muller(iseed);
            particle->slope.z = rffocus[jp] * particle->slope.z + quantumV_slope * box_muller(iseed); 
          }
        }/*** End of loop over jp ***/
//       }
    }
      
    //**********************************************************************************//  
    if(track.TrackPlane[HOR])     
    { 
      double PsiHj, cosHj, sinHj, amh11j, amh12j, amh13j, amh21j, amh22j, amh23j;
      double dQHj = 0;
      double EVj, EHj;

      PsiH0   = 2.0*M_PI*ring->QH0; 
//       PsiH0   = 2.0*M_PI*ring->QH0/fNDlong; 
//       for(is=0; is<NDlong; is++)
//       {
        /*** Transform the CM of kb-th bunch to the is-th observation point: ***/
        for(jp = 0; jp < bunch->Np; jp++)
        {
          particle_t * particle = &(bunch->particles[jp]);

	  if (ring->AmpDependdQV) {
	    EVj = ring->gamma1[VER]*particle0[jp].pos.z*particle0[jp].pos.z
	      +2*ring->alpha1[VER]*particle0[jp].pos.z*particle0[jp].slope.z
	      +ring->beta1[VER]*particle0[jp].slope.z*particle0[jp].slope.z;
	    EHj = ring->gamma1[VER]*particle0[jp].pos.x*particle0[jp].pos.x
	      +2*ring->alpha1[VER]*particle0[jp].pos.x*particle0[jp].slope.x
	      +ring->beta1[VER]*particle0[jp].slope.x*particle0[jp].slope.x;
	    dQHj = ring->CHH*EHj+ring->CHV*EVj;
	  }
       
          PsiHj  = PsiH0 * (1.0 + ring->Gzix * particle0[jp].slope.xtau + dQHj);
          cosHj  = cos(PsiHj);  
          sinHj = sin(PsiHj);          
          amh11j = cosHj + ring->alpha1[HOR]*sinHj;
          amh12j =         ring->beta1[HOR]*sinHj;
          amh13j = (1.0 - amh11j)*ring->dispH1 - amh12j*ring->disppH1;
          amh21j =       - ring->gamma1[HOR]*sinHj;
          amh22j = cosHj - ring->alpha1[HOR]*sinHj;
          amh23j =       -amh21j *ring->dispH1 + (1.0 - amh22j)*ring->disppH1;
          
          particle->pos.x = amh11j*particle0[jp].pos.x + amh12j*particle0[jp].slope.x  + amh13j*particle->slope.xtau;
          particle->slope.x = amh21j*particle0[jp].pos.x + amh22j*particle0[jp].slope.x + amh23j*particle->slope.xtau;
          
          if(track.EnableQuantum)
          {
            particle->pos.x = particle->pos.x + quantumH_pos * box_muller(iseed);
            particle->slope.x = rffocus[jp] * particle->slope.x + quantumH_slope * box_muller(iseed); 
          } 
        }/*** End of loop over jp ***/
         
//       } 
  }
  return 1;
}



//***************************************************************//
//*************** Selffield resonator section *******************//
//**************************************************************//

bool
selffield_model_init(selffield_model_t * SelfFieldModel, double sgm_xtau, weak_bunch_t * bunch)
{
    SelfFieldModel->PlaneL = 0;
    SelfFieldModel->PlaneH = 0;
    SelfFieldModel->PlaneV = 0;
    bunch->N_trash_low = 0;
    bunch->N_trash_high = 0;

  if(SelfFieldModel->Ncell > 0)
  {

    SelfFieldModel->Gl =
      (double *) calloc(SelfFieldModel->Ncell, sizeof(double));
    if (SelfFieldModel->Gl == NULL)
      return false;

    SelfFieldModel->Gh =
      (double *) calloc(SelfFieldModel->Ncell, sizeof(double));
    if (SelfFieldModel->Gh == NULL)
      return false;

    SelfFieldModel->Gv =
      (double *) calloc(SelfFieldModel->Ncell, sizeof(double));
    if (SelfFieldModel->Gv == NULL)
      return false;
    
    SelfFieldModel->dT = 2.0 * SelfFieldModel->Nsigma / ((double) SelfFieldModel->Ncell - 1);
    SelfFieldModel->sigma_tau = sgm_xtau;
    
  }
  else
  {
    SelfFieldModel->Gl = NULL;
    SelfFieldModel->Gh = NULL;
    SelfFieldModel->Gv = NULL;
  }

  return true;
}


void
selffield_model_destroy(selffield_model_t * SelfFieldModel)
{
  SelfFieldModel->Ncell = 0;

  if(SelfFieldModel->Gl != NULL)
    free(SelfFieldModel->Gl);
  SelfFieldModel->Gl = NULL;

  if(SelfFieldModel->Gh != NULL)
    free(SelfFieldModel->Gh);
  SelfFieldModel->Gh = NULL;

  if(SelfFieldModel->Gv != NULL)
    free(SelfFieldModel->Gv);
  SelfFieldModel->Gv = NULL;
  
  SelfFieldModel->dT = 0.0;
  SelfFieldModel->sigma_tau = 0.0;
  
  SelfFieldModel->PlaneL = 0;
  SelfFieldModel->PlaneH = 0;
  SelfFieldModel->PlaneV = 0; 
}




/* Calculates the greensfuctions from all resonator impedancens */
/*---------------------------------------------------------------------------------------------*/
/*  The way how the Greens function and the wake potential are computed follows completely     */
/*  that developed in the longitudinal tracking code of G. Besnier.      (April 2005)          */
/*                                                                                             */
/*  Modified the routine to handle different impedance components:       (September 2005)      */
/*       - Purely resistive                                                                    */
/*       - Purely inductive                                                                    */
/*       - BBR                                 3    69     69     1.000000e+00     -1.109145e+00                                                */
/*  The Green's function constructed here is the sum of all components, and the way it is      */
/*  defined for different types of impedance follows that of G. Besnier                        */
/*  (its correctness confirmed).                                                               */
/*  Impedance components treated here follow "sbtrack_Zlong.inp"                               */
/*                                                                                             */
/* <NB> To handle the sum of many different impedance components, the definition of Gl[]       */
/*      has been changed (The original Gl[] is multiplied by RL/(sgm_xtau/fnan), to make the   */
/*      multiplication factor factG1 independent of the impedance (19/09/2005)                 */
/*---------------------------------------------------------------------------------------------*/
int
construct_greensfunc_resonator(selffield_model_t * SelfFieldModel, 
                               const ring_t * ring)
{
  const double sgm_xtau = SelfFieldModel->sigma_tau; /* [s] */
  unsigned j;

  for (j = 0; j < ring->resonators_size; j++)
  {
    const resonator_t ResImp = ring->resonators[j];
    if (track.TrackPlane[ResImp.plane])
    { 
      double R = ResImp.Rs;
      double Q = ResImp.Qfactor;
      double w = ResImp.wr;                /* [rad/s] */
      const plane_t Pl = ResImp.plane;

      const double dTau = SelfFieldModel->dT;
      double resol = w * sgm_xtau * dTau;
  
      int nc;
      double fnc, arg, argE;
      switch (Pl)
      {
        case LON:
        {
          double sq = sqrt(1.0 - 1.0 / (4.0 * Q * Q));  
          double ar = 0.5 * resol / Q;
          double arq = sq * resol;
          double fac = 0.5 / (Q * sq);
          double ampl_long = w * R / Q;
          SelfFieldModel->Gl[0] += 0.5 * ampl_long;
          for (nc = 1; nc < SelfFieldModel->Ncell; nc++)
          {
            fnc = ((double) nc);
            arg = fnc * arq;
            argE = fnc * ar;            
            SelfFieldModel->Gl[nc] +=
              ampl_long * exp(-argE) * (cos(arg) - fac * sin(arg));  
          }
          SelfFieldModel->PlaneL ++;
          break;
        } /* end LON */

        case VER:
        {
          if(Q > 0.5)
          {
            double sq = sqrt(4.0 * Q * Q - 1.0);  
            double ar = 0.5 * resol / Q;
            double arq = sq * resol / (2 * Q);
            double ampl_trans = 2 * R * w / sq;
            for (nc = 0; nc < SelfFieldModel->Ncell; nc++)
            {
              fnc = ((double) nc);
              arg = fnc * arq;
              argE = fnc * ar;
              SelfFieldModel->Gv[nc] += ampl_trans * exp(-argE) * sin(arg);
            }
          }
          if(Q < 0.5)
          {
            double sq = sqrt(1.0 - 4.0 * Q * Q);  
            double ar = 0.5 * resol / Q;
            double arq = sq * resol / (2 * Q);
            double ampl_trans = 2 * R * w / sq;
            for (nc = 0; nc < SelfFieldModel->Ncell; nc++)
            {
              fnc = ((double) nc);
              arg = fnc * arq;
              argE = fnc * ar;
              SelfFieldModel->Gv[nc] += ampl_trans * exp(-argE) * sinh(arg);
            }
          }
          SelfFieldModel->PlaneV ++;
          break;
        } /* end VER */

        case HOR:
        {
          if(Q > 0.5)
          {
            double sq = sqrt(4.0 * Q * Q - 1.0);  
            double ar = 0.5 * resol / Q;
            double arq = sq * resol / (2 * Q);
            double ampl_trans = 2 * R * w / sq;
            for (nc = 0; nc < SelfFieldModel->Ncell; nc++)
            {
              fnc = ((double) nc);
              arg = fnc * arq;
              argE = fnc * ar;
              SelfFieldModel->Gh[nc] += ampl_trans * exp(-argE) * sin(arg);
            }
          }
          if(Q < 0.5)
          {
            double sq = sqrt(1.0 - 4.0 * Q * Q);  
            double ar = 0.5 * resol / Q;
            double arq = sq * resol / (2 * Q);
            double ampl_trans = 2 * R * w / sq;
            for (nc = 0; nc < SelfFieldModel->Ncell; nc++)
            {
              fnc = ((double) nc);
              arg = fnc * arq;
              argE = fnc * ar;
              SelfFieldModel->Gh[nc] += ampl_trans * exp(-argE) * sinh(arg);
            }
          }
          SelfFieldModel->PlaneH ++;
          break;
        } /* end HOR */
      } /* end switch */
    } /* end if */
  } /* end for */
  return 1;
} /* end fuction */



int
construct_greensfunc_RW(selffield_model_t * SelfFieldModel, 
                               const ring_t * ring)
{
  int nc;

  int ncellmax   = SelfFieldModel->Ncell;
  double sgmatau = SelfFieldModel->sigma_tau;
  double dTau    = SelfFieldModel->dT;
  double beffL1       = ring->beffL[0];      /*** Representative value of the chamber radius and its potentials*/
  double beffL2       = pow(ring->beffL[1], 2);
  double beffV3       = pow(ring->beffV[0], 3);      
  double beffH3       = pow(ring->beffH[0], 3);      
  double beffV4       = pow(ring->beffV[1], 4);      
  double beffH4       = pow(ring->beffH[1], 4);
  double sigmarw = 1.0/ring->rhorw;  /*** Representative value of the chamber conductivity [1/(ohm*m)]   ***/
  double s0;
  double pi = M_PI;
  double t_const = sgmatau * dTau;
    
  if (track.TrackPlane[LON] && track.EnableRW_short_LON)
  {    
    double ARW = sqrt(Z_0 / (sigmarw*C_LIGHT*pi)) / (4*beffL1*pi) * ring->Lc;   
    for(nc=1; nc<ncellmax; nc++) 
    {
      double tRW    = (nc+0.5)*t_const;
      SelfFieldModel->Gl[nc] -= ARW/(tRW*sqrt(tRW));     
    }
    double * RW0_val;
    RW0_val = (double *) malloc(1 * sizeof(double));
    s0 = pow((2 * beffL2 * ring->rhorw / Z_0), 1.0/3.0);
    if(RW_table_LON(sgmatau*dTau * C_LIGHT/s0, RW0_val) < 0)
      ERROR("RW_table_LON", return false);

    SelfFieldModel->Gl[0] += *(RW0_val) * 4*Z_0*C_LIGHT / (pi*beffL2) * ring->Lc;    
    SelfFieldModel->PlaneL ++;
   } /* end LON */

   if (track.TrackPlane[VER])
   {
    double ARW = sqrt(Z_0*C_LIGHT/sigmarw/pi)/(pi*beffV3)*ring->Lc;
    double * RW0_val;
    RW0_val = (double *) malloc(1 * sizeof(double));
    s0 = pow((2 * ring->beffV[1] * ring->beffV[1] * ring->rhorw / Z_0), 1.0/3.0);
    if(RW_table_TRANS(sgmatau*dTau * C_LIGHT/s0, RW0_val) < 0)
      ERROR("RW_table_VER", return -1);
    SelfFieldModel->Gv[0] += *(RW0_val) * 8*Z_0*C_LIGHT / (pi*beffV4) * s0 * ring->Lc;
    for(nc=1; nc<ncellmax; nc++) 
    {
      double tRW    = (nc+0.5)*t_const;
      SelfFieldModel->Gv[nc] += ARW/sqrt(tRW);
    }   
    SelfFieldModel->PlaneV++;
  }

  if (track.TrackPlane[HOR])
  {   
    double ARW = sqrt(Z_0*C_LIGHT/sigmarw/pi)/(pi*beffH3)*ring->Lc;
    double * RW0_val;
    RW0_val = (double *) malloc(1 * sizeof(double));
    s0 = pow((2 * ring->beffH[1] * ring->beffH[1] * ring->rhorw / Z_0), 1.0/3.0);
    if(RW_table_TRANS(sgmatau*dTau * C_LIGHT/s0, RW0_val) < 0)
      ERROR("RW_table_HOR", return -1);
    SelfFieldModel->Gh[0] += *(RW0_val) * 8*Z_0*C_LIGHT / (pi*beffH4) * s0 * ring->Lc;
    for(nc=1; nc<ncellmax; nc++) 
    {
      double tRW    = (nc+0.5)*t_const;
      SelfFieldModel->Gh[nc] += ARW/sqrt(tRW);
     }
     SelfFieldModel->PlaneH ++;
  } /* end HOR */
  return 1;
} /* end fuction */




int
transform_weak_bunch_selffield(weak_bunch_t * bunch,
                               const selffield_model_t SelfFieldModel, 
                               const cyclic_array_t * all_moments, const ring_t *ring,
                               long unsigned int rev, FILE * fp, double scan_val, int kb, 
			       e_beam_t * ebeam, double * phasor_end,  double * fnp_ring)

{
  unsigned jp;
  int icell, k;  
  int * mapcell = (int *) malloc(bunch->Np * sizeof(int));  
  double sgm_xtau = SelfFieldModel.sigma_tau; /* [s] */
  double sgmatau2 = sgm_xtau * sgm_xtau; //for purely inductive
  double dTau2    = SelfFieldModel.dT * SelfFieldModel.dT; //for purely inductive
  double factG = -ring->T0 * bunch->Ib / (bunch->Np * ring->E0 * FGIGA); /* Wakes are in [V/C] or [V/Cm] -> per unit length & for a test charge with Q = 1 C */
  double Tau, Gn;
  double fnp[SelfFieldModel.Ncell], GL1[SelfFieldModel.Ncell];
  double dipoleV[SelfFieldModel.Ncell], GlambdaV[SelfFieldModel.Ncell];
  double dipoleH[SelfFieldModel.Ncell], GlambdaH[SelfFieldModel.Ncell];
  double lr_wake[SelfFieldModel.Ncell];
  
  /*---------------------------------------------------------------------------------------------------*/
  /*      The way how the Greens function and the wake potential are computed follows completely       
   *      that developed in the longitudinal tracking code of G. Besnier.      (April 2005)            
   *      Modified the routine to handle different impedance components:       (September 2005)        
   *         - Purely resistive                                                                        
   *         - Purely inductive                                                                        
   *         - BBR                                                                                     
   *      The Green's function constructed is the sum of all components, and the way it is             
   *      defined for different types of impedance, as well as its convolution with beam follows       
   *      that of G. Besnier (its correctness confirmed).                                              
   *      Impedance components treated here follow "sbtrack_Zlong.inp"                                 
   *      <NB> To handle the sum of many different impedance components, the definition of Gl[]        
   *           has been changed (The original Gl[] is multiplied by RL/sgm_xtau, to make the           
   *           multiplication factor factG1 independent of the impedance (19/09/2005)                  */
  /*---------------------------------------------------------------------------------------------------*/
  /*    Update of the density distribution:                                                            */
  /*---------------------------------------------------------------------------------------------------*/
  
  memset(fnp, 0.0, sizeof(fnp));
  memset(dipoleV, 0.0, sizeof(dipoleV));
  memset(dipoleH, 0.0, sizeof(dipoleH));
  memset(GL1, 0.0, sizeof(GL1));
  memset(GlambdaV, 0.0, sizeof(GlambdaV));
  memset(GlambdaH, 0.0, sizeof(GlambdaH));
  memset(lr_wake, 0.0, sizeof(lr_wake)); 
  
  /* Assigning a bin to every particle */
  int icellmin, icellmax;
  particle_t * particle;
  
  for (jp = 0; jp < bunch->Np; jp++)
  {
    particle = &(bunch->particles[jp]);
    Tau = -particle->pos.xtau / sgm_xtau; // units of time!
    /*******  <NB> Corrected the sign xtau0[jp] --> -xtau0[jp] in Sep98  ******/
    icell = (int) ((Tau + SelfFieldModel.Nsigma)/SelfFieldModel.dT + 1);
    
    if(icell < 0) 
    {
      icell = 0;  // particles out of binning are put in first bin
      bunch->N_trash_low ++;
    }
    else if(icell >= SelfFieldModel.Ncell)
    {
      icell = SelfFieldModel.Ncell - 1;  // particles out of binning are put in last bin
      bunch->N_trash_high ++;
    }
    
    mapcell[jp] = icell; //asinges a cell to every macroparticle      
    fnp[icell] += 1.0; //counts macroparticles per cell  
    
    dipoleV[icell] += particle->pos.z;
    dipoleH[icell] += particle->pos.x;
  }
  
  icellmin = SelfFieldModel.Ncell-1;
  icellmax = 0;
  for (icell = 0; icell < icellmin; icell++)
  {
    if(fnp[icell]>0)
    {
      icellmin = 1*icell;
      break;
    } 
  }     
  
  for (icell = SelfFieldModel.Ncell-1; icell > icellmax; icell--)
  {
    if(fnp[icell] > 0)
    {
      icellmax = 1*icell;
      break;
    }  
  }
  
  //**********************************************************************************//
  if(SelfFieldModel.PlaneL > 0)
  {
    /*  Construction of wake potential:      */
    double ntmoy, ntmoy0, ntmoy1; /* smoothen values for local density by averaging over neichboring cells */
    for (icell = icellmin; icell <= SelfFieldModel.Ncell; icell++)
    {
      /***  Contribution of BBR impedances  ***/
      Gn = 0.0;
      for (k = icell; k >= icellmin; k--)
      {
        Gn += fnp[k] * SelfFieldModel.Gl[icell - k];
      }
      GL1[icell] = 1*Gn;
    }
    
    /***  Contribution of purely resistive impedances  ***/
    if(SelfFieldModel.RsisL > 0)
      for (icell = icellmin+1; icell <= icellmax; icell++)
      {
        ntmoy = (0.5 * fnp[icell-1] + fnp[icell] + 0.5 * fnp[icell+1]) / 2.0;
        GL1[icell] += SelfFieldModel.RsisL * ntmoy/(sgm_xtau * SelfFieldModel.dT);
      }
      /***  Contribution of purely inductive impedances  ***/
      if(SelfFieldModel.aindL > 0)
        for (icell = icellmin+1; icell <= icellmax; icell++)
        {
          ntmoy0 = (0.25 * fnp[icell-3] + 0.5 * fnp[icell-2] + fnp[icell-1] + 0.5 * fnp[icell]  + 0.25 * fnp[icell+1]) / 2.5;
          ntmoy1 = (0.25 * fnp[icell-1] + 0.5 * fnp[icell]   + fnp[icell+1] + 0.5 * fnp[icell+2] + 0.25 * fnp[icell+3]) / 2.5;
          GL1[icell] += (ntmoy1 - ntmoy0) / (2.0 * sgmatau2 * dTau2) * SelfFieldModel.aindL;
        }
  }
  
      //Construction of the longrange wake potential:
  if(ring->longrange_resonators_size > 0)
    construct_wake_phasor(lr_wake, phasor_end, &SelfFieldModel,
                          ring, kb, fnp_ring, ebeam, bunch);    
  
  if(ring->longrange_resonators_size + SelfFieldModel.PlaneL > 0)
  {
    // Total effect of wake potentials
    for (jp = 0; jp < bunch->Np; jp++)
      bunch->particles[jp].slope.xtau += factG * GL1[mapcell[jp]] + lr_wake[mapcell[jp]];
  }
  
  
  //**********************************************************************************//
  if(SelfFieldModel.PlaneV > 0)
  {
    for (icell = 0; icell < SelfFieldModel.Ncell; icell++)
    {
      Gn = 0.0;
      for (k = icell; k >= 0; k--)
        Gn += dipoleV[k] * SelfFieldModel.Gv[icell - k];
      GlambdaV[icell] = 1*Gn;
    }    
    for (jp = 0; jp < bunch->Np; jp++)
      bunch->particles[jp].slope.z += -factG * GlambdaV[mapcell[jp]];
  }
  
  //**********************************************************************************//
  if(SelfFieldModel.PlaneH > 0)
  {
    for (icell = 1; icell < SelfFieldModel.Ncell; icell++)
    {
      Gn = 0.0;
      for (k = icell; k >= 0; k--)
        Gn += dipoleH[k] * SelfFieldModel.Gh[icell - k];
      GlambdaH[icell] = 1*Gn;
    }    
    for (jp = 0; jp < bunch->Np; jp++)
      bunch->particles[jp].slope.x += -factG * GlambdaH[mapcell[jp]];
  }
  
  //**********************************************************************************//
  
  if(bunch->kb_out == 1 && (rev+1)%track.NrevMon == 0 && rev+1>=track.NrevPotentialsOut+(rev/track.NrevScan)*track.NrevScan) /* Output potentials */
  {
    double rftest, taucell;
    const double dTau = SelfFieldModel.dT;
    double active, FacI, ideal;
    double eps_const = ring->Vrf0 / (ring->E0 * FKILO);
    int j;
    active_HC_t * aHC;
    if(track.EnableIdealHC == 0)
      FacI = 0.;
    else
      FacI = eps_const *  ring->HC_k;
    
    
    for (icell = 0; icell < SelfFieldModel.Ncell; icell++)
    {
      taucell = -(-SelfFieldModel.Nsigma + dTau * (icell - 0.5))* sgm_xtau;
      rftest = eps_const * sin(ring->wrf * taucell + ring->phai0) - ring->Urad;
      // adding potential of active HCs
      active = 0.0;
      for(j = 0; j < ring->active_HC_size; j++)
      {
        aHC = &(ring->active_HC[j]);
        active += aHC->Vpeak/(ring->E0 * FGIGA) * sin(ring->wrf * aHC->nHC * taucell + aHC->phi_aHC + aHC->nHC * kb * 2 * M_PI);
      }
      ideal = FacI * sin(ring->m_aHC * ring->wrf * taucell + ring->m_aHC * ring->phi_n);
      fprintf(fp, "\n %ld   %e   %e   %e    %e   %e   %e   %e   %e",rev+1, scan_val, taucell,  factG * GL1[icell], lr_wake[icell], rftest, active, ideal, fnp[icell]);
      if(SelfFieldModel.PlaneV > 0)
        fprintf(fp, "   %e   %e", -factG * GlambdaV[icell], dipoleV[icell]);
      if(SelfFieldModel.PlaneH > 0)
        fprintf(fp, "   %e   %e", -factG * GlambdaH[icell], dipoleH[icell]);
    }
    fprintf(fp, "\n");
  }
  free(mapcell);
  return 1;
}



//****************************************************************//
//************** RW longrange section ****************************//
//****************************************************************//

int
transform_weak_bunch_RW_longrange_cyclic(const int in, const int bpos,
                                  cyclic_array_t * all_moments,
                                  const plane_t plane,
                                  const bunch_macroparticle_model_t * MPmodel,
                                  e_beam_t * ebeam,
                                  weak_bunch_t * bunch,
                                  ring_t * ring,
                                  tracking_t * track)
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
  
  const int h = ring->Nharm;
  double beff3 = plane == 1 ? ring->beffH[0] : ring->beffV[0];
  const int Ntrigger = 50;
  if(track->triggerRW && in < Ntrigger)
      beff3 /=  2;
  
  int m, k, in_mod;
  double tt;
  const double * moments;
  double wake_voltage = 0.;
  double deltab = ring->T0 / ring->h;
  double RWconst = ring->T0 / ring->E0 / MPmodel->fNp / FTERA * ring->Lc / (M_PI * pow(beff3,3)) * sqrt(Z_0 * C_LIGHT * ring->rhorw / M_PI);  
  
  for (m = 1; m < ((track->Nmlt+1)*ring->Nharm  + 1); m++)
  {
    k = bpos + m; // Bunch index in cyclic array (not bounded)           
    // Get dipole moment of bunch kb + m which affects bunch kb 
    in_mod = k < 0 ? (k%h + h)%h : k%h;    
    if(ebeam->nfFill[in_mod])
    {
      moments = cyclic_array_get(k, all_moments);
      tt = m * deltab;        
      wake_voltage += MPmodel->fNp * moments[plane-1] * ring->Ibunch[in_mod] / sqrt(tt);
    }  
  }
  wake_voltage *= RWconst;  

  
  int jp;
  for(jp = 0; jp < bunch->Np; jp++)
      bunch->particles[jp].slope.v[plane] += wake_voltage;

  return 1;
}



 /* determines macroparicles per bin before first turn  for bunch kb*/
 int
 fnp_ring_update(ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel, weak_bunch_t * bunch, int kb)
 {
   int jp, icell;
   int Nbin = SelfFieldModel->Ncell;
   double Tau;
   
   for(icell=0; icell < Nbin; icell++)
     fnp_ring[kb*Nbin + icell] = 0.0;     
     

   particle_t * particle;
   for (jp = 0; jp < bunch->Np; jp++)
   {
     particle = &(bunch->particles[jp]);
     Tau = -particle->pos.xtau / SelfFieldModel->sigma_tau; // units of time!
     /*******  <NB> Corrected the sign xtau0[jp] --> -xtau0[jp] in Sep98  ******/
     icell = (int) ((Tau + SelfFieldModel->Nsigma)/SelfFieldModel->dT + 1);
     
     if(icell >= 0 &&  icell < Nbin)
     {       
        fnp_ring[kb*Nbin + icell] += 1.0; //counts macroparticles per cell     
     }
   } 
   return 1;
 }
 

 void
 wake_phasor_init(ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel, 
		  weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam)
 {
   int i, j, k, l, m;
   double V_old[2], V_new[2], progress[2], C[2], prog2beam[2], progb2beam[2];
   
   // Determines the phasor after all bunches have Nturn-times passed
   for(l = 0; l < ring->longrange_resonators_size; l++)
   {
     
     LR_resonator_t * lr_res = &(ring->longrange_resonators[l]);   
     double alpha =  lr_res->wr * 0.5 / lr_res->Qfactor;
     double Amp = 2 * alpha * lr_res->Rs * ring->T0 / ring->E0 / FGIGA;
     double dTau = SelfFieldModel->dT * SelfFieldModel->sigma_tau; // Bin-width
     double tbucket = (ring->T0 / ring->h - SelfFieldModel->Ncell * dTau); // distance between two bunches
     int Nbin = SelfFieldModel->Ncell;
     double fac = Amp / bunch->Np * FMILLI;
     int Nturns = (int)(ring->Nbumax / ring->h); // 10 damping times
     
     C[0] = -alpha;
     C[1] = lr_res->wr;
     
     // Decay and oscillation of the phasors during the passage of 1 bin with length dTau
     progress[0] = exp(C[0] * dTau) * cos(C[1] * dTau);
     progress[1] = exp(C[0] * dTau) * sin(C[1] * dTau);
     
     V_old[0] = 0.0;
     V_old[1] = 0.0;
     V_new[0] = phasor_end[l*2];
     V_new[1] = phasor_end[l*2 + 1];

     double taubeam=0;
     double taub2beam = 0;
     if (ring->has_rf_feedback && ring->rf_feedback->lr_resonator==l+1) {
       taubeam = SelfFieldModel->Nsigma*SelfFieldModel->sigma_tau+dTau*1.5;
       taub2beam = -SelfFieldModel->Nsigma*SelfFieldModel->sigma_tau+dTau*0.5-tbucket;
       //if (Nbin/2==Nbin/2.0) taubeam = dTau/2.0*(Nbin+1);
       //else taubeam = dTau*(Nbin/2+1);
     }
     prog2beam[0] = exp(C[0] * taubeam) * cos(C[1] * taubeam); // Decay and rotation of phasor until synchronous phase
     prog2beam[1] = exp(C[0] * taubeam) * sin(C[1] * taubeam); // (for RF feedback)
     progb2beam[0] = exp(C[0] * taub2beam) * cos(C[1] * taub2beam); // Decay and rotation of phasor until synchronous phase
     progb2beam[1] = exp(C[0] * taub2beam) * sin(C[1] * taub2beam); // (for RF feedback)
     
     for(k=0; k<Nturns; k++)
     {
       for(m=0; m<ring->h; m++)
       {
         i = ring->h - m - 1;
	 if (ring->has_rf_feedback && ring->rf_feedback->lr_resonator==l+1)
	   rffb_get_vrf_phi(ring->rf_feedback,V_new[0]*prog2beam[0]-V_new[1]*prog2beam[1],V_new[0]*prog2beam[1]+V_new[1]*prog2beam[0]);

         if(ebeam->nfFill[i] == 1)
         {
           for(j=0; j<Nbin; j++)      
           {        
             V_new[0] = (V_old[0] * progress[0] - V_old[1] * progress[1]) - ring->Ibunch[i] * fnp_ring[i*Nbin + j] * fac;
             V_new[1] = (V_old[0] * progress[1] + V_old[1] * progress[0]);

             V_old[0] = V_new[0];
             V_old[1] = V_new[1]; 
           }
           // Decay and oscillation of the phasors inbetween two bunches
           V_new[0] = (V_old[0] * exp(C[0]*tbucket)*cos(C[1]*tbucket) - V_old[1] * exp(C[0]*tbucket)*sin(C[1]*tbucket));
           V_new[1] = (V_old[0] * exp(C[0]*tbucket)*sin(C[1]*tbucket) + V_old[1] * exp(C[0]*tbucket)*cos(C[1]*tbucket));
           
           V_old[0] = V_new[0];
           V_old[1] = V_new[1];
         }
         
         else
         { // Decay and rotation of phasor during the passage of an empty buncket
           V_new[0] = (V_old[0] * exp(C[0]*ring->T0 / ring->h)*cos(C[1]*ring->T0 / ring->h) - V_old[1] * exp(C[0]*ring->T0 / ring->h)*sin(C[1]*ring->T0 / ring->h));
           V_new[1] = (V_old[0] * exp(C[0]*ring->T0 / ring->h)*sin(C[1]*ring->T0 / ring->h) + V_old[1] * exp(C[0]*ring->T0 / ring->h)*cos(C[1]*ring->T0 / ring->h));

           V_old[0] = V_new[0];
           V_old[1] = V_new[1];
         }
	 if (ring->has_rf_feedback && ring->rf_feedback->lr_resonator==l+1)
	   rffb_get_vrf_phi(ring->rf_feedback,V_new[0]*progb2beam[0]-V_new[1]*progb2beam[1],V_new[0]*progb2beam[1]+V_new[1]*progb2beam[0]);
       }
     }
     phasor_end[l*2] = V_new[0];
     phasor_end[l*2 + 1] = V_new[1];    
   }
 }
 
 
 /**
  * 
  * @param lr_wake Stores wake function of all bins from several longrange resonators
  * @param phasor_end Stores phasors at the end of the last turn
  */
 void
 construct_wake_phasor(double * lr_wake, double * phasor_end, 
                       const selffield_model_t * SelfFieldModel,
                       const ring_t * ring, int kb, const double * fnp, 
		       e_beam_t * ebeam, weak_bunch_t * bunch)
 {
   int i, j, l, m;
   double progress[2];
   double C[2];
   double V_old[2];
   double V_new[2];
   double prog2beam[2];
   double progb2beam[2];
   
   // Determines the phasor after all bunches have passed and stores the sum in lr_wake
   for(l = 0; l < ring->longrange_resonators_size; l++)
   {     
     LR_resonator_t * lr_res = &(ring->longrange_resonators[l]);
     double alpha =  lr_res->wr * 0.5 / lr_res->Qfactor;
     double Amp = 2 * alpha * lr_res->Rs * ring->T0 / ring->E0 / FGIGA;
     double dTau = SelfFieldModel->dT * SelfFieldModel->sigma_tau; // Bin-width
     double tbucket = (ring->T0 / ring->h - SelfFieldModel->Ncell * dTau); // distance between two bunches
     int Nbin = SelfFieldModel->Ncell;   
     double fac = Amp / bunch->Np * FMILLI;
     
     C[0] = -alpha;
     C[1] = lr_res->wr;
     progress[0] = exp(C[0] * dTau) * cos(C[1] * dTau); // Decay and rotation of phasor during one bin
     progress[1] = exp(C[0] * dTau) * sin(C[1] * dTau);
     
     V_old[0] = phasor_end[l*2];
     V_old[1] = phasor_end[l*2 + 1];
     V_new[0] = phasor_end[l*2];
     V_new[1] = phasor_end[l*2 + 1];
     double taubeam = 0;
     double taub2beam = 0;
     if (ring->has_rf_feedback && ring->rf_feedback->lr_resonator==l+1) {
       taubeam = SelfFieldModel->Nsigma*SelfFieldModel->sigma_tau+dTau*1.5;
       taub2beam = -SelfFieldModel->Nsigma*SelfFieldModel->sigma_tau+dTau*0.5-tbucket;
       //if (Nbin/2==Nbin/2.0) taubeam = dTau/2.0*(Nbin+1);
       //else taubeam = dTau*(Nbin/2+1);
     }
     prog2beam[0] = exp(C[0] * taubeam) * cos(C[1] * taubeam); // Decay and rotation of phasor until synchronous phase
     prog2beam[1] = exp(C[0] * taubeam) * sin(C[1] * taubeam); // (for RF feedback)
     progb2beam[0] = exp(C[0] * taub2beam) * cos(C[1] * taub2beam); // Decay and rotation of phasor until synchronous phase
     progb2beam[1] = exp(C[0] * taub2beam) * sin(C[1] * taub2beam); // (for RF feedback)
     
     for(m=0; m<ring->h; m++)
     {
       i = ring->h - m - 1;
       if (ring->has_rf_feedback && ring->rf_feedback->lr_resonator==l+1)
	 rffb_get_vrf_phi(ring->rf_feedback,V_new[0]*prog2beam[0]-V_new[1]*prog2beam[1],V_new[0]*prog2beam[1]+V_new[1]*prog2beam[0]);

       if(ebeam->nfFill[i] == 1)
       {
         for(j=0; j<Nbin; j++)
         {

           V_new[0] = (V_old[0] * progress[0] - V_old[1] * progress[1]) - ring->Ibunch[i] * fnp[i*Nbin + j] * fac;
           V_new[1] = (V_old[0] * progress[1] + V_old[1] * progress[0]);

           if(i == kb)
           {
             lr_wake[j] += V_new[0]+0.5*ring->Ibunch[i]*fnp[i*Nbin+j]*fac; // Phasor of actual resonator l is added 
           }

           V_old[0] = V_new[0];
           V_old[1] = V_new[1];
         }
         V_new[0] = (V_old[0] * exp(C[0]*tbucket)*cos(C[1]*tbucket) - V_old[1] * exp(C[0]*tbucket)*sin(C[1]*tbucket)); // Decay and rotation of phasor between bunches
         V_new[1] = (V_old[0] * exp(C[0]*tbucket)*sin(C[1]*tbucket) + V_old[1] * exp(C[0]*tbucket)*cos(C[1]*tbucket));

         V_old[0] = V_new[0];
         V_old[1] = V_new[1]; 
       }
       else
       { // Decay and rotation of phasor during the passage of an empty buncket
         V_new[0] = (V_old[0] * exp(C[0]*ring->T0 / ring->h)*cos(C[1]*ring->T0 / ring->h) - V_old[1] * exp(C[0]*ring->T0 / ring->h)*sin(C[1]*ring->T0 / ring->h));
         V_new[1] = (V_old[0] * exp(C[0]*ring->T0 / ring->h)*sin(C[1]*ring->T0 / ring->h) + V_old[1] * exp(C[0]*ring->T0 / ring->h)*cos(C[1]*ring->T0 / ring->h));

         V_old[0] = V_new[0];
         V_old[1] = V_new[1];
       }
       if (ring->has_rf_feedback && ring->rf_feedback->lr_resonator==l+1)
	 rffb_get_vrf_phi(ring->rf_feedback,V_new[0]*progb2beam[0]-V_new[1]*progb2beam[1],V_new[0]*progb2beam[1]+V_new[1]*progb2beam[0]);
     }    
     phasor_end[l*2] = V_new[0];
     phasor_end[l*2 + 1] = V_new[1];   
   }
 }

void
fnp_ring_destroy(double * fnp_ring)
 {
   free(fnp_ring);
   fnp_ring = NULL;
 }
 
 
 
 
 
 
 
 
 
