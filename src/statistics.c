#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "statistics.h"
#include "bunch.h"

vector_t
vector_new(const double xtau, const double x, const double z)
{
  const vector_t v = {{xtau, x, z}};
  return v;
}

vector_t
vector_add(const vector_t a, const vector_t b)
{
  const vector_t c = {{a.xtau + b.xtau, a.x + b.x, a.z + b.z}};
  return c;
}

vector_t
vector_mul_scalar(const vector_t a, const double k)
{
  const vector_t v = {{k * a.xtau, k * a.x, k * a.z}};
  return v;
}

vector_t
vector_square(const vector_t a)
{
  const vector_t v = {{a.xtau * a.xtau, a.x * a.x, a.z * a.z}};
  return v;
}

vector_t
vector_variance(const vector_t a, const vector_t b_sqr, const long int N)
{
  vector_t a_sqr = vector_square(a);
  a_sqr = vector_mul_scalar(a_sqr, -1.0/((double)N));
  vector_t v = vector_add(a_sqr, b_sqr);
  v = vector_mul_scalar(v, 1.0/((double)(N-1)));
  return v;
}


double weighted_rms(const double * in, const int * weights, const size_t n, const double avr)
{
  unsigned i;
  double sqsum = 0.0;
  int weightsum = 0;
  for(i = 0; i < n; i++)
  {
    sqsum += (in[i] - avr)*(in[i] - avr) * weights[i];
    weightsum += weights[i];
  }
  return sqrt(sqsum / ((double) weightsum));
}

int histogram_generate(const double * input, const unsigned n,
                       const double middle, const unsigned m, const double width,
                       unsigned * hist)
{
  if(input == NULL || hist == NULL)
  {
    return -1;
  }
  else
  {
    const double hist_left = middle - width * ((double) m) / 2.0;
    unsigned ibin = 0;
    for(ibin = 0; ibin < m; ibin++)
    {
      const double bin_left = hist_left + ((double) ibin) * width;
      const double bin_right = hist_left + ((double) ibin + 1) * width;
      unsigned i = 0;
      hist[ibin] = 0;
      for(i = 0; i < n; i++)
      {
        if(input[i] >= bin_left && input[i] < bin_right)
        {
          hist[ibin] += 1;
        }
      }
    }
    return 1;
  }
}

void
weak_bunch_mean_rms(weak_bunch_t * bunch, unsigned int plane)
{
  unsigned int jp;
  bunch_stats_t * bstats = &(bunch->stats);
  double pos_sum = 0.0;
  double slope_sum = 0.0;
  double pos_ave;
  double slope_ave;
  particle_t * particle;

  /* Average */
  pos_sum = 0.0;
  slope_sum = 0.0;
  for(jp = 0; jp < bunch->Np; jp++)
  {
    particle = &(bunch->particles[jp]);
    pos_sum += particle->pos.v[plane];
    slope_sum += particle->slope.v[plane];
  }   
  pos_ave = pos_sum/((double) bunch->Np);
  slope_ave = slope_sum/((double) bunch->Np);
  
  bstats->pos.v[plane] = pos_ave;  
  bstats->slope.v[plane] = slope_ave;
  
  /* RMS */
  pos_sum = 0.0;
  slope_sum = 0.0;
  for(jp = 0; jp < bunch->Np; jp++)
  {
    particle = &(bunch->particles[jp]);
    pos_sum += pow(pos_ave - particle->pos.v[plane], 2);
    slope_sum += pow(slope_ave - particle->slope.v[plane], 2);
  }
  
  bstats->pos_sigma.v[plane] = sqrt(pos_sum / ((double) bunch->Np));
  bstats->slope_sigma.v[plane] = sqrt(slope_sum / ((double) bunch->Np));
}

void
weak_bunch_calc_statistics(weak_bunch_t * bunch,  const tracking_t * track)
{
  unsigned int plane;
  
  if(bunch == NULL || bunch->particles == NULL)
    return;

  bunch_stats_t * bstats = &(bunch->stats);
  
  for(plane = 0; plane < 3; plane ++)
  {
    if(track->TrackPlane[plane])
    {
      weak_bunch_mean_rms(bunch,plane);

      bstats->sum_pos.v[plane] += bstats->pos.v[plane];
      bstats->sum_pos_sigma.v[plane] += bstats->pos_sigma.v[plane];
      bstats->sum_slope.v[plane] += bstats->slope.v[plane];
      bstats->sum_slope_sigma.v[plane] += bstats->slope_sigma.v[plane];
        
      bstats->sum_pos_sqr.v[plane] += pow(bstats->pos.v[plane], 2.);
      bstats->sum_pos_sigma_sqr.v[plane] += pow(bstats->pos_sigma.v[plane], 2.);
      bstats->sum_slope_sqr.v[plane] += pow(bstats->slope.v[plane], 2.);
      bstats->sum_slope_sigma_sqr.v[plane] += pow(bstats->slope_sigma.v[plane], 2.);
    }
  }
}




void
weak_bunch_calc_ampinv(weak_bunch_t * bunch,  const tracking_t * track)
{
  // SNAPSHOT of ampinv
  
  unsigned int jp, plane;
  double ampinv_sum;
  
  if(bunch == NULL || bunch->particles == NULL)
    return;

  bunch_stats_t * bstats = &(bunch->stats);
  particle_t * particle;
  
  for(plane = 0; plane < 3; plane ++)
  {
    if(track->TrackPlane[plane])
    { 
      /* RMS, ampinv */
      ampinv_sum = 0.0;
      for(jp = 0; jp < bunch->Np; jp++)
      {
        particle = &(bunch->particles[jp]);
        ampinv_sum += ring.gamma1[plane] * pow(particle->pos.v[plane], 2)
                  + 2.0 * ring.alpha1[plane] * particle->pos.v[plane] * particle->slope.v[plane]
                  + ring.beta1[plane] * pow(particle->slope.v[plane], 2);
      }
      bstats->ampinv.v[plane] = ampinv_sum/((double) bunch->Np);
      bstats->ampinv_cm.v[plane] = ring.gamma1[plane] * pow(bstats->pos.v[plane], 2)
          + 2.0 * ring.alpha1[plane] * bstats->pos.v[plane] * bstats->slope.v[plane]
          + ring.beta1[plane] * pow(bstats->slope.v[plane], 2);
    }
  }
}







void
weak_bunch_add_statistics(bunch_stats_t * allstats, const bunch_stats_t * addstats)
{
  allstats->sum_pos = vector_add(allstats->sum_pos, addstats->pos);
  allstats->sum_pos_sigma = vector_add(allstats->sum_pos_sigma, addstats->pos_sigma);
  allstats->sum_slope = vector_add(allstats->sum_slope, addstats->slope);
  allstats->sum_slope_sigma = vector_add(allstats->sum_slope_sigma, addstats->slope_sigma);
        
  vector_t v_sqr;
  v_sqr = vector_square(addstats->pos);
  allstats->sum_pos_sqr = vector_add(allstats->sum_pos_sqr, v_sqr);
  v_sqr = vector_square(addstats->pos_sigma);
  allstats->sum_pos_sigma_sqr = vector_add(allstats->sum_pos_sigma_sqr, v_sqr);
  v_sqr = vector_square(addstats->slope);
  allstats->sum_slope_sqr = vector_add(allstats->sum_slope_sqr, v_sqr);
  v_sqr = vector_square(addstats->slope_sigma);
  allstats->sum_slope_sigma_sqr = vector_add(allstats->sum_slope_sigma_sqr, v_sqr);
}
  
void
weak_bunch_update_statistics(bunch_stats_t * bstats, const long int Nstat)
{
  if(Nstat > 1)
  {
    /* variances */
    bstats->sum_pos_sqr = vector_variance(bstats->sum_pos, bstats->sum_pos_sqr, Nstat);
    bstats->sum_pos_sigma_sqr = vector_variance(bstats->sum_pos_sigma, bstats->sum_pos_sigma_sqr, Nstat);
    bstats->sum_slope_sqr = vector_variance(bstats->sum_slope, bstats->sum_slope_sqr, Nstat);
    bstats->sum_slope_sigma_sqr = vector_variance(bstats->sum_slope_sigma, bstats->sum_slope_sigma_sqr, Nstat);
  }
  else
  {
    bstats->sum_pos_sqr = vector_new(0.0, 0.0, 0.0);
    bstats->sum_pos_sigma_sqr = vector_new(0.0, 0.0, 0.0);
    bstats->sum_slope_sqr = vector_new(0.0, 0.0, 0.0);
    bstats->sum_slope_sigma_sqr = vector_new(0.0, 0.0, 0.0);  
  }
    /* means */
    bstats->sum_pos = vector_mul_scalar(bstats->sum_pos, 1.0/((double) Nstat));
    bstats->sum_pos_sigma = vector_mul_scalar(bstats->sum_pos_sigma, 1.0/((double) Nstat));
    bstats->sum_slope = vector_mul_scalar(bstats->sum_slope, 1.0/((double) Nstat));
    bstats->sum_slope_sigma = vector_mul_scalar(bstats->sum_slope_sigma, 1.0/((double) Nstat));
}


void
weak_bunch_reset_statistics(bunch_stats_t * bstats)
{    
  bstats->sum_pos = vector_new(0.0, 0.0, 0.0);
  bstats->sum_pos_sigma = vector_new(0.0, 0.0, 0.0);
  bstats->sum_slope = vector_new(0.0, 0.0, 0.0);
  bstats->sum_slope_sigma = vector_new(0.0, 0.0, 0.0);
        
  bstats->sum_pos_sqr = vector_new(0.0, 0.0, 0.0);
  bstats->sum_pos_sigma_sqr = vector_new(0.0, 0.0, 0.0);
  bstats->sum_slope_sqr = vector_new(0.0, 0.0, 0.0);
  bstats->sum_slope_sigma_sqr = vector_new(0.0, 0.0, 0.0);
}

