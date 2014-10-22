#include <float.h>
#include "test.h"
#include "bunch.h"

/* Global variables */
ring_t ring;
tracking_t track;
grid_t fbii_grid;

int main()
{
  weak_bunch_t bunch;
  const unsigned int Np = 1234;
  const double Ib = 1.0;
  unsigned int jp;
  test_init(true);
  
  ASSERT_TRUE(weak_bunch_create(&bunch, Np, Ib));
  ASSERT_EQUAL(Np, bunch.Np);
  ASSERT_EQUAL(Ib, bunch.Ib);
  
  for(jp = 0; jp < bunch.Np; jp++)
  {
    bunch.particles[jp].pos.x = 5.0;
    bunch.particles[jp].pos.z = 6.0;
    bunch.particles[jp].pos.xtau = 7.0;
    bunch.particles[jp].slope.x = 8.5;
    bunch.particles[jp].slope.z = 9.5;
    bunch.particles[jp].slope.xtau = 10.5;
  }
  
  weak_bunch_update_statistics(&bunch, HOR);
  weak_bunch_update_statistics(&bunch, VER);
  weak_bunch_update_statistics(&bunch, LON);
  
  ASSERT_EQUAL_FLOAT(bunch.stats.pos.x, 5.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.pos.z, 6.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.pos.xtau, 7.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.slope.x, 8.5, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.slope.z, 9.5, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.slope.xtau, 10.5, LDBL_EPSILON);
  
  ASSERT_EQUAL_FLOAT(bunch.stats.pos_sigma.x, 0.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.pos_sigma.z, 0.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.pos_sigma.xtau, 0.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.slope_sigma.x, 0.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.slope_sigma.z, 0.0, LDBL_EPSILON);
  ASSERT_EQUAL_FLOAT(bunch.stats.slope_sigma.xtau, 0.0, LDBL_EPSILON);

  return test_backend();
}
