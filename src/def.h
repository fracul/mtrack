/**
 * @file
 * Commonly used definitions (constants, vector structure)
 */

#ifndef MBTRACK_DEF_H
#define MBTRACK_DEF_H

/* SI prefixes (floating point number format) */
#define FPICO  1.0e-12
#define FNANO  1.0e-09
#define FMICRO 1.0e-06
#define FMILLI 1.0e-03

#define FKILO  1.0e+03 /* replace fmil */
#define FMEGA  1.0e+06 /* replace fmll */
#define FGIGA  1.0e+09 /* replace fnan */
#define FTERA  1.0e+12 /* replace fpco */

/* Physics constants */
#define C_LIGHT  2.997925e+08  /* Speed of Light [m.s-1]. Replace clight. */
#define EPSILON0 8.854187817620e-12 /* Vacuum permittivity = 1/(C_LIGHT*C_LIGHT*mu0) [F·m−1] */
#define E_CHARGE 1.602e-19     /* Electron charge [C] */
#define E_MASS   9.1e-31       /* mass of the electron [kg] */
#define E_RADIUS 2.819e-15     /* Classical eletron radius */
#define MP_O_EM  1.836e3       /* Ratio of proton mass to eletron mass */
#define Z_0 376.73             /* Impedance of free space */

/* Fillings, do not change numerical values */
typedef enum
{
  UNIFORM = 0,
  ONETHIRD = 1,
  TWOTHIRD = 2,
  SINGLE = 3,
  SINGLE2 = 4,
  SINGLE3 = 5,
  THREEFOURTH = 6,
  ONEFOURTH = 7,
  SINGLE10 = 8,
  SINGLE20 = 9,
  SINGLE30 = 10,
  EVENHALF = 11,
  EVENFOURTH = 12
}
filling_t;

/* Filling constants, for compatibility */

#define uniform       UNIFORM
#define onefourth     ONEFOURTH
#define onethird      ONETHIRD
#define twothird      TWOTHIRD
#define threefourth   THREEFOURTH
#define single        SINGLE
#define single2       SINGLE2
#define single3       SINGLE3
#define single10      SINGLE10
#define single20      SINGLE20
#define single30      SINGLE30
#define evenhalf      EVENHALF
#define evenfourth    EVENFOURTH

/* Planes */
typedef enum
{
  LON = 0,
  HOR = 1,
  VER = 2,
  TRA = 1
}
plane_t;


/**
 * \brief Vector of coordinates
 */
typedef union
{
  struct /** 3D vector as a structure */
  {
    double xtau; /**< Longitudinal axis */
    double x; /**< Horizontal axis */
    double z; /**< Vertical axis */
  };
  double v[3]; /** 3D vector as an array */
} vector_t;



/* Error and warning macros */
#include <errno.h>
#include <string.h>

#define ERROR(function_name, instr) \
  { fprintf(stderr, "ERROR: %s failed (%s)\n", function_name, strerror(errno));\
    instr; }

#define WARNING(function_name) \
  fprintf(stderr, "WARNING: %s failed (%s)\n", function_name, strerror(errno));


#endif /* MBTRACK_TYPES_H */
