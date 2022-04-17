/*************************************************************************
 file name: wls.h     
*************************************************************************/

#ifndef __WLS_H__
#define __WLS_H__

/*************************************************************************
 defines    
*************************************************************************/
#define WLS_MAX_RECEIVERS   (50 + 1)
#define WLS_LAMDA_TOLERANCE (0.1)
#define WLS_MAX_ITERATIONS  (100)

#define TDOA_FT_PER_NS      (0.98357105643L)
#define TDOA_M_PER_NS       (0.299792458L)

/*************************************************************************
 structure definitions     
*************************************************************************/     
struct vector3
{
  double x;   /* x coordinate */
  double y;   /* y coordinate */
  double dt;  /* time difference with respect to reference */  
};

struct vector4
{
  double x;   /* x coordinate */
  double y;   /* y coordinate */
  double z;   /* z coordinate */
  double dt;  /* time difference with respect to reference */  
};

/*************************************************************************
 function prototypes     
*************************************************************************/
void wlsInit3(void);

int wlsComputePosition3(int m,double xn,double yn,double zn,
                        struct vector4 *v,
                        double *x,double *y, double *z);

double wlsComputeNomRange3(double xn,double xi,double x0,
                           double yn,double yi,double y0,
                           double zn,double zi,double z0);
			 
void wlsComputeDirection3(double xn,double xi,double x0,
                          double yn,double yi,double y0,
                          double zn,double zi,double z0,
			  double *cx,double *cy,double *cz);
			 
int wlsComputeDeltaXyz(int m,double **c,double **dR,double **dXyz);

#endif /* __WLS_H__ */
