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

/*************************************************************************
 function prototypes     
*************************************************************************/
void wlsInit(void);

int wlsComputePosition(int m,double xn,double yn,struct vector3 *v,
                       double *x,double *y);

double wlsComputeNomRange(double xn,double xi,double x0,
                          double yn,double yi,double y0);
			 
void wlsComputeDirection(double xn,double xi,double x0,
                         double yn,double yi,double y0,
			 double *cx,double *cy);
			 
int wlsComputeDeltaXy(int m,double **c,double **dR,double **dXy);

#endif /* __WLS_H__ */
