/*************************************************************************
 file name: Trilaterate.h     
*************************************************************************/

#ifndef __TRI_H__
#define __TRI_H__

/*************************************************************************
 defines    
*************************************************************************/
#define TRI_MAX_RECEIVERS   (50)
#define TRI_LAMDA_TOLERANCE (0.1)
#define TRI_MAX_ITERATIONS  (100)

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
void triInit(void);

int triComputePosition(int m,double xn,double yn,struct vector3 *v,
                       double *x,double *y);

double triComputeNomRange(double xn,double xi,double x0,
                          double yn,double yi,double y0);
			 
void triComputeDirection(double xn,double xi,double x0,
                         double yn,double yi,double y0,
			 double *cx,double *cy);
			 
int triComputeDeltaXy(int m,double **c,double **dR,double **dXy);

#endif /* __TRI_H__ */
