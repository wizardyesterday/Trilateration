/*************************************************************************
 file name: EmitterLocator.h     
*************************************************************************/

#ifndef __EMITTERLOCATOR_H__
#define __EMITTERLOCATOR_H__



/*************************************************************************
 defines    
*************************************************************************/
#define MAX_RECEIVERS   (50)
#define LAMDA_TOLERANCE (0.1)
#define MAX_ITERATIONS  (100)

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
class EmitterLocator
{
  public:
  EmitterLocator(void);
  ~EmitterLocator(void);

  int computePosition(int m,double xn,double yn,struct vector3 *v,
                      double *x,double *y);

  private:

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Utility functions.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  double computeNominalRange(double xn,double xi,double x0,
                             double yn,double yi,double y0);
			 
  void computeDirection(double xn,double xi,double x0,
                        double yn,double yi,double y0,
                        double *cx,double *cy);
			 
  int computeDeltaXy(int m,double **c,double **dR,double **dXy);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Attributes.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  int spaceIsInFeet;

  double dXym[3][2];
  dRm[TRI_MAX_RECEIVERS][2];
  double *dXy[3];
  double *dR[TRI_MAX_RECEIVERS];

  double cm[TRI_MAX_RECEIVERS][3];
  double *c[TRI_MAX_RECEIVERS];

  double cTm[3][TRI_MAX_RECEIVERS];
  double cTxcm[2][2];
  double cTxcInvm[2][2];
  doublecTxdRm[3][2];
  double *cT[3];
  double *cTxc[3];
  double *cTxcInv[3];
  double *cTxdR[3];
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
};

#endif /* __EMITTERLOCATOR_H__ */
