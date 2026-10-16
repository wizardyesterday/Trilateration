/*************************************************************************
 file name: EmitterLocator.h     
*************************************************************************/

#include <math.h>

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
  EmitterLocator(int n,int spaceIsInFeet);
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
			 
  int computeDeltaXy(int m,double **c,double **deltaR,double **deltaXy);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Attributes.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //-----------------------------------------------
  // Configuration information.
  //-----------------------------------------------
  // meters - 0, feet -1
  int spaceIsInFeet;
  //-----------------------------------------------

  //-----------------------------------------------
  // Distance corrections from  nominal values,
  // x and y, in the linearized distance formula.
  //-----------------------------------------------
  // The matrix.
  double deltaXym[3][2];

  // The pointer to the matrix.
  double *deltaXy[3];
  //-----------------------------------------------

  //-----------------------------------------------
  // Information for difference in differential
  // pseudorange.
  //-----------------------------------------------
  // The matrix.
  double deltaRm[MAX_RECEIVERS][2];

  // The pointer to the matrix.
  double *deltaR[MAX_RECEIVERS];
  //-----------------------------------------------

  //-----------------------------------------------
  // Differential direction cosine information. 
  //-----------------------------------------------
  // The matrix.
  double cm[MAX_RECEIVERS][2];

  // The pointer to the matrix.
  double *c[MAX_RECEIVERS];
  //-----------------------------------------------

  //-----------------------------------------------
  // [C]' information.
  //-----------------------------------------------
  // The matrix.
  double cTm[3][MAX_RECEIVERS];

  // The pointer to the matrix.
  double *cT[3];
  //-----------------------------------------------

  //-----------------------------------------------
  // [C]'[C] information.
  //-----------------------------------------------
  // The matrix.
  double cTxcm[2][2];

  // The pointer to the matrix.
  double *cTxc[2];
  //-----------------------------------------------

  //-----------------------------------------------
  // Inverse([C]'[C]) information.
  //-----------------------------------------------
  // The matrix.
  double cTxcInvm[2][2];

  // The pointer to the matrix.
  double *cTxcInv[2];
  //-----------------------------------------------

  double cTxdeltaRm[3][2];
  double *cTxdeltaR[3];

  //-----------------------------------------------
  // Pseudorange information.
  //-----------------------------------------------
  double *r0i;
  double *rn0i;
  //-----------------------------------------------
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
};

#endif /* __EMITTERLOCATOR_H__ */
