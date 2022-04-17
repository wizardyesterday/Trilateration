#include <stdio.h>
#include <math.h>
#include "wls.h"

/*************************************************************************

  Name: gendata
  
  Purpose: The purpose of this program is to generate time difference of
  arrival (TDOA) information and data necessary for testing of the wls
  program that is used for RF emitter location determination.  The
  position parameters have units of feet or meters, and the time parameters
  have units of nanoseconds.
  
  Calling Sequence: wls < [-f] inputFile >outputFile
  
  Inputs:

    An optional flag (the -f) which indicates that position parameters are
    in feet, otherwise, they are in units of meters.
       
    A text file with the following format:
      n
      x y z
      x0 y0 z0
      x1 y1 z1
      x2 y2 z2
        .
	.
      xn-1 yn-1 zn-1
      	     
      The tuples, [xi yi zi]:i=0,1,2,...n-1 correspond to the locations
      of receivers. Note that the receiver, located at (x0,y0,z0) is designated
      as the reference receiver for which all time differences are computed.
      The tuple [x y z] is the location of the RF emitter. The value n indicates
      the number of receivers that are used for the collection of
      measurements.      
  
  Outputs:
  
    A text file with the following format:
      n
      x0 y0 z0 t00
      x1 y1 z1 t01
      x2 y2 z2 t02
        .
	.
      xm-1 ym-1 zm-1 t0m-1
      	     
      The tuples, [xi yi zi t0i]:i=0,1,2,...m-1 correspond to the locations
      of receivers, to the differences of time of arrival between the
      reference receiver, located at (x0,y0,z0), and a receiver at location
      (xi,yi,zi). Note that the receiver, located at (x0,y0,z0) is designated
      as the reference receiver for which all time differences are computed.
      The value n indicates the number of receivers that were used for the
      collection of measurements      

*************************************************************************/
int main(int argc, char **argv)
{
  char buffer[80];
  int i, n;
  double x, y, z;
  double xo, yo, zo;
  double deltax, deltay, deltaz;
  double t[50];
  int spaceIsInFeet;

  if (argc == 2)
  {
    spaceIsInFeet = 1;
  } /* if */
  else
  {
    spaceIsInFeet = 0;
  } /* else */
  
  /* read number of receivers */
  fgets(buffer,sizeof(buffer),stdin );
  sscanf(buffer,"%d",&n);
       
  /* get position of source */
  fgets(buffer,sizeof(buffer),stdin);
  sscanf(buffer,"%lf %lf %lf",&x,&y,&z);

  /* indicate number of receivers */
  printf("%d\n",n);
  
  for (i = 0; i < n; i++)
  {
    /* get location of observer i */
    fgets(buffer,sizeof(buffer),stdin);
    sscanf(buffer,"%lf %lf %lf",&xo,&yo,&zo);
  
    /* compute difference vector */
    deltax = x - xo;
    deltay = y - yo;
    deltaz = z - zo;
  
    /* compute distance */
    t[i] = (deltax * deltax) + (deltay * deltay) + (deltaz * deltaz);
    t[i] = sqrt(t[i]);
    
    /* convert to time */
    if (spaceIsInFeet)
    {
      t[i] /= TDOA_FT_PER_NS;
    } /* if */
    else
    {
      t[i] /= TDOA_M_PER_NS;
    } /* else */
            
    /* output position and time difference results */
    printf("%f %f %f %20.14f\n",xo,yo,zo,(t[i]-t[0]));
  } /* for */

  exit(0);

} /* main */

