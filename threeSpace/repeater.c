#include <stdio.h>
#include <math.h>
#include "wls.h"

/*************************************************************************

  Name: repeater
  
  Purpose: The purpose of this program is to modify time difference of
  arrival (TDOA) information necessary for use of the wls program that
  is used for RF emitter location determination.  The position parameters
  have units of feet or meters, and the time parameters have units of
  nanoseconds. It should be noted that only the time difference information
  is modified. All other parameters of the file are unchanged.
  
  Calling Sequence: repeater [-a] [-f] < inputFile >outputFile
  
  Inputs:
  
    a - A flag which indicates that repeater path delay should be added
    to the time differences.  If this flag is not specified, repeater
    path delay will be subtracted from the time differences.
    
    An optional flag (the -f) which indicates that position parameters are
    in feet, otherwise, they are in units of meters.
       
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
  
  Outputs:
  
    A text file with the same format as the input file with adjusted
    time differences.

*************************************************************************/
int main(int argc, char **argv)
{
  char buffer[80];
  int i, n;
  double x[50], y[50], z[50], dt;
  double deltax, deltay, deltaz;
  double t;
  int addRepeaterPathDelay;
  int spaceIsInFeet;

  switch (argc)
  {
    case 3:
      addRepeaterPathDelay = 1;
      spaceIsInFeet = 1;
      break;
    
    case 2:
      if (argv[1][1] == 'd')
      {
        addRepeaterPathDelay = 1;
	spaceIsInFeet = 0;
      } /* if */
      else
      {
        addRepeaterPathDelay = 0;
        spaceIsInFeet = 1;
      } /* else */
      break;
      
    default:
      addRepeaterPathDelay = 0;
      spaceIsInFeet = 0;
      break;
  } /* case */
  
  /* read number of receivers, and echo the line */  
  fgets(buffer,sizeof(buffer),stdin);
  sscanf(buffer,"%d",&n);
  printf("%s",buffer);
        
  for (i = 0; i < n; i++)
  {
    /* get location and time difference of observer i */
    fgets(buffer,sizeof(buffer),stdin);
    sscanf(buffer,"%lf %lf %lf %lf\n",&x[i],&y[i], &z[i],&dt);
      
    /* compute difference vector */
    deltax = x[i] - x[0];
    deltay = y[i] - y[0];
    deltaz = z[i] - z[0];
  
    /* compute distance */
    t = (deltax * deltax) + (deltay * deltay) + (deltaz * deltaz);
    t = sqrt(t);

    /* convert to time */
        if (spaceIsInFeet)
    {
      t /= TDOA_FT_PER_NS;
    } /* if */
    else
    {
      t /= TDOA_M_PER_NS;
    } /* else */

    if (addRepeaterPathDelay == 1)
    {
      /* add repeater path delay for testing purposes */
      dt += t;
    } /* if */
    else
    {
      /* subtract repeater path delay for testing purposes */
      dt -= t;
    } /* else */
            
    /* output position and corrected time difference results */
    printf("%f %f %f %20.14f\n",x[i],y[i],z[i],dt);
  } /* for */

  exit(0);

} /* main */
