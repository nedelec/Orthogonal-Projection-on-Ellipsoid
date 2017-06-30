
/**
 These two C-routines calculate the projection P of a point W on the ellipse,
 or on the surface of an Ellipsoid in 3D.
 The ellipse/ellipsoid is aligned with the principal axes (X, Y, Z),
 and has radii (lenX, lenY, lenZ).
 
 Francois Nedelec. Copyright 2007-2017 EMBL.
 
 This code is Open Source covered by the GNU GPL v3.0 License
 */

#include "project_ellipse.h"

#include <math.h>
//#include <stdio.h>

/**
 projectEllipse2D calculates the projection P = (pX, pY) of the point W = (wX, wY)
 on the ellipse that is aligned with the X and Y axis, and has radii (lenX, lenY).
 
 Method:
 
 The normal to the ellipse at position P is N = ( pX / lenX^2, pY / lenY^2 ),
 and we can thus write W = P + h * N, which leads to:
 @code
 pX = wX * lenX^2 / ( lenX^2 + h );
 pY = wY * lenY^2 / ( lenY^2 + h );
 @endcode
 if wX and wY are not both null.
 
 Moreover, the projection should be on the ellipse and thus `h` should be a zero of:
 @code
 F(h) = ( pX / lenX )^2 + ( pY / lenY )^2 - 1
 @endcode
 We follow Newton's rule to find the root of F(h), and use the formula above to
 calculate the projection.
 */

void projectEllipse(double *  pX, double * pY,
                    double    wX, double   wY,
                    double  lenX, double lenY)
{
    // handle special cases:
    if ( wX == 0 )
    {
        *pX = 0;
        *pY = ( wY > 0 ) ? lenY : -lenY;
        return;
    }
    if ( wY == 0 )
    {
        *pX = ( wX > 0 ) ? lenX : -lenX;
        *pY = 0;
        return;
    }
    
    double aa = lenX * lenX;
    double bb = lenY * lenY;
    
    double h = 0, h_old;
    
    // we derive a lower limit for 'h' from  pX^2 + pY^2 > max(lenX,lenY)^2
    double RR = ( bb < aa ) ? aa : bb;
    // 'hmin' is the minimum value that 'h' could have
    double hmin = sqrt( ( wX*wX*aa*aa + wY*wY*bb*bb ) / RR ) - RR;
    
    // we derive another lower limit for 'h' from  |pX| < lenX
    h = ( fabs(wX) - lenX ) * lenX;
    if ( h > hmin )
        hmin = h;

    // we derive another lower limit for 'h' from  |pY| < lenY
    h = ( fabs(wY) - lenY ) * lenY;
    if ( h > hmin )
        hmin = h;

    // if the point is outside, then 'h' should be positive:
    if ( wX*wX/aa + wY*wY/bb > 1 )
    {
        if ( hmin < 0 )
            hmin = 0;
    }
    
    h = hmin;

    // follow Newton's iteration to find the root
    unsigned cnt = 0;
    do {
        double aah = aa + h;
        double bbh = bb + h;

        double waX = wX / aah;
        double waY = wY / bbh;
        
        double pXX = waX * waX * aa;
        double pYY = waY * waY * bb;

        h_old = h;
        
        double F    = 1 - ( pXX         + pYY       );
        double dF   = 2 * ( pXX / aah   + pYY / bbh );

        // Newtons' method
        h -= F / dF;
        
        //fprintf(stderr, "  %i : h %+f  F %+20.16f  dF %+20.16f  dh %e\n", cnt, h, F, dF, h-h_old);
        
        if ( h < hmin )
        {
            h = 0.5 * ( h_old + hmin );
            continue;
        }
        
#if ( 0 )
        if ( cnt > 16 )
            fprintf(stderr, "projectEllipse fails %u :  h %+f  F %+e  dh %e\n", cnt, h, F, h-h_old);
#endif

        if ( ++cnt > 20 )
            break;
        
    } while ( h > h_old );

    // calculate the projection from h
    *pX = wX * aa / ( aa + h );
    *pY = wY * bb / ( bb + h );
    
#if ( 0 )
    // verify that projection is on ellipse:
    double F = 1 - ( pX*pX/aa + pY*pY/bb );
    fprintf(stderr, " %2i  >>> h %12.8f  F  %+e\n", cnt, h, F);
#endif
}






/**
 projectEllipsoid calculates the projection P = (pX, pY, pZ) of the point W = (wX, wY, wZ)
 on the ellipse that is aligned with the X and Y axis, and has radii (lenX, lenY, lenZ).
 
 Method:
 
 The normal to the ellipse at position P is N = ( pX / lenX^2, pY / lenY^2, pZ / lenZ^2 ),
 and we can thus write W = P + h * N, for some scalar `h' which leads to:
 @code
 pX = wX / ( 1 + h / lenX^2 );
 pY = wY / ( 1 + h / lenY^2 );
 pZ = wZ / ( 1 + h / lenZ^2 );
 @endcode
 
 Moreover, the projection should be on the ellipse and thus `h` should be a zero of:
 @code
 F(h) = ( pX / lenX )^2 + ( pY / lenY )^2 + ( pZ / lenZ )^2 - 1
 @endcode
 We follow Newton's rule to find the root of F(h), and use the formula above to
 calculate the projection.
 */
void projectEllipsoid(double  p[3],
                      const double w[3],
                      const double len[3])
{
    // handle special cases:
    if ( w[0] == 0 )
    {
        p[0] = 0;
        projectEllipse(p+1, p+2, w[1], w[2], len[1], len[2]);
        return;
    }
    if ( w[1] == 0 )
    {
        p[1] = 0;
        projectEllipse(p+0, p+2, w[0], w[2], len[0], len[2]);
        return;
    }
    if ( w[2] == 0 )
    {
        p[2] = 0;
        projectEllipse(p+0, p+1, w[0], w[1], len[0], len[1]);
        return;
    }

    double aa = len[0] * len[0];
    double bb = len[1] * len[1];
    double cc = len[2] * len[2];
    
    double h = 0, h_old;

    // we derive a lower limit for 'h' from  pX^2 + pY^2 + pZ^2 < max(lenX,lenY,lenZ)^2
    double RR = ( bb < aa ) ? ( cc < aa ? aa : cc ) : ( cc < bb ? bb : cc );
    // 'hmin' is the minimum value that 'h' can have
    double hmin = sqrt( ( w[0]*w[0]*aa*aa + w[1]*w[1]*bb*bb + w[2]*w[2]*cc*cc ) / RR ) - RR;

    // we derive another lower limit for 'h' from  |pX| < lenX
    h = ( fabs(w[0]) - len[0] ) * len[0];
    if ( h > hmin )
        hmin = h;

    // we derive another lower limit for 'h' from  |pY| < lenY
    h = ( fabs(w[1]) - len[1] ) * len[1];
    if ( h > hmin )
        hmin = h;
    
    // we derive another lower limit for 'h' from  |pZ| < lenZ
    h = ( fabs(w[2]) - len[2] ) * len[2];
    if ( h > hmin )
        hmin = h;

    if ( w[0]*w[0]/aa + w[1]*w[1]/bb + w[2]*w[2]/cc > 1 )
    {
        // if the point is outside, then 'h' should be positive:
        if ( hmin < 0 )
            hmin = 0;
    }

    h = hmin;
    //fprintf(stderr, "----- h %+f\n", h);

    /*
     Follow Newton's iteration to find the largest root.
     We start with h>0, and h should only increase
     */
    unsigned cnt = 0;
    do {
        double aah = aa + h;
        double bbh = bb + h;
        double cch = cc + h;

        double waX = w[0] / aah;
        double waY = w[1] / bbh;
        double waZ = w[2] / cch;
        
        double pXX = waX * waX * aa;
        double pYY = waY * waY * bb;
        double pZZ = waZ * waZ * cc;

        h_old = h;

        double   F = 1 - ( pXX         + pYY         + pZZ       );
        double  dF = 2 * ( pXX / aah   + pYY / bbh   + pZZ / cch );

        // Newton's method
        h -= F / dF;
        
        //fprintf(stderr, "  %i : h %+f  F %+e dh %+.20f\n", cnt, h_old, F, h-h_old);
        //fprintf(stderr, "       %+.10f   %+.10f   %+.10f   %+.10f\n", F, F/dF, ddF/dF, dddF/dF);

        if ( h < hmin )
        {
            h = 0.5 * ( h_old + hmin );
            continue;
        }

#if ( 0 )
        if ( cnt > 16 )
        {
            fprintf(stderr, "projectEllipsoid fails %u :  h %+f  F %.6e dh %.6e\n", cnt, h_old, F, h-h_old);
            //fprintf(stderr, "    pos  %+.10f     %+.10f       %+.10f\n", w[0], w[1], w[2]);
            //fprintf(stderr, "    F    %+.10f  dF %+.10f   ddF %+.10f\n", F, dF, ddF);
        }
#endif

        if ( ++cnt > 20 )
            break;
        
    } while ( h > h_old );

    // calculate the projection from h
    p[0] = w[0] * aa / ( aa + h );
    p[1] = w[1] * bb / ( bb + h );
    p[2] = w[2] * cc / ( cc + h );
    
#if ( 0 )
    // verify that projection is on ellipse
    double F = 1 - ( p[0]*p[0]/aa + p[1]*p[1]/bb + p[2]*p[2]/cc );
    fprintf(stderr, " %2i  >>> h %12.8f  F  %+e\n", cnt, h, F);
#endif
}


