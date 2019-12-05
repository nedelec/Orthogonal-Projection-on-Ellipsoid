
/**
 These two C-routines calculate the projection P of a point W on the ellipse,
 or on the surface of an Ellipsoid in 3D.
 The ellipse/ellipsoid is aligned with the principal axes (X, Y, Z),
 and has radii (radX, radY, radZ).
 
 Francois Nedelec. Copyright 2007-2017 EMBL. 2019 Cambridge University
 
 This code is Open Source covered by the GNU GPL v3.0 License

 For mathematical background, check
 ALGORITHMS OF PROJECTION OF A POINT ONTO AN ELLIPSOID
 Yu. N. Kiseliov
 Lithuanian Mathematical Journal, Vol. 34, No. 2, 1994
 */

#include "project_ellipse.h"
#include <assert.h>
#include <math.h>

/**
 Calculate the projection P = (pX, pY) of the point W = (wX, wY) on the ellipse that
 is aligned with the X and Y axis, and has radii (radX, radY).
 
 Method:
 
 A vector orthogonal to the ellipse at position (X, Y) is
 
     N = ( X / radX^2, Y / radY^2 ),
 
 and we can thus write W = P + h * N, leading to:

     pX = wX * radX^2 / ( radX^2 + h );
     pY = wY * radY^2 / ( radY^2 + h );

 if wX and wY are not both null.
 
 Moreover, the projection should be on the ellipse and thus `h` should be a zero of:
 
     F(h) = ( pX / radX )^2 + ( pY / radY )^2 - 1
 
 We follow Newton's rule to find the root of F(h), and use the formula above to
 calculate the projection.
 */
void projectEllipse(real&   pX, real&  pY,
                    real    wX, real   wY,
                    real  radX, real radY)
{
    // handle the pathological cases:
    if ( wX == 0 )
    {
        pX = 0;
        pY = copysign(radY, wY);
        return;
    }
    if ( wY == 0 )
    {
        pX = copysign(radX, wX);
        pY = 0;
        return;
    }
    
    real aa = radX * radX;
    real bb = radY * radY;
    
    // we derive a lower limit for 'h' from  pX^2 + pY^2 > max(radX,radY)^2
    real RR = fmax(aa, bb);
    // 'hmin' is the minimum value that 'h' could have
    real hmin = sqrt( ( wX*wX*aa*aa + wY*wY*bb*bb ) / RR ) - RR;
    
    // we derive another lower limit for 'h' from  |pX| < radX
    hmin = fmax(hmin, ( fabs(wX) - radX ) * radX);

    // we derive another lower limit for 'h' from  |pY| < radY
    hmin = fmax(hmin, ( fabs(wY) - radY ) * radY);

    // if the point is outside, then 'h' should be positive:
    if ( wX*wX/aa + wY*wY/bb > 1  &&  hmin < 0 )
        hmin = 0;
    
    real h_old, h = hmin;

    //fprintf(stderr, " <<< %+.10f  %+.10f    hmin %+10.4f", wX, wY, hmin);

    // follow Newton's iteration to find the root
    unsigned cnt = 0;
    do {
        real aah = aa + h;
        real bbh = bb + h;

        real waX = wX / aah;
        real waY = wY / bbh;
        
        real pXX = waX * waX * aa;
        real pYY = waY * waY * bb;

        h_old = h;
        
        real F    = 1 - ( pXX         + pYY       );
        real dF   = 2 * ( pXX / aah   + pYY / bbh );

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
    pX = wX * aa / ( aa + h );
    pY = wY * bb / ( bb + h );
    
#if ( 0 )
    // verify that projection is on ellipse:
    real F = 1 - ( pX*pX/aa + pY*pY/bb );
    fprintf(stderr, " %2i  >>> h %12.8f  F  %+e\n", cnt, h, F);
#endif
}


/**
 Calculates the projection P = (pX, pY, pZ) of the point W = (wX, wY, wZ) on the ellipse that
 is aligned with the X and Y axis, and has radii (radX, radY, radZ).
 
 Method:
 
 A vector orthogonal to the ellipse at position ( X, Y, Z ) is
 
    N = ( X / radX^2, Y / radY^2, Z / radZ^2 ),
 
 and we can thus write W = P + h * N, for some scalar `h` leading to:

    pX = wX / ( 1 + h / radX^2 );
    pY = wY / ( 1 + h / radY^2 );
    pZ = wZ / ( 1 + h / radZ^2 );
 
 Moreover, the projection should be on the ellipse and thus `h` should be a zero of:

     F(h) = ( pX / radX )^2 + ( pY / radY )^2 + ( pZ / radZ )^2 - 1
 
 We follow Newton's rule to find the root of F(h), and use the formula above to
 calculate the projection.
 */
void projectEllipsoid(real  p[3],
                      const real w[3],
                      const real rad[3])
{
    assert( rad[0]==rad[0] && rad[0] > 0 );
    assert( rad[1]==rad[1] && rad[1] > 0 );
    assert( rad[2]==rad[2] && rad[2] > 0 );
    
    // handle the pathological cases:
    if ( w[0] == 0 )
    {
        p[0] = 0;
        projectEllipse(p[1], p[2], w[1], w[2], rad[1], rad[2]);
        return;
    }
    if ( w[1] == 0 )
    {
        p[1] = 0;
        projectEllipse(p[0], p[2], w[0], w[2], rad[0], rad[2]);
        return;
    }
    if ( w[2] == 0 )
    {
        p[2] = 0;
        projectEllipse(p[0], p[1], w[0], w[1], rad[0], rad[1]);
        return;
    }

    real aa = rad[0] * rad[0];
    real bb = rad[1] * rad[1];
    real cc = rad[2] * rad[2];

    // we derive a lower limit for 'h' from  pX^2 + pY^2 + pZ^2 < max(radX,radY,radZ)^2
    real RR = fmax(aa, fmax(bb,cc));
    // 'hmin' is the minimum value that 'h' can have
    real hmin = sqrt( ( w[0]*w[0]*aa*aa + w[1]*w[1]*bb*bb + w[2]*w[2]*cc*cc ) / RR ) - RR;

    // we derive another lower limit for 'h' from  |pX| < radX
    hmin = fmax(hmin, ( fabs(w[0]) - rad[0] ) * rad[0]);

    // we derive another lower limit for 'h' from  |pY| < radY
    hmin = fmax(hmin, ( fabs(w[1]) - rad[1] ) * rad[1]);
    
    // we derive another lower limit for 'h' from  |pZ| < radZ
    hmin = fmax(hmin, ( fabs(w[2]) - rad[2] ) * rad[2]);

    if ( w[0]*w[0]/aa + w[1]*w[1]/bb + w[2]*w[2]/cc > 1  &&  hmin < 0 )
    {
        // if the point is outside, then 'h' should be positive:
        hmin = 0;
    }

    real h_old, h = hmin;
    //fprintf(stderr, "----- h %+f\n", h);

    /*
     Follow Newton's iteration to find the largest root.
     We start with h>0, and h should only increase
     */
    unsigned cnt = 0;
    do {
        real aah = aa + h;
        real bbh = bb + h;
        real cch = cc + h;

        real waX = w[0] / aah;
        real waY = w[1] / bbh;
        real waZ = w[2] / cch;
        
        real pXX = waX * waX * aa;
        real pYY = waY * waY * bb;
        real pZZ = waZ * waZ * cc;

        h_old = h;

        real   F = 1 - ( pXX         + pYY         + pZZ       );
        real  dF = 2 * ( pXX / aah   + pYY / bbh   + pZZ / cch );

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
    real F = 1 - ( p[0]*p[0]/aa + p[1]*p[1]/bb + p[2]*p[2]/cc );
    fprintf(stderr, " %2i  >>> h %12.8f  F  %+e\n", cnt, h, F);
#endif
}


