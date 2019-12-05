
/**
 These two C-routines calculate the projection P of a point W on the ellipse,
 or on the surface of an Ellipsoid in 3D.
 The ellipse/ellipsoid is aligned with the principal axes (X, Y, Z),
 and has radii (radX, radY, radZ).

 Francois Nedelec. Copyright 2007-2017 EMBL.
 
 This code is Open Source covered by the GNU GPL v3.0 License
 */
 

#ifndef PROJECT_ELLIPSE_H
#define PROJECT_ELLIPSE_H


/// set the floating point precision here:
typedef double real;


/// calculate `(pX, pY)`, the projection of `(wX, wY)` on the ellipse of axes `radX, radY`
void projectEllipse(real* pX,  real* pY,
                    real wX,   real wY,
                    real radX, real radY);


/// calculate `p`, the projection of a 3D point `w` on the ellipse of axes given in `rad[]`
void projectEllipsoid(real p[3],
                      const real w[3],
                      const real rad[3]);



#endif
