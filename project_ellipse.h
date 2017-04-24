
/**
 These two C-routines calculate the projection P of a point W on the ellipse,
 or on the surface of an Ellipsoid in 3D.
 The ellipse/ellipsoid is aligned with the principal axes (X, Y, Z),
 and has radii (lenX, lenY, lenZ).

 Francois Nedelec. Copyright 2007-2017 EMBL.
 
 This code is Open Source covered by the GNU GPL v3.0 License
 */
 

#ifndef PROJECT_ELLIPSE_H
#define PROJECT_ELLIPSE_H

/// calculate `(pX, pY)`, the projection of `(wX, wY)` on the ellipse of axes `lenX, lenY`
void projectEllipse2D(double* pX,        double* pY,
                      const double wX,   const double wY,
                      const double lenX, const double lenY,
                      const double precision);


/// calculate `p`, the projection of a 3D point `w` on the ellipse of axes given in `len[]`
void projectEllipse3D(double p[3],
                      const double w[3],
                      const double len[3],
                      const double precision);



#endif
