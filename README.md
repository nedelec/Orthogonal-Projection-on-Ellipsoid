# Orthogonal projection of a point onto an ellipse (2D) or ellipsoid (3D)

We provide two C-routines to calculate the projection P of a point W on an ellipse aligned with the principal axes.
The projection P is the point of the ellipse that is closest to W. 
The result is calculated with the precision of the machine.

----

/// set the floating point precision here:
typedef double real;


----


/// calculate `(pX, pY)`, the projection of `(wX, wY)` on the ellipse of axes `radX, radY`
void projectEllipse(real* pX,  real* pY,
                    real wX,   real wY,
                    real radX, real radY);



----


/// calculate `p`, the projection of a 3D point `w` on the ellipse of axes given in `rad[]`
void projectEllipsoid(real p[3],
                      const real w[3],
                      const real rad[3]);

