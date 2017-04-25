# Orthogonal-Projection-on-Ellipsoid
Calculate the projection of a 3D-point on the surface of an ellipsoid,
or the projection of a 2D-point on the surface of an ellipse.


Two C-routines to calculate the projection P of a point W on the ellipse,

----
/// calculate `(pX, pY)`, the projection of `(wX, wY)` on the ellipse of axes `lenX, lenY`

  void projectEllipse2D(double* pX,        double* pY,
                        const double wX,   const double wY,
                        const double lenX, const double lenY,
                        const double precision);


/// calculate `p`, the projection of a 3D point `w` on the ellipsoid of axes given in `len[]`

  void projectEllipse3D(double p[3],
                        const double w[3],
                        const double len[3],
                        const double precision);
                      
----
