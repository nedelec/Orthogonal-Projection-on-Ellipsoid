# Orthogonal projection of a point on an ellipse (2D) or ellipsoid (3D)

We provide two C-routines to calculate the projection P of a point W on an ellipse aligned with the principal axes.
The projection P is the point of the ellipse that is closest to W. 

----
/// projectEllipse2D calculates `(pX, pY)`, the projection of `(wX, wY)` on the ellipse of axes `lenX, lenY`

  void projectEllipse2D(double* pX,        double* pY,
                        const double wX,   const double wY,
                        const double lenX, const double lenY,
                        const double precision);

----

/// projectEllipse3D calculates `p`, the projection of a 3D point `w` on the ellipsoid of axes given in `len[]`

  void projectEllipse3D(double p[3],
                        const double w[3],
                        const double len[3],
                        const double precision);
                      
