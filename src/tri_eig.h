/* The eigensolver is taken from http://www.mymathlib.com/ on 30.10.2011
 * At that moment the code was under public domain
 */

#ifndef TRI_EIG_H
#define TRI_EIG_H

#include <eigen2/Eigen/Core>

//eigenvalues in diagonal
static double Calculate_Shift(double d2, double d1, double off) 
{
   double h;

   h = ( d2 - d1 ) / ( off + off );
   off = sqrt( h * h + 1.0 );
   d2 = ( h < 0.0 ) ? h - off : h + off;
   return d1 - off / d2;
}

int QL_Tridiagonal_Symmetric_Matrix(Eigen::VectorXd  &diagonal, Eigen::VectorXd &p_off,
                                   int n, int max_iteration_count)
{
   int i, j, k;
   int iteration = 0;
   double epsilon; 

   double s,c,g,h,q;
   double shift;
   double dum;

   p_off[n-1] = 0.0;
   for (i = 0; i < n; i++) {
      for (iteration = 0; iteration < max_iteration_count; iteration++) { 
         for (k = i; k < (n - 1); k++) {
            epsilon = DBL_EPSILON * ( fabs(diagonal[k]) + fabs(diagonal[k+1]) );
            if ( fabs(p_off[k]) <= epsilon ) break;
         }
         if ( k == i ) break;
         shift = Calculate_Shift(diagonal[i+1], diagonal[i], p_off[i]);
         q = diagonal[k] - shift;
         c = 1.0;
         s = 1.0;
         shift = 0.0;
         for (j = k - 1; j >= i; j--) {
            h = c * p_off[j];
            g = s * p_off[j];
            if ( fabs( g ) >= fabs( q ) ) {
               c = q / g;
               dum = sqrt( c * c + 1.0 );
               p_off[j+1] = g * dum;
               s = 1.0 / dum;
               c *= s;
            }
            else {
               s = g / q;
               dum = sqrt( s * s + 1.0 );
               p_off[j+1] = q * dum;
               c = 1.0 / dum;
               s *= c;
            }
            q = diagonal[j + 1] - shift;
            dum = s * (diagonal[j] - q) + 2.0 * c * h;
            shift = s * dum;
            diagonal[j+1] = q + shift;
            q = c * dum - h;
         }
         diagonal[i] -= shift;
         p_off[i] = q;
         p_off[k] = 0.0;
      }
   }
   if ( iteration >= max_iteration_count ) return -1;    
   return 0;
}
#endif /* TRI_EIG_H*/
