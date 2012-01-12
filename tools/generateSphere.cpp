/*
 * ParOSol: a parallel FE solver for trabecular bone modeling
 * Copyright (C) 2011, Cyril Flaig
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//generates a sphere in ascii format

//g++ -Wall -O2 -o Sphere generateSphere.cpp 

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
  fstream ofst("sphere.txt", fstream::out);
  int sizex, sizey, sizez;
  int s = 101;
  double cx, cy;
  cx = cy = ((double)s-1)/2;
  double cz = cx -2;
  sizex = sizey = s;
  sizez = s - 4;
  int m = 3;
  //cout << "m: " << m << endl;
  double *elas = new double[m];
  elas[0] = 0;
  elas[1] = 5.5;
  elas[2] = 1.1;

  short *imageneu = new short[sizex*sizey*sizez];

  int pos;
  double r;

  for (long z =0; z <sizez;z++) {
    for (long y =0; y <sizey;y++) {
      for (long x =0; x <sizex;x++) {
		  r = (x-cx)*(x-cx)+(y-cy)*(y-cy)+(z-cz)*(z-cz);
          pos = (z*sizey+y)*sizex+x;


		  if ( r < cx*cx) {
			  if ( r < (cx*cx)*0.64) {
				  if (r < (cx*cx)*0.25) {
                      imageneu[pos] = 0;
                  } else { 
                      imageneu[pos] = 1;
                  }
              } else {
                  imageneu[pos] = 2;
              }
          } else {
              imageneu[pos] = 0;
          }
      }
    }
  }


  ofst << sizex <<" ";
  ofst << sizey <<" ";
  ofst << sizez << endl;
  ofst << m << " ";
  //cout << "m: " << m << endl;
  for( int i=0; i < m; i++) {
    ofst << elas[i] << " ";
  }
  ofst << endl;

  for (long z =0; z <sizez;z++) {
    for (long y =0; y <sizey;y++) {
      for (long x =0; x <sizex;x++) {
        ofst << imageneu[(z*sizey+y)*sizex+x] << " ";
      }
      ofst << endl;
    }
  }


  return 0;
}
