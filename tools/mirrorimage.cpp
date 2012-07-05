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

//g++ -Wall mirrorimage.cpp -o mirrorimage

//3Dmirroring of a textimage. 

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cout << "usage: " << argv[0] << "factor infile outfile";
    return 0;
    }
  int scale = atoi(argv[1]);
  fstream ifst(argv[2], fstream::in);
  fstream ofst(argv[3], fstream::out);
  int sizex, sizey, sizez;
  ifst >> sizex;
  ifst >> sizey;
  ifst >> sizez;
  int m;
  ifst >> m;
  //cout << "m: " << m << endl;
  double *elas = new double[m];
  for( int i=0; i < m; i++) {
    ifst >> elas[i];
    }

  short *image = new short[sizex*sizey*sizez];
  for (long z =0; z <sizez;z++) {
    for (long y =0; y <sizey;y++) {
      for (long x =0; x <sizex;x++) {
        ifst >> image[(z*sizey+y)*sizex+x];
      }
    }
  }
  short *imageneu = new short[sizex*sizey*sizez*scale*scale*scale];

  int pos;
  int sz = scale * sizez;
  int sy = scale * sizey;
  int sx = scale * sizex;
  int nx,ny,nz;

  for (long z =0; z <sizez;z++) {
    for (long y =0; y <sizey;y++) {
      for (long x =0; x <sizex;x++) {
        for(int i =0; i < scale; i++) {
          for(int j =0; j < scale; j++) {
            for(int k =0; k < scale; k++) {
              nx = x;
              if ( k%2==1 )
                nx = sizex - 1 -x;
              nx += sizex*k;

              ny = y;
              if ( j%2==1 )
                ny = sizey - 1 -y;
              ny += sizey*j;

              nz = z;
              if ( i%2==1 )
                nz = sizez - 1 -z;
              nz += sizez*i;

              pos = (nz*sy+ny)*sx+nx;
              imageneu[pos] = image[(z*sizey+y)*sizex+x];
            }
          }
        }
      }
    }
  }


  ofst << sizex*scale <<" ";
  ofst << sizey*scale <<" ";
  ofst << sizez*scale << endl;
  ofst << m << " ";
  //cout << "m: " << m << endl;
  for( int i=0; i < m; i++) {
    ofst << elas[i] << " ";
  }
  ofst << endl;

  for (long z =0; z <sizez*scale;z++) {
    for (long y =0; y <sizey*scale;y++) {
      for (long x =0; x <sizex*scale;x++) {
        ofst << imageneu[(z*sizey*scale+y)*sizex*scale+x] << " ";
      }
      ofst << endl;
    }
  }


  return 0;
}
