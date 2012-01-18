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

//g++ scaleimage.cpp -o scaleimage


#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 6) {
    cout << "usage: " << argv[0] << "scalex scaley scalez infile outfile";
    return 0;
    }
  int scalex = atoi(argv[1]);
  int scaley = atoi(argv[2]);
  int scalez = atoi(argv[3]);
  fstream ifst(argv[4], fstream::in);
  fstream ofst(argv[5], fstream::out);
  int sizex, sizey, sizez;
  ifst >> sizex;
  ifst >> sizey;
  ifst >> sizez;
  int m;
  ifst >> m;
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
  short *imageneu = new short[sizex*scalex*sizey*scaley*sizez*scalez];

  int pos;
  int sz = scalez * sizez;
  int sy = scaley * sizey;
  int sx = scalex * sizex;
  int zold, yold, xold;

  for (long z =0; z <sz;z++) {
    for (long y =0; y <sy;y++) {
      for (long x =0; x <sx;x++) {
          zold = z/scalez;
          yold = y/scaley;
          xold = x/scalex;
          imageneu[(z*sy+y)*sx+x] = image[(zold*sizey+yold)*sizex+xold];
      }
    }
  }


  ofst << sizex*scalex <<" ";
  ofst << sizey*scaley <<" ";
  ofst << sizez*scalez << endl;
  ofst << m << " ";
  //cout << "m: " << m << endl;
  for( int i=0; i < m; i++) {
    ofst << elas[i] << " ";
  }
  ofst << endl;

  for (long z =0; z <sizez*scalez;z++) {
    for (long y =0; y <sizey*scaley;y++) {
      for (long x =0; x <sizex*scalex;x++) {
        ofst << imageneu[(z*sizey*scaley+y)*sizex*scalex+x] << " ";
      }
      ofst << endl;
    }
  }


  return 0;
}
