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

//g++ -Wall -O2 -o extrudeiny extrudeiny.cpp

// extrude an 2d ascii file into a 3d asciifile

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "usage: " << argv[0] << " infile outfile";
        return 0;
    }
    fstream ifst(argv[1], fstream::in);
    fstream ofst(argv[2], fstream::out);
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
        for (long x =0; x <sizex;x++) {
            ifst >> image[z*sizex+x];
        }
    }
    short *imageneu = new short[sizex*sizey*sizez];

    int pos;

    for (long z =0; z <sizez;z++) {
        for (long x =0; x <sizex;x++) {
            for (long y =0; y <sizey;y++) {
                pos = (z*sizey+y)*sizex+x;
                imageneu[pos] = image[z*sizex+x];
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
