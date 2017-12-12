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

// g++ -Wall -O2 -o statsofascii statsofascii.cpp

//This program counts full and empty voxel

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 1)
        cout << "usage: " << argv[0] << "infile ";
    fstream ifst(argv[1], fstream::in);
    long sizex, sizey, sizez;
    ifst >> sizex;
    ifst >> sizey;
    ifst >> sizez;
    int m;
    long mat1 = 0, mat2 =0;
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
                if (image[(z*sizey+y)*sizex+x] == 0)
                    mat1++;
                else
                    mat2++;
            }
        }
    }
    cout << "empty: :" << mat1 << " bone: " << mat2 << " ratio: " << ((double) mat2)/ (mat1+mat2);
    cout << endl;
    return 0;
}
