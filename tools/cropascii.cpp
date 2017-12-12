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

// g++ -Wall -O2 -o crop cropascii.cpp 

//crops a file to a sub cube of dim x dim x dim

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "usage: " << argv[0] << " dim infile outfile";
        return -1;
    }
    int dim = atoi(argv[1]);
    fstream ifst(argv[2], fstream::in);
    fstream ofst(argv[3], fstream::out);

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
    short *imageneu = new short[dim*dim*dim];

    int pos;

    for (long z =0; z <dim;z++) {
        for (long y =0; y <dim;y++) {
            for (long x =0; x <dim;x++) {
                pos = (z*dim+y)*dim+x;
                imageneu[pos] = image[(z*sizey+y)*sizex+x];
            }
        }
    }


    ofst << dim <<" ";
    ofst << dim <<" ";
    ofst << dim << endl;
    ofst << m << " ";
    for( int i=0; i < m; i++) {
        ofst << elas[i] << " ";
    }
    ofst << endl;

    for (long z =0; z <dim;z++) {
        for (long y =0; y <dim;y++) {
            for (long x =0; x <dim;x++) {
                ofst << imageneu[(z*dim+y)*dim+x] << " ";
            }
            ofst << endl;
        }
    }


    return 0;
}
