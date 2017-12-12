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

//g++ FractalGenerator.cpp -o FractalGenerator

//Generates a fractal with a base pattern that is read from a small image


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

struct image
{
    image(int x,int y, int z): dimx(x), dimy(y), dimz(z), _data(x*y*z)
    {
    }

    union{
        struct{ short dimx, dimy, dimz;};
        short dim[3];
    };

    vector<short> _data;

    void read(fstream &ifst) {
        for (long z =0; z <dimz;z++) {
            for (long y =0; y <dimy;y++) {
                for (long x =0; x <dimx;x++) {
                    ifst >> _data[(z*dimy+y)*dimx+x];
                }
            }
        }
    }

    inline void put(int x, int y, int z, short value)
    {
        _data[(z*dimy+y)*dimx+x]=value; 
    }

    inline short get(int x, int y, int z)
    {
        return _data[(z*dimy+y)*dimx+x]; 
    }

    void write(fstream &ofst) {
        for (long z =0; z <dimz;z++) {
            for (long y =0; y <dimy;y++) {
                for (long x =0; x <dimx;x++) {
                    ofst <<_data[(z*dimy+y)*dimx+x] << " ";
                }
                ofst << endl;
            }
        }
    }

    void putimage(image &im, int posx, int posy, int posz) {
        short v;
        for (long z =0; z <im.dimz;z++) {
            for (long y =0; y <im.dimy;y++) {
                for (long x =0; x <im.dimx;x++) {
                    v = im.get(x,y,z);
                    put(posx+x,posy+y, posz+z,v);
                }
            }
        }
    }

};

void fractal(int it, image &base, image &out, int posx, int posy, int posz)
{
    if (it < 1) {
        out.putimage(base, posx,posy,posz);
        return;
    }
    for (long z =0; z <base.dimz;z++) {
        for (long y =0; y <base.dimy;y++) {
            for (long x =0; x <base.dimx;x++) {
                if (base.get(x,y,z) != 0) {
                    fractal(it-1, base, out, (posx+x)*base.dimx, (posy+y)*base.dimy, (posz+z)*base.dimz);
                }
            }
        }
    }
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "usage: " << argv[0] << " NumberRecursion infile outfile";
        return 0;
    }
    int maxiter = atoi(argv[1]);
    fstream ifst(argv[2], fstream::in);
    fstream ofst(argv[3], fstream::out);

    int sizex, sizey, sizez;
    ifst >> sizex;
    ifst >> sizey;
    ifst >> sizez;
    int osizex = pow(sizex, maxiter);
    int osizey = pow(sizey, maxiter);
    int osizez = pow(sizez, maxiter);
    int m;
    ifst >> m;
    double *elas = new double[m];
    for( int i=0; i < m; i++) {
        ifst >> elas[i];
    }

    image inpic(sizex, sizey, sizez);
    inpic.read(ifst);

    image outpic(osizex, osizey, osizez);
    fractal(maxiter-1, inpic, outpic, 0,0,0);

    ofst << osizex <<" ";
    ofst << osizey <<" ";
    ofst << osizez << endl;
    ofst << m << " ";
    for( int i=0; i < m; i++) {
        ofst << elas[i] << " ";
    }
    ofst << endl;

    outpic.write(ofst);

    return 0;
}
