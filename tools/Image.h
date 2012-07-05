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

#ifndef IMAGE_H
#define IMAGE_H
#include <istream>
#include <ostream>
#include <vector>

struct Image
{
    Image(int x,int y, int z): dimx(x), dimy(y), dimz(z), _data(x*y*z)
    {
    }

    union{
        struct{ short dimx, dimy, dimz;};
        short dim[3];
    };

    std::vector<short> _data;

    void read(std::istream &ist) {
        for (long z =0; z <dimz;z++) {
            for (long y =0; y <dimy;y++) {
                for (long x =0; x <dimx;x++) {
                    ist >> _data[(z*dimy+y)*dimx+x];
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

    void write(std::ostream &ost) {
        for (long z =0; z <dimz;z++) {
            for (long y =0; y <dimy;y++) {
                for (long x =0; x <dimx;x++) {
                    ost <<_data[(z*dimy+y)*dimx+x] << " ";
                }
                ost << std::endl;
            }
        }
    }

    void putimage(Image &im, int posx, int posy, int posz) {
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

#endif
