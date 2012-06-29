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
