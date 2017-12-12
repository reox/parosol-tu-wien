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

//g++ -o perlin -O3 -lnoise perlin.cpp

// perlin generates an artificial bone from pseudo random number with
// perlin noise. 
// Most important options:
// --seed changes the seed
// --lowerthres and --upperthres should be symmetric around 0. With a small
// band you get a sparse bone.
// --fsize set the feature size. A big value means a zoom in. A small value
// a zoom out.
// --size[z,y,x] image size 
// To generate the same image with different resolution fsize/size have to be
// constant.


#include "Image.h"
#include "CheckForBiggestComponent.h"
#include <noise/noise.h>
#include <fstream>
#include <iostream>
#include <string>

int main (int argc, char** argv)
{
    std::string usage = std::string("usage: \n") + argv[0]
        + std::string(" [--octaves arg (1)]\n")
        + std::string(" [--seed arg (0)]\n")
        + std::string(" [--upperthres arg (0.23)]\n")
        + std::string(" [--lowerthres arg (-0.23)]\n")
        + std::string(" [--fsize arg (59)]\n")
        + std::string(" [--sizex arg (160)]\n")
        + std::string(" [--sizey arg (160)]\n")
        + std::string(" [--sizez arg (160)]\n")
        + std::string(" filename\n");

    if (argc < 2) {
        std::cout << usage;
        return -1;
    }

    int freq = 1;
    int octaves = 1;
    int seed = 0;
    int sizex = 160;
    int sizey = 160;
    int sizez = 160;
    float upperthres = 0.255;
    float lowerthres = -0.255;
    float fsize=59;
    std::string file;

    int i = 1;
    std::string param;
    while(i < argc) {
        param = std::string(argv[i]);
        if (param.compare("--freq") == 0) {
            i++;
            freq = atoi(argv[i]);
        } else if(param.compare("--octaves") == 0) {
            octaves = atoi(argv[++i]);
        } else if(param.compare("--seed") == 0) {
            seed = atoi(argv[++i]);
        } else if(param.compare("--sizex") == 0) {
            sizex = atoi(argv[++i]);
        } else if(param.compare("--sizey") == 0) {
            sizey = atoi(argv[++i]);
        } else if(param.compare("--sizez") == 0) {
            sizez = atoi(argv[++i]);
        } else if(param.compare("--upperthres") == 0) {
            upperthres = atof(argv[++i]);
        } else if(param.compare("--lowerthres") == 0) {
            lowerthres = atof(argv[++i]);
        } else if(param.compare("--fsize") == 0) {
            fsize = atof(argv[++i]);
        } else {
            file = param;
        }
        i++;
    }

    std::cout << "used parameter:\n"
        << " octaves:       \t" << octaves << "\n"
        << " seed:          \t" << seed << "\n"
        << " upperthreshold:\t" << upperthres << "\n"
        << " lowerthreshold:\t" << lowerthres << "\n"
        << " fsize:         \t" << fsize<< "\n";

    std::fstream ofst(file.c_str(), std::fstream::out);
    noise::module::Perlin perlin;
    perlin.SetNoiseQuality(noise::QUALITY_BEST);
    perlin.SetFrequency(freq);
    perlin.SetOctaveCount(octaves);
    perlin.SetSeed(seed);

    int out;
    Image bone(sizex, sizey, sizez);

    for(int z=0; z<sizez;z++) 
        for(int y=0; y<sizey;y++) 
            for(int x=0; x<sizex;x++) {
                //Take two Perlin random numbers and use the intersection
                double value = perlin.GetValue (x/fsize+0.5, y/fsize+0.5, z/fsize+0.5);
                double value2 = perlin.GetValue (x/fsize*0.95+33.3, y/fsize*0.99+42.42, z/fsize*0.9+314.314);
                out=0;
                if (value > lowerthres && value < upperthres)
                    if (value2 > lowerthres && value2 < upperthres)
                        out =1;
                bone.put(x,y,z,out);
            }

    CheckForBiggestComponent comp;
    comp(bone);

    {
        int i;
        for(i=0; i<sizey*sizex;i++)
            if (bone._data[i] != 0)
                break;
        if (i==sizey*sizex)
            std::cout << "Warning no bone at z=0\n";

        for(i=sizex*sizey*(sizez-1); i<sizez*sizey*sizex;i++)
            if (bone._data[i] != 0)
                break;
        if (i==sizez*sizey*sizex)
            std::cout << "Warning no bone at z=max\n";
    }


    ofst << sizex<<" ";
    ofst << sizey <<" ";
    ofst << sizez << std::endl;
    int m = 2;
    float elas[] = {0.0, 1.0};
    ofst << m << " ";
    for( int i=0; i < m; i++) {
        ofst << elas[i] << " ";
    }
    ofst << std::endl;
    bone.write(ofst);
    return 0;
}
