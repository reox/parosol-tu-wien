#ifndef CHECKFORBIGGESTCOMPONENT_H
#define CHECKFORBIGGESTCOMPONENT_H

#include "Image.h"
#include <set>
#include <algorithm>
#include <iostream>

struct coordinate{
	coordinate(short lx=0, short ly=0, short lz=0):x(lx), y(ly), z(lz)
	{
	}
	short x,y,z;

	bool operator<(const coordinate &rhs) const {
		if (z < rhs.z)
			return true;
		if (z == rhs.z) {
			if (y < rhs.y)
				return true;
			if ( y == rhs.y) {
				if (x < rhs.x)
					return true;
			}
		}
		return false;
	}
};


class CheckForBiggestComponent {
public:
    CheckForBiggestComponent():NumberMaterials(0), _bone(0)
    {
    }
    int NumberMaterials;
    Image *_bone;

    long FindFirstBoneElement() {
        unsigned long i;
        for(i = 0; i <_bone->_data.size(); i++)
            if (_bone->_data[i] > 0 && _bone->_data[i] < NumberMaterials)
                return i;
        return i; 
    }

    void ComputeCoordinateFromIndex(long i, short &x, short &y, short &z)
    {
        x = i % _bone->dimx;
        i = i/ _bone->dimx;
        y = i % _bone->dimy;
        i = i/ _bone->dimy;
        z = i % _bone->dimz;
        return;
    }

    bool TestAndFlagAsVisited(short x, short y, short z, short component_nr) {
        short mat_of_active = _bone->get(x, y,z);
        if (mat_of_active > 0 && mat_of_active < NumberMaterials) {
                _bone->put(x,y,z,-(component_nr*NumberMaterials+mat_of_active));
                return true;
        }
        return false;
    }

    void VisitElement(short x, short y, short z, short component_nr, long &count, std::set<coordinate> &activated) {
        if (TestAndFlagAsVisited(x, y, z,component_nr)) {
            activated.insert(coordinate(x,y,z));
            count++;
        }
    }

    void operator()(Image &bone) {
        std::vector<int> components;
        _bone = &bone;
        NumberMaterials = *max_element(bone._data.begin(), bone._data.end())+1;
        unsigned long i = FindFirstBoneElement();

        short x, y, z;
        coordinate tmp;
        long count = 0;
        short iter =0;
        while (i < _bone->_data.size()) {
            std::set<coordinate> activated;
            count =0;
            ComputeCoordinateFromIndex(i,x,y,z);
            VisitElement(x, y, z,iter, count, activated);
            while(!activated.empty()) {
                tmp = *activated.begin();
                activated.erase(activated.begin());
                x = tmp.x;
                y = tmp.y;
                z = tmp.z;
                for(int n = -1; n <2; n+=2)
                    if (x+n >= 0 && x+n < bone.dimx)
                        VisitElement(x+n, y, z,iter, count, activated);
                for(int n = -1; n <2; n+=2) {
                    if (y+n >= 0 && y+n < bone.dimy)
                        VisitElement(x, y+n, z,iter, count, activated);
                }
                for(int n = -1; n <2; n+=2)
                    if (z+n >= 0 && z+n < bone.dimz)
                        VisitElement(x, y, z+n, iter, count, activated);
            }

            components.push_back(count);
        
            i = FindFirstBoneElement();
            iter++;
        }

        int biggest_region =  max_element(components.begin(), components.end())-components.begin();
        std::cout << "Nr. of components:      " << components.size() << "\n";
        std::cout << "Biggest region:         " << biggest_region << "\n";
        std::cout << "Size of biggest region: " << components[biggest_region] << "\n";

        short tmp_mat;
        for(unsigned long i = 0; i <_bone->_data.size(); i++) {
            tmp_mat = -_bone->_data[i];

            if (tmp_mat > biggest_region * NumberMaterials &&
                tmp_mat < (biggest_region+1)*NumberMaterials)

                _bone->_data[i] = tmp_mat % NumberMaterials;
            else
                _bone->_data[i] = 0;
        }
        return;
    }
};

#endif
