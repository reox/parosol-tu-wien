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



#ifndef _BOUNDARY_CONDITON_H_
#define _BOUNDARY_CONDITON_H_
#include <iostream>
#include <vector>

#include "CPULayout.h"
#include "Config.h"


//! This class holds the boundary conditions.

/*! The BoundaryCondition holds the indexes and the values of the boundary nodes. 
*/
class BoundaryCondition
{

	public:
		//!Constructor
		BoundaryCondition()
		{
		}
		//!Destructor
		~BoundaryCondition()
		{
		}

		//!it holds the index of the BC to access it in the displacements array
		std::vector<t_index> FixedNodes_Ind;
		//!FixedNodes hold the value of the BC.
		std::vector<double> FixedNodes;

		//! \brief Generate BC
		/*! It reads in the boundary nodes list in and converts it to the local data structure.  
		 */
		void GenerateBC(std::vector<t_boundary_node> &list)
		{
			std::vector< t_boundary_node >::iterator it;
			for(it = list.begin(); it != list.end(); it++) {

					FixedNodes_Ind.push_back(it->first*3+it->second.dir);
					FixedNodes.push_back(it->second.disp);
			}

			return;
		}
		
		//! It prints some information of the  boundary condition. 
		friend std::ostream& operator<<(std::ostream& stream, const BoundaryCondition &bc)
		{
			return bc.print(stream);
		}

	private:
		std::ostream& print(std::ostream& stream) const;
};

std::ostream& BoundaryCondition::print(std::ostream& stream) const
{
	stream << "Boundarycondition:\n";
	stream << "   length of the fixed nodes: " << FixedNodes_Ind.size()  << "\n";
	
	return stream;
}

#endif
