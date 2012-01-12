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

#include "Toolbox.h"

void setcoord(double a[3], double *coord) {
  setcoord(a[0],a[1],a[2], coord);
  return;
}
void setcoord(double a, double b, double c, double *coord) {
  coord[0]= 0;
  coord[1]= 0;
  coord[2]= 0;

  coord[3]= a;
  coord[4]= 0;
  coord[5]= 0;

  coord[6]= a;
  coord[7]= b;
  coord[8]= 0;

  coord[9]= 0;
  coord[10]= b;
  coord[11]= 0;

  coord[12]= 0;
  coord[13]= 0;
  coord[14]= c;

  coord[15]= a;
  coord[16]= 0;
  coord[17]= c;

  coord[18]= a;
  coord[19]= b;
  coord[20]= c;

  coord[21]= 0;
  coord[22]= b;
  coord[23]= c;

  return;
}
