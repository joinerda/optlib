//  OptLib version 1.0
//    Copyright (C) 2012  David A. Joiner
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "ol_transform.h"
#include <math.h>
#define PI 3.141592653589793

double ab_transform(double x,double a, double b) {
    return tan(-PI/2+(x-a)/(b-a)*PI);
}
double ab_invtransform(double x, double a, double b) {
    return a + (atan(x)+PI/2)/PI*(b-a);
}
double log_ab_transform(double x,double a, double b) {
    return ab_transform(log(x),log(a),log(b));
}
double log_ab_invtransform(double x, double a, double b) {
    return exp(ab_invtransform(x,log(a),log(b)));
}

