//
//  min.h
//
//
//  Created by Johnathan Gross
//
//

#ifndef ____min_JLG_
#define ____min_JLG_

#include "type.h"
#include <functional>
#include <cmath>

std::array<Double,6> bracket(const std::function<Double(const Double)>&,//func
                             Double,//xa
                             Double,//xb
                             bool=false,//pos
                             Double=
                             std::numeric_limits<Double>::quiet_NaN(),//fa
                             Double=
                             std::numeric_limits<Double>::quiet_NaN());//fb
/*
 Find three points that bracket a local minimum of func with guesses at
 xa and xb
 If pos==true, then variable must be positive, otherwise can be negative
 Output: return {left,min,right,func(left),func(min),func(right)}
 where left<min<right and func(left)>func(min)<func(right)
 */

std::array<Double,6> Brentmin(const std::function<Double(const Double)>&,//func
                              Double,//a
                              Double,//b
                              Double=std::numeric_limits<Double>::min(),//prec
                              bool=false,//rel
                              bool=false,//pos
                              bool=false,//disp
                              Double=//dx
                              std::sqrt(std::numeric_limits<Double>
                                        ::epsilon()),
                              Double=//fa
                              std::numeric_limits<Double>::quiet_NaN(),
                              Double=//fb
                              std::numeric_limits<Double>::quiet_NaN());
/*
 Find the x value that minimizes func with guesses a and b to precision prec.
 If rel==true, precision is relative, otherwise precision is absolute
 If pos==true, then variable must be positive, otherwise can be negative
 Output: returns {left,min,right,func(left),func(min),func(right)}
 where left<min<right and func(left)>func(min)<func(right)
 */

#endif //____min_JLG_
