#ifndef __UTILS_H // Header guards
#define __UTILS_H

#include "GeometryLibrary.hpp"

using namespace std;

namespace GeometryLibrary {

bool importData(const string& path, Fractures& fract);

bool NearFractures(const Fractures& frc, unsigned int id_fract1, unsigned int id_fract2);

bool IntersectionFractures(Fractures& frc, unsigned int id_fract1, unsigned int id_fract2);

} //namespace

#endif
