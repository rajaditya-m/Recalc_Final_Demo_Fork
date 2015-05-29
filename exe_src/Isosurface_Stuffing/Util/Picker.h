#ifndef PICKER_H
#define PICKER_H

#include "vector.h"

void pickFromXYPlane(Vector result, int x, int y, double depth = 0.0);

void convertToScreen(double& x, double& y, Vector object_coords);

#endif