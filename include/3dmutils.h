#pragma once

#include "curve.h"
#include "surface.h"
#include <opennurbs.h>

Curve* ONCurveToGoCurve(const ON_Curve* curve);
Surface* ONSurfaceToGoSurface(const ON_Surface* surf);
