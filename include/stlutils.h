#pragma once

#include "surface.h"
#include "volume.h"

#include <fstream>

int WriteSurfaceSTL(std::ofstream& stl_file, Surface* surface,
                    bool ascii, int res[2]);

int WriteVolumeSTL(std::ofstream& stl_file, Volume* volume,
                   bool ascii, int res[3]);
