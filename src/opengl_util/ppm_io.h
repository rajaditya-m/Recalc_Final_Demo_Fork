#pragma once
#include <vector>
namespace dj {
bool ReadPPM(const char *filename, std::vector<unsigned char>& image, int &width, int &height);
}

