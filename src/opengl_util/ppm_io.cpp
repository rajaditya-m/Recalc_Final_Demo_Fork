#include <stdio.h>
#include "ppm_io.h"

namespace dj {
bool ReadPPM(const char *filename, std::vector<unsigned char>& image, int &width, int &height)
{
  FILE *fp = fopen(filename, "rb");
  if (!fp) return false;

  /* output header */
  float endian;
  char type[128];

  int ret = fscanf(fp, "%s\n", type);
  //To remove the comment made by Irfanview.
  //This is for Irfanview only.
  while (!feof(fp)) {
    char c = fgetc(fp);
    if (c == '\r' || c == '\n')	break;
  }

  ret = fscanf( fp, "%d %d\n", &width, &height);
  ret = fscanf( fp, "%f\n", &endian);
  ret = ret;

  image.resize(width * height * 3);
  if (strcmp(type, "P6") == 0) {
    for (int j = height - 1; j >= 0; j--)
      for (int i = 0; i < width; i++) {
        if (feof(fp)) printf("end of file already.\n");
        image[(j * width + i) * 3 + 0] = fgetc(fp);
        image[(j * width + i) * 3 + 1] = fgetc(fp);
        image[(j * width + i) * 3 + 2] = fgetc(fp);
      }
  }
  fclose( fp );
  return fp != NULL;
}

} // namespace dj
