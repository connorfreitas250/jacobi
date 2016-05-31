#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  double original;
  double new;

  FILE* originalfile = fopen(argv[1], "r");
  FILE* newfile = fopen(argv[2], "r");

  for (int i = 0; i < 2048; ++i) {
    for (int j = 0; j < 2048; ++j) {
      fscanf(originalfile, "%lf", &original);
      fscanf(newfile, "%lf", &new);
      if (fabsf(original - new) > 0.0001) {
        printf("error at i %d j %d\n", i, j);
        exit(1);
      }
    }
  }
  fclose(originalfile);
  fclose(newfile);
}