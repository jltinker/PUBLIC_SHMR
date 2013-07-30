#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

void deconvolved_smf(void);
void bias_from_simulation(void);

void test(int argc, char **argv)
{
  shambias_loop();
  bias_from_simulation();
  deconvolved_smf();

}
