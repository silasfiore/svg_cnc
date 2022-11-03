//libraries
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//header files
#include "graphic.h"
#include "svg_path.h"

int main(int argc, char **argv) {
    if (argc > 1) {
        initgraph(1200, 800);

        machine_state ms1;
        machine_state ms2;
        float vmax = 50;

        if (argc > 2) {
            float p1, p2, v1, v2;
            p1 = pathlen(argv[1]);
            p2 = pathlen(argv[2]);

            v1 = vmax * p1 / (max(p1, p2));
            v2 = vmax * p2 / (max(p1, p2));

            initialize(&ms1, v1, argv[1]);
            initialize(&ms2, v2, argv[2]);
        } else {
            initialize(&ms1, vmax, argv[1]);
            initialize(&ms2, vmax, argv[1]);
        }
        superloop(&ms1, &ms2);
        terminate();
    }

    return 0;
}