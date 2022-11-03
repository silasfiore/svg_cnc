#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bezdist.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define sgn(a) ((a < 0) ? -1 : 1)
#define round(x) ((x < 0) ? (x - 0.5) : (x + 0.5));

#define SPLITCBNMAX 17
#define SPLITCBERROR 0.0001f  //0.001f
#define EVALQBMAXITER 5

typedef struct elliptic_arc {
    int fs;
    float a;
    float b;
    float th1;
    float th2;
    float complex c;
    float complex rot;
} elliptic_arc;

float arclenqb(float complex *q);
void evalqb(float complex *pos, float *u, float s, float complex *q);

float arclencb(int *pidx, float complex *qctrl, float complex *cctrl);
void evalcb(float complex *pos, float *u, float s, int iend, float complex *q);

float arclenel(int *pidx, float complex *qctrl, float complex start, float complex end, float rx, float ry, float phi, int fa, int fs);
void evalel(float complex *pos, float *u, float s, int iend, float complex *q, float arclen);
/*
float arclencb(int *iq, float complex *q, float complex *c);
void evalcb(float complex *pos, float *u, float s, int iend, float complex *q);
*/