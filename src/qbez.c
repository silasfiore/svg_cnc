#include "qbez.h"
static void split(float complex *P, float complex *Q, float complex *R, int n);
static void splitcb(int nsplits, int *iq, float complex *q, float complex *c);

int ipow(int base, int exp) {
    int result = 1;
    for (;;) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

float arclenqb(float complex *q) {
    float complex a, b;
    float arclen, na, nb, d, e, f, g, sqrt1, sqrt2;

    a = q[0] - 2 * q[1] + q[2];
    b = (q[1] - q[0]);

    na = cabsf(a);
    nb = cabsf(b);

    e = nb * nb / na / na;
    d = (crealf(a) * crealf(b) + cimagf(a) * cimagf(b)) / na / na;

    f = e - d * d;

    g = 1.0f + d;

    if (na < 0.0001) {
        arclen = 2 * nb;

    } else {
        if (f < 0.0001) {
            arclen = na * (sgn(g) * g * g - sgn(d) * d * d);

        } else {
            sqrt1 = sqrtf(g * g + f);
            sqrt2 = sqrtf(e);
            arclen = na * (sqrt1 + d * (sqrt1 - sqrt2) + f * logf((g + sqrt1) / (d + sqrt2)));
        }
    }

    return arclen;
}

void evalqb(float complex *pos, float *u, float s, float complex *q) {
    int iter;
    float complex a, b;
    float na, nb, d, e, f, dg, sqrt1, sqrt2, tmp1, tmp2;
    static float g;

    a = q[0] - 2 * q[1] + q[2];
    b = (q[1] - q[0]);

    na = cabsf(a);
    nb = cabsf(b);

    e = nb * nb / na / na;
    d = (crealf(a) * crealf(b) + cimagf(a) * cimagf(b)) / na / na;

    f = e - d * d;

    if (na < 0.0001) {
        *u = s / (2 * nb);
    } else {
        if (f < 0.0001) {
            if (d > 0.0f) {
                *u = sqrtf(d * d + s / na) - d;

            } else {
                if ((s / na) < (d * d)) {
                    *u = -d - sqrtf(d * d - s / na);
                } else {
                    *u = sqrtf(s / na - d * d) - d;
                }
            }
        } else {
            sqrt2 = sqrtf(e);
            tmp1 = d + sqrt2;
            tmp2 = s / na + d * sqrt2;

            g = *u + d;

            iter = 0;
            do {
                sqrt1 = sqrtf(g * g + f);

                dg = 0.5 * g + (f * logf((g + sqrt1) / tmp1) - tmp2) / (2 * sqrt1);

                g -= dg;
                iter += 1;
            } while (fabsf(dg) > FLT_EPSILON && iter < EVALQBMAXITER);
            *u = g - d;
        }
    }

    tmp1 = 1.0f - *u;

    a = tmp1 * q[0] + *u * q[1];
    b = tmp1 * q[1] + *u * q[2];

    *pos = tmp1 * a + *u * b;

    return;
}

float arclencb(int *pidx, float complex *qctrl, float complex *cctrl) {
    float arclen = 0.0f;
    int i;

    splitcb(1, pidx, qctrl, cctrl);
    for (i = 0; i < (*pidx) / 2; i++) {
        arclen += arclenqb(qctrl + 2 * i);
    }
    return arclen;
}

void evalcb(float complex *pos, float *u, float s, int iend, float complex *q) {
    int k;
    float sk;

    for (k = 0; k < iend / 2; k++) {
        sk = arclenqb(q + k * 2);
        if (s > sk) {
            s -= sk;
        } else {
            break;
        }
    }

    evalqb(pos, u, s, q + k * 2);
}

static void split(float complex *P, float complex *Q, float complex *R, int n) {
    int i, j, k;
    float complex qtemp1, qtemp2;

    for (i = 0; i < n; i++) {
        Q[i] = P[i];
        R[i] = 0.0f;
    }
    R[n - 1] = Q[n - 1];
    for (j = 1; j < n; j++) {
        qtemp2 = Q[j - 1];
        for (k = j; k < n; k++) {
            qtemp1 = qtemp2;
            qtemp2 = (Q[k - 1] + Q[k]) / 2;
            Q[k - 1] = qtemp1;
        }
        Q[n - 1] = qtemp2;
        R[n - j - 1] = qtemp2;
    }
    return;
}

// TODO: change this function such that cctrl doesnt have to be a malloc array
// TODO: change iend to counter

// use the same array for cctrl and qctrl!!!!

static void splitcb(int n, int *pidx, float complex *qctrl, float complex *cctrl) {
    float complex s1, s2;
    double complex w[3];
    double d1, d2;

    w[0] = cctrl[0];
    w[1] = (3 * cctrl[1] - cctrl[0] + 3 * cctrl[2] - cctrl[3]) / 4;  //intersect(cctrl[0], cctrl[1], cctrl[2], cctrl[3]);
    w[2] = cctrl[3];

    s1 = (8 * cctrl[0] + cctrl[3] + 12 * cctrl[1] + 6 * cctrl[2]) / 27;
    s2 = (cctrl[0] + 8 * cctrl[3] + 6 * cctrl[1] + 12 * cctrl[2]) / 27;

    d1 = bezdist(s1, 3, w);
    d2 = bezdist(s2, 3, w);

    if (max(d1, d2) > SPLITCBERROR && (1 + ipow(2, n + 1)) <= SPLITCBNMAX) {
        float complex *botctrl = (float complex *)malloc(4 * sizeof(float complex));
        float complex *topctrl = (float complex *)malloc(4 * sizeof(float complex));
        split(cctrl, botctrl, topctrl, 4);
        free(cctrl);
        splitcb((n + 1), pidx, qctrl, botctrl);
        splitcb((n + 1), pidx, qctrl, topctrl);
        return;
    } else {
        //float s = quad_bez_len(w[0], w[1], w[2]);
        qctrl[++(*pidx)] = w[1];
        qctrl[++(*pidx)] = w[2];
        free(cctrl);
        return;
    }
}

float angle(float ux, float uy, float vx, float vy) {
    float sign;
    float tmp;

    tmp = ux * vy - uy * vx;

    sign = (float)sgn(tmp);

    tmp = (ux * vx + uy * vy) / (sqrtf(ux * ux + uy * uy) * sqrtf(vx * vx + vy * vy));

    return sign * acosf(tmp);
}

float complex intersect(float complex q0, float complex d0, float complex q2, float complex d2) {
    float x0, y0, x1, y1, x2, y2, x3, y3, x, y, d;
    x0 = creal(q0);
    y0 = cimag(q0);
    x1 = creal(q0 + d0);
    y1 = cimag(q0 + d0);
    x2 = creal(q2 - d2);
    y2 = cimag(q2 - d2);
    x3 = creal(q2);
    y3 = cimag(q2);

    d = (x0 - x1) * (y2 - y3) - (y0 - y1) * (x2 - x3);

    if (d == 0.0f) {
        x = (x0 + x3) / 2;
        y = (y0 + y3) / 2;
    } else {
        x = ((x0 * y1 - y0 * x1) * (x2 - x3) - (x0 - x1) * (x2 * y3 - y2 * x3)) / d;
        y = ((x0 * y1 - y0 * x1) * (y2 - y3) - (y0 - y1) * (x2 * y3 - y2 * x3)) / d;
    }
    return CMPLXF(x, y);
}

float arclenel(int *pidx, float complex *qctrl, float complex start, float complex end, float a, float b, float phi, int fa, int fs) {
    float complex tmp, xd, cd, c, dlast, dnext;
    float ftmp, th1, th2, ps1, ps2, dth, dps, arclen;

    xd = cexpf(I * phi) * (start - end) / 2;

    tmp = CMPLXF(a * cimag(xd) / b, -b * creal(xd) / a);

    ftmp = sqrtf((a * a * b * b - a * a * cimag(xd) * cimag(xd) - b * b * creal(xd) * creal(xd)) / (a * a * cimag(xd) * cimag(xd) + b * b * creal(xd) * creal(xd)));

    if (fa == fs) {
        cd = -ftmp * tmp;
    } else {
        cd = ftmp * tmp;
    }

    c = cexpf(-I * phi) * cd + (start + end) / 2;

    tmp = xd - cd;
    th1 = atan2f(cimagf(tmp) / b, crealf(tmp) / a);
    tmp = -xd - cd;
    th2 = atan2f(cimagf(tmp) / b, crealf(tmp) / a);

    //atan2(sin(theta1)*test.a,cos(theta1)*test.b)
    ps1 = atan2f(sinf(th1) * a, cosf(th1) * b);
    ps2 = atan2f(sinf(th2) * a, cosf(th2) * b);

    dps = ps2 - ps1;

    dth = th2 - th1;

    if ((fs == 0) && (dth > 0)) {
        dth -= 2 * M_PI;
    }
    if ((fs == 1) && (dth < 0)) {
        dth += 2 * M_PI;
    }
    if ((fs == 0) && (dps > 0)) {
        dps -= 2 * M_PI;
    }
    if ((fs == 1) && (dps < 0)) {
        dps += 2 * M_PI;
    }

    //theta = atan2(test.b*sin(psi+psi1),test.a*cos(psi+psi1))

    ftmp = dps / 8;
    //ftmp = dth / 1000;

    qctrl[0] = cexpf(-I * phi) * CMPLXF(a * cosf(th1), b * sinf(th1)) + c;
    dlast = cexpf(-I * phi) * CMPLXF(-a * sinf(th1), b * cosf(th1));

    for (int i = 1; i < 9; i++) {
        float theta;

        theta = atan2f(b * sinf(ps1 + ftmp * i), a * cosf(ps1 + ftmp * i));

        qctrl[i * 2] = cexpf(-I * phi) * CMPLXF(a * cosf(theta), b * sinf(theta)) + c;
        dnext = cexpf(-I * phi) * CMPLXF(-a * sinf(theta), b * cosf(theta));

        qctrl[i * 2 - 1] = intersect(qctrl[(i - 1) * 2], dlast, qctrl[i * 2], dnext);
        dlast = dnext;
    }

    *pidx = 16;
    arclen = 0.0f;
    for (int i = 0; i < (*pidx) / 2; i++) {
        arclen += arclenqb(qctrl + 2 * i);
    }
    return arclen;
}
