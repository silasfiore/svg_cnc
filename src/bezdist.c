#include "bezdist.h"

double complex bez(double u, int n, double complex *ctrl) {
    if (n > 1) {
        return (1.0 - u) * bez(u, n - 1, ctrl) + u * bez(u, n - 1, ctrl + 1);
    } else {
        return *ctrl;
    }
}

void clip(int n, double *out, double ubot, double utop, double *ctrl) {
    static int j, k, deg;
    static double temp;
    k = deg = n - 1;
    temp = 1 - ubot;
    out[deg] = ctrl[deg];
    for (j = 0; j < k; j++) {
        out[j] = temp * ctrl[j] + ubot * ctrl[j + 1];
    }
    while (k--) {
        for (j = 0; j < k; j++) {
            out[j] = temp * out[j] + ubot * out[j + 1];
        }
    }
    utop = (utop - ubot) / temp;
    temp = 1 - utop;
    for (k = 0; k < deg; k++) {
        for (j = deg; j > k; j--) {
            out[j] = temp * out[j - 1] + utop * out[j];
        }
    }
}

void udome(double *umin, double *umax, int n, double *ctrl, int kmin, int kmax,
           unsigned int mask) {
    static double delta, d, dmax;
    static int k, kdmax;
    static unsigned long int flag;
    delta = ctrl[kmin] - ctrl[kmax];
    dmax = 0;
    kdmax = kmax;
    flag = 1;
    for (k = 0; k < n; k++) {
        if (mask & flag) {
            d = delta * (k - kmin) + (kmax - kmin) * (ctrl[k] - ctrl[kmin]);
            if (d < 0.0) {
                mask |= flag;
                mask ^= flag;
            } else if (d > dmax) {
                dmax = d;
                kdmax = k;
            }
        }
        flag <<= 1;
    }
    if (dmax > 0) {
        udome(umin, umax, n, ctrl, kmin, kdmax, mask);
        udome(umin, umax, n, ctrl, kdmax, kmax, mask);
    } else {
        d = ctrl[kmin] / delta;
        if (d <= 1 && d >= 0) {
            delta = ((1 - d) * kmin + d * kmax) / (n - 1);
            if (delta < *umin) {
                *umin = delta;
            }
            if (delta > *umax) {
                *umax = delta;
            }
        }
    }
}

void ppa(int n, int *ridx, double *roots, double *work, double *ctrl,
         double ubot, double utop) {
    static double umin, umid, umax, temp;
    static int i, v;
    if ((utop - ubot) <= TOL) {
        roots[(*ridx)++] = (utop + ubot) / 2;
    } else {
        clip(n, work, ubot, utop, ctrl);
        umin = 1;
        umax = 0;
        v = 0;
        for (i = 1; i < n; i++) {
            if (signbit(work[i - 1]) != signbit(work[i])) {
                v++;
            };
        }
        if (v) {
            udome(&umin, &umax, n, work, 0, n - 1, pow(2, n) - 1);
            udome(&umin, &umax, n, work, n - 1, 0, pow(2, n) - 1);
            if (umin <= umax) {
                temp = utop - ubot;
                ubot += umin * temp;
                utop += (umax - 1) * temp;
                if (v == 1) {
                    ppa(n, ridx, roots, work, ctrl, ubot, utop);
                } else {
                    umid = (utop + ubot) / 2;
                    ppa(n, ridx, roots, work, ctrl, ubot, umid);
                    ppa(n, ridx, roots, work, ctrl, umid, utop);
                }
            }
        }
    }
}

void rfcn(int n, double complex q, double complex *ctrl, double *rootfcn) {
    int i, j;
    double fi, fj;
    double complex a, b;
    for (i = 0; i < (2 * n - 2); i++) {
        rootfcn[i] = 0.0;
    }

    fi = 1.0;
    fj = 1.0;
    for (i = 0; i < (n - 1); i++) {
        if (i > 0) {
            fi = fi * (i - n + 1) / (i - 2 * n + 2);
            fj = fi;
        }
        for (j = 0; j < n; j++) {
            if (j > 0) {
                fj = fj * ((i + j) * (j - n)) / (j * (i + j - 2 * n + 2));
            }
            a = ctrl[i + 1] - ctrl[i];
            b = q - ctrl[j];
            rootfcn[i + j] += creal(a) * creal(b) * fj;
            rootfcn[i + j] += cimag(a) * cimag(b) * fj;
        }
    }
}

void lsqbezmat(int m, int n, double *b, double *u) {
    // m >= n
    int i, j, k;
    double compu, temp1, temp2;

    for (i = 0; i < m * n; i++) {
        b[i] = 0;
    }
    for (i = 0; i < m; i++) {
        b[i] = 1.0;
        compu = (1 - *u);
        for (j = 0; j < (n - 1); j++) {
            temp1 = 0.0;
            for (k = 0; k < (j + 2); k++) {
                temp2 = b[i + k * m] * *u;
                b[i + k * m] *= compu;
                b[i + k * m] += temp1;
                temp1 = temp2;
            }
        }
    }
}

double bezdist(double complex q, int n, double complex *ctrl) {
    int ridx;
    double roots[2 * n - 3];
    double work[2 * n - 2];
    double rootfcn[2 * n - 2];

    double dist, d, u;
    int j;

    rfcn(n, q, ctrl, rootfcn);
    for (j = 0; j < (2 * n - 3); j++) {
        roots[j] = -1.0;
    }
    ridx = 0;
    ppa((2 * n - 2), &ridx, roots, work, rootfcn, 0, 1);
    u = 0.0;
    d = cabs(ctrl[0] - q);

    dist = cabs(ctrl[n - 1] - q);
    if (dist < d) {
        u = 1.0;
        d = dist;
    }

    for (j = 0; j < 2 * n - 3; j++) {
        if (roots[j] > 0.0 && roots[j] < 1.0) {
            dist = cabs(bez(roots[j], n, ctrl) - q);
            if (dist < d) {
                u = roots[j];
                d = dist;
            }
        }
    }
    return d;
}