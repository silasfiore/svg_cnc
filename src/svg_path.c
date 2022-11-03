//standard libraries
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//header files
#include "bezdist.h"
#include "graphic.h"
#include "qbez.h"
#include "svg_path.h"

void visualize(machine_state *ms, unsigned long color) {
    drawpixel((int)(creal(ms->q) * 4), (int)(cimag(ms->q) * 4), color);
}

void high_lvl_ctrl(machine_state *ms) {
    float t;  // t_star, t_end;

    t = (float)(clock() - ms->ts) / CLOCKS_PER_SEC;
    //t_end = ms->curr.arclen / ms->vel;

    // timing law
    //TODO: test timing 20.410714,71.726189 63.5,7.559524 54.428566,6.803572 27.59227,3.590774 -8.31548,25.513391  law once more modes are implemented
    // STOP and GO must be properties of the segments

    /*
    if (t < t_end) {
        tmp = ms->acc;
        if (cabsf(ms->next.tstart - ms->curr.tend) < 0.017455f) {
            ms->acc = GO;
        } else {
            ms->acc = STOP;
        }

        if (tmp == STOP && ms->acc == STOP) {
            t_star = -2 * t * t * t / t_end / t_end + 3 * t * t / t_end;
        }

        if (tmp == STOP && ms->acc == GO) {
            t_star = t * t * t * t / t_end / t_end / t_end - 3 * t * t * t / t_end / t_end + 3 * t * t / t_end;
        }

        if (tmp == GO && ms->acc == STOP) {
            t_star = -t * t * t * t / t_end / t_end / t_end + t * t * t / t_end / t_end + t;
        }

        if (tmp == GO && ms->acc == GO) {
            t_star = t;
        }

    } else {
        t_star = t;
    }
    */

    ms->s = t * ms->vel;
    return;
}

void eval_line(machine_state *ms) {
    float u;

    if (ms->s > ms->curr.arclen) return;

    u = sat(ms->s / ms->curr.arclen, 0, 1);

    ms->q = (1.0f - u) * ms->curr.ctrl[0] + u * ms->curr.ctrl[1];
    return;
}

/*

    float quad_bez_len(float complex q0, float complex q1, float complex q2) {
    float complex a, b, c;
    float na, nb, nc, kna, k, tmp, s;

    a = q0 - 2 * q1 + q2;
    b = (q1 - q0);
    c = (q2 - q1);

    nc = cabsf(c);
    na = cabsf(a);
    nb = cabsf(b);

    kna = (crealf(a) * crealf(b) + cimagf(a) * cimagf(b));
    k = kna / na;

    if (fabsf(nb * nb - k * k) < 0.00001) {
        k /= na;
        tmp = 1.0f + k;

        s = na * (sgn(tmp) * tmp * tmp - sgn(k) * k * k);
        // D = kna / na/na
        //printf("fuck\n");
    } else {
        s = (nc + (nc - nb) * k / na + (nb * nb - k * k) * logf((na + k + nc) / (k + nb)) / na);
    }

    return s;
}

void eval_quad_bez(machine_state *ms) {
    int iter;
    float complex a, b;
    float na, nb, c, d, e, f, g, eta, deta, u, tmp;

    if (ms->s > ms->curr.arclen) return;

    a = ms->curr.ctrl[0] - 2 * ms->curr.ctrl[1] + ms->curr.ctrl[2];
    b = (ms->curr.ctrl[1] - ms->curr.ctrl[0]);

    na = cabsf(a);
    nb = cabsf(b);

    e = nb * nb / na / na;
    d = (crealf(a) * crealf(b) + cimagf(a) * cimagf(b)) / na / na;

    f = e - d * d;

    if (fabsf(f) < 0.00001) {
        u = sqrt(fabsf(ms->s / na + sgn(d) * d * d)) - d;
    } else {
        c = ms->s / na + d * sqrtf(e);

        g = (d + sqrtf(e));

        eta = d;
        deta = 1.0f;
        iter = 0;
        while (fabs(deta) > 0.00001 && ++iter < MAXITER) {
            tmp = sqrtf(eta * eta + f);
            deta = 0.5 * eta + (f * logf((eta + tmp) / g) - c) / (2 * tmp);
            eta -= deta;
        }

        u = sat(eta - d, 0, 1);
    }
    tmp = 1.0f - u;

    a = tmp * ms->curr.ctrl[0] + u * ms->curr.ctrl[1];
    b = tmp * ms->curr.ctrl[1] + u * ms->curr.ctrl[2];

    ms->q = tmp * a + u * b;
    return;
}

void eval_cubic_bez(machine_state *ms) {
    int k, iter;
    float complex a, b, q0, q1, q2;
    float na, nb, c, d, e, f, g, eta, deta, u, tmp, sk, swork;

    if (ms->s > ms->curr.arclen) return;

    swork = ms->s;

    for (k = 0; k < ms->curr.iend / 2; k++) {
        q0 = ms->curr.ctrl[k * 2];
        q1 = ms->curr.ctrl[k * 2 + 1];
        q2 = ms->curr.ctrl[k * 2 + 2];
        sk = arclenqb(ms->curr.ctrl + 2 * k);  //quad_bez_len(q0, q1, q2);
        if (swork > sk) {
            swork -= sk;
        } else {
            break;
        }
    }

    a = q0 - 2 * q1 + q2;
    b = q1 - q0;

    na = cabsf(a);
    nb = cabsf(b);

    e = nb * nb / na / na;
    d = (crealf(a) * crealf(b) + cimagf(a) * cimagf(b)) / na / na;

    f = e - d * d;

    if (fabsf(f) < 0.00001) {
        u = sqrt(fabsf(swork / na + sgn(d) * d * d)) - d;
    } else {
        c = swork / na + d * sqrtf(e);

        g = (d + sqrtf(e));

        eta = d;
        deta = 1.0f;
        iter = 0;
        while (fabs(deta) > 0.00001 && iter < MAXITER) {
            iter++;
            tmp = sqrtf(eta * eta + f);
            deta = 0.5 * eta + (f * logf((eta + tmp) / g) - c) / (2 * tmp);
            eta -= deta;
        }
        printf("%d\n", iter);
        u = sat(eta - d, 0, 1);
    }

    tmp = 1.0f - u;

    a = tmp * q0 + u * q1;
    b = tmp * q1 + u * q2;
    float complex prev = ms->q;
    ms->q = tmp * a + u * b;

    if (cabsf(ms->q - prev) > 1) {
        printf("jump\n");
    }

    return;
}

*/

void low_lvl_ctrl(machine_state *ms) {
    if (ms->s <= ms->curr.arclen) {
        switch (ms->curr.type) {
            case NOCOM:
                break;
            case M:
            case m:
            case L:
            case l:
            case H:
            case h:
            case V:
            case v:
            case Z:
            case z:
                eval_line(ms);
                break;
            case Q:
            case q:
            case T:
            case t:
                evalqb(&ms->q, &ms->u, ms->s, ms->curr.ctrl);
                break;
            case C:
            case c:
            case S:
            case s:
            case A:
            case a:
                evalcb(&ms->q, &ms->u, ms->s, ms->curr.iend, ms->curr.ctrl);
                break;
            default:
                break;
        }
    }
}
/*
float complex intersect(float complex q0, float complex q1, float complex q2, float complex q3) {
    float x0, y0, x1, y1, x2, y2, x3, y3, x, y, d;
    x0 = creal(q0);
    y0 = cimag(q0);
    x1 = creal(q1);
    y1 = cimag(q1);
    x2 = creal(q2);
    y2 = cimag(q2);
    x3 = creal(q3);
    y3 = cimag(q3);

    d = (x0 - x1) * (y2 - y3) - (y0 - y1) * (x2 - x3);

    if (d == 0.0f) {
        x = (x0 + x3) / 2;
        y = (y0 + y3) / 2;
    } else {
        x = ((x0 * y1 - y0 * x1) * (x2 - x3) - (x0 - x1) * (x2 * y3 - y2 * x3)) / d;
        y = ((x0 * y1 - y0 * x1) * (y2 - y3) - (y0 - y1) * (x2 * y3 - y2 * x3)) / d;
    }

    //d == 0 the two derivatives are parallel or one of the two is zero

    return CMPLXF(x, y);
}
*/
/*
void split(float complex *P, float complex *Q, float complex *R, int n) {
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
*/
/*
float cubic_bez_len(int n, int *pidx, float complex *qctrl, float complex *cctrl) {
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

    if (max(d1, d2) > PRECISION && (n + 1) < NQUBICSMAX) {
        float complex *botctrl = (float complex *)malloc(4 * sizeof(float complex));
        float complex *topctrl = (float complex *)malloc(4 * sizeof(float complex));
        split(cctrl, botctrl, topctrl, 4);
        free(cctrl);
        return (cubic_bez_len((n + 1), pidx, qctrl, botctrl) + cubic_bez_len((n + 1), pidx, qctrl, topctrl));
    } else {
        float s = quad_bez_len(w[0], w[1], w[2]);
        qctrl[++(*pidx)] = w[1];
        qctrl[++(*pidx)] = w[2];
        free(cctrl);
        return s;
    }
}
*/

bool update_M(machine_state *ms, bool relative) {
    int offset;
    float endx, endy;
    float complex tmp;
    if (sscanf(ms->d + ms->inext, "%f,%f%n", &endx, &endy, &offset) == 2) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = M;

        ms->next.ctrl[0] = tmp;
        ms->next.ctrl[1] = endx + I * endy;

        if (relative) {
            ms->next.type = m;
            ms->next.ctrl[1] += tmp;
        }

        ms->next.iend = 1;

        ms->next.ders[0] = ms->next.ctrl[1] - ms->next.ctrl[0];

        ms->next.arclen = cabsf(ms->next.ders[0]);

        //ms->next.arclen = 0.0f;

        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

bool update_L(machine_state *ms, bool relative) {
    int offset;
    float endx, endy;
    float complex tmp;
    if (sscanf(ms->d + ms->inext, "%f,%f%n", &endx, &endy, &offset) == 2) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = L;

        ms->next.ctrl[0] = tmp;
        ms->next.ctrl[1] = endx + I * endy;

        if (relative) {
            ms->next.type = l;
            ms->next.ctrl[1] += tmp;
        }

        ms->next.iend = 1;

        ms->next.ders[0] = ms->next.ctrl[1] - ms->next.ctrl[0];

        ms->next.arclen = cabsf(ms->next.ders[0]);

        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

bool update_H(machine_state *ms, bool relative) {
    int offset;
    float endx;
    float complex tmp;
    if (sscanf(ms->d + ms->inext, "%f%n", &endx, &offset) == 1) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = h;

        ms->next.ctrl[0] = tmp;
        ms->next.ctrl[1] = endx + cimagf(tmp) * I;

        if (relative) {
            ms->next.type = h;
            ms->next.ctrl[1] = tmp + endx;
        }

        ms->next.iend = 1;

        ms->next.ders[0] = ms->next.ctrl[1] - ms->next.ctrl[0];

        ms->next.arclen = cabsf(ms->next.ders[0]);

        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

bool update_V(machine_state *ms, bool relative) {
    int offset;
    float endy;
    float complex tmp;
    if (sscanf(ms->d + ms->inext, "%f%n", &endy, &offset) == 1) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = V;

        ms->next.ctrl[0] = tmp;
        ms->next.ctrl[1] = crealf(tmp) + I * endy;

        if (relative) {
            ms->next.type = v;
            ms->next.ctrl[1] = tmp + I * endy;
        }

        ms->next.iend = 1;

        ms->next.ders[0] = ms->next.ctrl[1] - ms->next.ctrl[0];

        ms->next.arclen = cabsf(ms->next.ders[0]);

        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

bool update_Z(machine_state *ms, bool relative) {
    float endx, endy;
    float complex tmp;
    if (sscanf(ms->d, "%*c %f,%f", &endx, &endy) == 2) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = Z;

        ms->next.ctrl[0] = tmp;
        ms->next.ctrl[1] = endx + I * endy;

        if (relative) {
            ms->next.type = z;
        }

        ms->next.iend = 1;

        ms->next.ders[0] = ms->next.ctrl[1] - ms->next.ctrl[0];

        ms->next.arclen = cabsf(ms->next.ders[0]);

        //ms->next.arclen = 0.0f;

        // increment the string index
        ms->inext = ms->dc;
        return 1;
    }
    return 0;
}

bool update_Q(machine_state *ms, bool relative) {
    int offset;
    float ctrlx, ctrly, endx, endy;
    float complex tmp;
    if (sscanf(ms->d + ms->inext, "%f,%f %f,%f%n", &ctrlx, &ctrly, &endx, &endy, &offset) == 4) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = Q;

        ms->next.ctrl[0] = tmp;
        ms->next.ctrl[1] = ctrlx + I * ctrly;
        ms->next.ctrl[2] = endx + I * endy;

        if (relative) {
            ms->next.type = Q;
            ms->next.ctrl[1] += tmp;
            ms->next.ctrl[2] += tmp;
        }

        ms->next.ders[0] = 2 * (ms->next.ctrl[1] - ms->next.ctrl[0]);
        ms->next.ders[1] = 2 * (ms->next.ctrl[2] - ms->next.ctrl[1]);

        ms->next.iend = 2;

        ms->next.arclen = arclenqb(ms->next.ctrl);

        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

bool update_T(machine_state *ms, bool relative) {
    int offset, prevtype;
    float endx, endy;
    float complex tmp, dtmp;
    if (sscanf(ms->d + ms->inext, "%f,%f%n", &endx, &endy, &offset) == 2) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        dtmp = ms->curr.ders[1];
        prevtype = ms->curr.type;
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = T;

        ms->next.ctrl[0] = tmp;
        ms->next.ctrl[2] = endx + I * endy;

        if (relative) {
            ms->next.type = t;
            ms->next.ctrl[2] += tmp;
        }

        if ((prevtype == T) | (prevtype == t) | (prevtype == Q) | (prevtype == q)) {
            ms->next.ctrl[1] = ms->next.ctrl[0] + dtmp / 2;
        } else {
            ms->next.ctrl[1] = ms->next.ctrl[0];
        }

        ms->next.ders[0] = 2 * (ms->next.ctrl[1] - ms->next.ctrl[0]);
        ms->next.ders[1] = 2 * (ms->next.ctrl[2] - ms->next.ctrl[1]);

        ms->next.iend = 2;

        ms->next.arclen = arclenqb(ms->next.ctrl);

        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

bool update_C(machine_state *ms, bool relative) {
    int offset;
    float ctrl1x, ctrl1y, ctrl2x, ctrl2y, endx, endy;
    float complex tmp, *cctrl;
    if (sscanf(ms->d + ms->inext, "%f,%f %f,%f %f,%f%n", &ctrl1x, &ctrl1y, &ctrl2x, &ctrl2y, &endx, &endy, &offset) == 6) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = C;

        cctrl = (float complex *)malloc(4 * sizeof(float complex));

        cctrl[0] = tmp;
        cctrl[1] = ctrl1x + I * ctrl1y;
        cctrl[2] = ctrl2x + I * ctrl2y;
        cctrl[3] = endx + I * endy;

        if (relative) {
            ms->next.type = c;
            cctrl[1] += tmp;
            cctrl[2] += tmp;
            cctrl[3] += tmp;
        }

        ms->next.ctrl[0] = tmp;
        ms->next.iend = 0;

        //ms->next.arclen = arclencb(&ms->next.iend, ms->next.ctrl, cctrl);

        ms->next.ders[0] = 3 * (cctrl[1] - cctrl[0]);
        ms->next.ders[1] = 6 * (cctrl[2] - cctrl[1]);
        ms->next.ders[2] = 3 * (cctrl[3] - cctrl[2]);

        ms->next.arclen = arclencb(&ms->next.iend, ms->next.ctrl, cctrl);  //cubic_bez_len(1, &ms->next.iend, ms->next.ctrl, cctrl);

        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

bool update_S(machine_state *ms, bool relative) {
    int offset, prevtype;
    float ctrl2x, ctrl2y, endx, endy;
    float complex tmp, dtmp, *cctrl;
    if (sscanf(ms->d + ms->inext, "%f,%f %f,%f%n", &ctrl2x, &ctrl2y, &endx, &endy, &offset) == 4) {
        // one set of arguments was successfully read!
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        tmp = ms->curr.ctrl[ms->curr.iend];
        dtmp = ms->curr.ders[2];
        prevtype = ms->curr.type;
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = S;

        cctrl = (float complex *)malloc(4 * sizeof(float complex));

        cctrl[0] = tmp;

        cctrl[2] = ctrl2x + I * ctrl2y;
        cctrl[3] = endx + I * endy;

        if (relative) {
            ms->next.type = s;
            cctrl[2] += tmp;
            cctrl[3] += tmp;
        }

        if ((prevtype == S) | (prevtype == s) | (prevtype == C) | (prevtype == c)) {
            cctrl[1] = cctrl[0] + dtmp / 3;
        } else {
            cctrl[1] = cctrl[0];
        }

        ms->next.ctrl[0] = tmp;
        ms->next.iend = 0;

        ms->next.ders[0] = 3 * (cctrl[1] - cctrl[0]);
        ms->next.ders[1] = 6 * (cctrl[2] - cctrl[1]);
        ms->next.ders[2] = 3 * (cctrl[3] - cctrl[2]);

        ms->next.arclen = arclencb(&ms->next.iend, ms->next.ctrl, cctrl);
        // increment the string index
        ms->inext += offset;
        return 1;
    }
    return 0;
}

bool update_A(machine_state *ms, bool relative) {
    //a 11.165707,5.9868603 0 0 1 10.537136,-3.333713
    int offset, fa, fs;
    float rx, ry, phi, endx, endy;
    float complex start, end;
    if (sscanf(ms->d + ms->inext, "%f,%f %f %d %d %f,%f%n", &rx, &ry, &phi, &fa, &fs, &endx, &endy, &offset) == 7) {
        memcpy(&ms->curr, &ms->next, sizeof(path_segment));
        start = ms->curr.ctrl[ms->curr.iend];
        // initialize next path segment
        memset(&ms->next, 0, sizeof(path_segment));
        ms->next.type = A;

        end = endx + I * endy;
        if (relative) {
            ms->next.type = a;
            end += start;
        }

        ms->next.iend = 0;

        ms->next.arclen = arclenel(&ms->next.iend, ms->next.ctrl, start, end, rx, ry, phi, fa, fs);
        // increment the string index
        ms->inext += offset;

        return 1;
    }
    return 0;
}

void update_next_command(machine_state *ms) {
    // if cmd is not NOCOM: try to read another set or arguments of the current command
    switch (ms->next.type) {
        case NOCOM:
            break;
        case M:
            if (update_M(ms, false)) return;
            break;
        case m:
            if (update_M(ms, true)) return;
            break;
        case L:
            if (update_L(ms, false)) return;
            break;
        case l:
            if (update_L(ms, true)) return;
            break;
        case H:
            if (update_H(ms, false)) return;
            break;
        case h:
            if (update_H(ms, true)) return;
            break;
        case V:
            if (update_V(ms, false)) return;
            break;
        case v:
            if (update_V(ms, true)) return;
            break;
        case Q:
            if (update_Q(ms, false)) return;
            break;
        case q:
            if (update_Q(ms, true)) return;
            break;
        case T:
            if (update_T(ms, false)) return;
            break;
        case t:
            if (update_T(ms, true)) return;
            break;
        case C:
            if (update_C(ms, false)) return;
            break;
        case c:
            if (update_C(ms, true)) return;
            break;
        case S:
            if (update_S(ms, false)) return;
            break;
        case s:
            if (update_S(ms, true)) return;
            break;
        case A:
            if (update_A(ms, false)) return;
            break;
        case a:
            if (update_A(ms, true)) return;
            break;
        default:
            break;
    }
    // if this succeeds the cmd stays the same but the index is incremented
    // else: continue to skip over all non-keyword chars until the next command is found
    while (ms->inext < ms->dc) {
        switch (ms->d[ms->inext++]) {
            case 'M':
                if (update_M(ms, false)) return;
            case 'm':
                if (update_M(ms, true)) return;
            case 'L':
                if (update_L(ms, false)) return;
            case 'l':
                if (update_L(ms, true)) return;
            case 'H':
                if (update_H(ms, false)) return;
            case 'h':
                if (update_H(ms, true)) return;
            case 'V':
                if (update_V(ms, false)) return;
            case 'v':
                if (update_V(ms, true)) return;
            case 'Z':
                if (update_Z(ms, false)) return;
            case 'z':
                if (update_Z(ms, true)) return;
            case 'Q':
                if (update_Q(ms, false)) return;
            case 'q':
                if (update_Q(ms, true)) return;
            case 'T':
                if (update_T(ms, false)) return;
            case 't':
                if (update_T(ms, true)) return;
            case 'C':
                if (update_C(ms, false)) return;
            case 'c':
                if (update_C(ms, true)) return;
            case 'S':
                if (update_S(ms, false)) return;
            case 's':
                if (update_S(ms, true)) return;
            case 'A':
                if (update_A(ms, false)) return;
            case 'a':
                if (update_A(ms, true)) return;
            default:
                break;
        }
    }
    // if no more keywords are found on this string

    //reset the next path segement and set its type to NOCOM
    memcpy(&ms->curr, &ms->next, sizeof(path_segment));
    memset(&ms->next, 0, sizeof(ms->next));
    ms->next.type = NOCOM;
    return;
}

float pathlen(const char *d) {
    float pathlen;
    machine_state ms = {0};
    ms.d = d;
    ms.dc = strlen(d);
    pathlen = 0.0f;
    do {
        update_next_command(&ms);
        pathlen += ms.curr.arclen;
    } while (ms.curr.type != NOCOM);
    return pathlen;
}

void initialize(machine_state *ms, float v, const char *d) {
    memset(ms, 0, sizeof(machine_state));
    // the problem is with the intersect method: it can cast w1 very far away -> check if w1 is in the vicinity
    //ms->d = "m 22.450724,37.124724 c 1.37e-4,-10.601784 6.147335,-18.887138 16.615085,-20.089781 10.467749,-1.202643 25.257195,4.677499 33.45336,10.869091 8.196165,6.191592 9.799798,12.695211 15.412755,11.314149 5.612958,-1.381062 15.234386,-10.646138 26.637686,-14.744147 11.4033,-4.098008 24.58897,-3.0289 32.07235,0.133654 7.48338,3.162553 9.26521,8.41895 11.18067,14.699882 1.91547,6.280932 3.96447,13.586082 -0.80179,19.599727 -4.76626,6.013646 -16.34766,10.735289 -28.1968,13.497142 -11.84913,2.761853 -23.96501,3.56364 -29.533271,8.730853 -5.568258,5.167214 -4.588293,14.699607 -12.517372,18.30785 C 78.844318,103.05139 62.006753,100.73511 53.23128,93.696791 44.455806,86.658476 43.743097,74.898774 37.150679,67.370865 30.55826,59.842956 18.08542,56.546562 15.145534,51.067315 12.205648,45.588068 18.798199,37.926454 24.499883,32.135711 30.201567,26.344967 35.012596,22.42487 40.001388,21.623137 c 4.988793,-0.801734 10.156314,1.514741 15.991646,3.652873 5.835333,2.138131 12.339043,4.098153 16.971569,7.661606 4.632525,3.563453 7.394383,8.7308 13.586401,9.7553 6.192017,1.024499 15.813536,-2.093586 23.163346,-6.058083 7.3498,-3.964496 12.42809,-8.775508 18.97594,-10.200844 6.54785,-1.425336 14.56636,0.534744 18.9762,6.102609 4.40983,5.567864 5.21166,14.744333 1.55904,21.203442 -3.65263,6.45911 -11.75959,10.200784 -21.55972,12.249927 -9.80013,2.049142 -21.29227,2.405487 -30.245871,5.389945 -8.953598,2.984457 -15.367915,8.596984 -17.238856,12.338941 -1.870941,3.741956 0.801636,5.61276 -0.178182,6.949056 -0.979819,1.336295 -5.612581,2.138119 -10.780118,0.08893 C 64.055246,88.707647 58.353867,83.808024 54.389329,77.616171 50.424791,71.424319 48.197555,63.940802 44.099595,59.219173 40.001636,54.497544 34.03248,52.537523 32.56255,48.973678 c -1.469931,-3.563845 1.559042,-8.730916 6.325191,-11.269822 4.76615,-2.538906 11.2701,-2.449811 16.08077,-1.068964 4.810669,1.380846 7.928939,4.053649 11.314322,7.439032 3.385383,3.385382 7.038091,7.483543 13.497288,8.107126 6.459198,0.623583 15.724506,-2.227281 22.896139,-5.523604 7.17164,-3.296322 12.24993,-7.03822 18.04059,-9.621708 5.79065,-2.583488 12.29459,-4.009009 16.08071,-2.984609 3.78611,1.024401 4.85524,4.499068 4.72165,8.196317 -0.13359,3.697248 -1.46987,7.617014 -5.52363,11.047181 -4.05376,3.430167 -10.82403,6.369893 -17.86235,7.349942 -7.03832,0.980049 -14.34334,1.07e-4 -20.356943,1.336382 -6.013603,1.336276 -10.735219,4.988847 -15.189839,8.374361 -4.45462,3.385513 -8.641583,6.503465 -13.140954,4.944263 -4.49937,-1.559201 -9.30988,-7.795048 -13.808891,-13.051318 -4.499011,-5.25627 -8.686259,-9.532609 -8.730644,-12.517182 -0.04439,-2.984574 4.05375,-4.677282 7.260717,-3.920101 3.206967,0.757182 5.523483,3.964665 8.552592,7.216501 3.02911,3.251836 6.770826,6.54811 15.234479,6.191704 8.463654,-0.356406 21.649011,-4.365467 30.348333,-7.7306 8.69931,-3.365132 13.03859,-6.167583 15.93386,-7.192016 2.89527,-1.024433 4.40991,-0.311662 4.12059,1.182195 -0.28932,1.493858 -2.29415,3.673025 -5.85802,5.544096 -3.56388,1.87107 -8.64138,3.385412 -13.54148,4.365447 -4.9001,0.980034 -9.621569,1.425456 -15.368055,3.474561 -5.746486,2.049106 -12.516903,5.701568 -17.283335,7.216188 -4.766431,1.51462 -7.528088,0.89102 -10.290065,-0.757257 -2.761978,-1.648276 -5.523569,-4.320783 -7.661743,-6.637141 -2.138175,-2.316359 -3.652699,-4.276332 -3.563561,-4.810834 0.08914,-0.534503 1.781876,0.356412 3.51906,1.380907 1.737184,1.024494 3.519086,2.18273 -0.400567,5.256135 C 53.989111,63.584266 44.367057,68.573479 36.48232,64.608794 28.597583,60.644109 22.450586,47.726508 22.450724,37.124724 Z";

    //maybe there is a -nan that is saturated to a 0 at the end of the segment
    //ms->d = "M 90.099935,98.557671 C 73.711391,82.688222 57.322842,66.818769 47.830804,54.842438 38.338766,42.866107 35.743296,34.78307 36.114182,28.516838 c 0.370886,-6.266232 3.707987,-10.715701 7.267412,-13.533503 3.559425,-2.817803 7.341649,-4.004383 11.568326,-4.412199 4.226677,-0.407816 8.898917,-0.037 12.235775,1.223584 3.336858,1.260588 5.339182,3.411233 6.562742,5.265062 1.22356,1.853829 1.668496,3.411105 2.150543,5.265126 0.482047,1.854021 1.001083,4.004313 1.520224,6.155037 M 77.419202,28.479945 c -0.296645,-3.633902 -0.59327,-7.267562 -0.444928,-10.270849 0.148342,-3.003288 0.741632,-5.376448 2.113488,-7.452671 1.371856,-2.0762232 3.522671,-3.8562086 6.006709,-4.7459657 2.484037,-0.8897572 5.302353,-0.8897572 7.489816,-0.5189835 2.187463,0.3707737 3.744891,1.1124058 5.250768,2.0573383 1.505877,0.9449324 2.869265,2.0212928 4.092795,3.4672659 1.22353,1.445973 2.26174,3.225766 3.04039,5.042592 0.77866,1.816826 1.29773,3.67065 1.5573,6.340399 0.25957,2.66975 0.25957,6.154765 -0.18538,10.53017 -0.44494,4.375406 -1.33475,9.640104 -2.78085,16.944766 -1.4461,7.304662 -3.44821,16.647825 -5.784152,25.101764 -2.335945,8.45394 -5.005484,16.017637 -7.675223,23.581898";

    ms->d = d;

    //"V 100 H 250.000000 M 250,50 c 0.025421,0.001054 0.050843,0.002108 -0.259995,0.037500 -0.310838,0.035392 -0.957921,0.105120 -1.860008,0.182500 -0.902072,0.077380 -2.059117,0.162410 -3.490001,0.275000 -1.430899,0.112590 -3.135622,0.252738 -5.094990,0.430000 -1.959369,0.177262 -4.173383,0.391638 -6.625011,0.652500 -2.451643,0.260862 -5.140901,0.568209 -8.047491,0.927500 -2.906576,0.359291 -6.030470,0.770526 -9.342506,1.242500 -3.312037,0.471974 -6.812215,1.004688 -10.472491,1.600000 -3.660291,0.595313 -7.480651,1.253225 -11.422500,1.977500 -3.941849,0.724275 -8.005201,1.514912 -12.152508,2.370000 -4.147321,0.855088 -8.378565,1.774628 -12.647494,2.755000 -4.268914,0.980372 -8.575484,2.021574 -12.865007,3.112500 -4.289523,1.090924 -8.562014,2.231571 -12.769997,3.397499 -4.207983,1.165928 -8.351430,2.357140 -12.287498,3.545001 -3.936067,1.187861 -7.664748,2.372377 -11.272497,3.267501 -3.607750,0.895126 -7.094554,1.500864 -10.615006,1.962500 -3.520444,0.461634 -7.074520,0.779163 -10.617495,0.987500 -3.542967,0.208335 -7.074811,0.307474 -10.577499,0.319999 -3.502682,0.012524 -6.976210,-0.061564 -10.397501,-0.215000 -3.421292,-0.153439 -6.790348,-0.386227 -10.087505,-0.687500 -3.297158,-0.301275 -6.522421,-0.671033 -9.657495,-1.104999 -3.135078,-0.433967 -6.179955,-0.932144 -9.117503,-1.484999 -2.937544,-0.552855 -5.767755,-1.160391 -8.464999,-1.832500 -2.697244,-0.672109 -5.261533,-1.408795 -7.712498,-2.205001 -2.450969,-0.796208 -4.788620,-1.651939 -7.000001,-2.555000 -2.211379,-0.903061 -4.296484,-1.853453 -6.247500,-2.807500 -1.951015,-0.954046 -3.767942,-1.911746 -5.440001,-2.887500 -1.672058,-0.975752 -3.199247,-1.969555 -4.577500,-2.957500 -1.378254,-0.987943 -2.607570,-1.970029 -3.670000,-2.904999 -1.062430,-0.934972 -1.957973,-1.822829 -2.722500,-2.745000 -0.764527,-0.922171 -1.398037,-1.878654 -1.745000,-2.460000 -0.346963,-0.581346 -0.407380,-0.787554 -0.495000,-1.065000 -0.087620,-0.277445 -0.202443,-0.626128 -0.232500,-0.912500 -0.030057,-0.286372 0.024653,-0.510433 0.062500,-0.675000 0.037847,-0.164567 0.058833,-0.269638 0.507500,-0.472500 0.448667,-0.202862 1.325017,-0.503514 2.425000,-0.742500 1.099983,-0.238987 2.423599,-0.416308 4.000000,-0.580000 1.576401,-0.163692 3.405589,-0.313755 5.460000,-0.440000 2.054411,-0.126245 4.334045,-0.228672 6.827500,-0.315000 2.493454,-0.086328 5.200731,-0.156557 8.102501,-0.210000 2.901770,-0.053443 5.998040,-0.090099 9.270000,-0.115000 3.271960,-0.024901 6.719615,-0.038048 10.322500,-0.040000 3.602881,-0.001952 7.360999,0.007292 11.250000,0.025000 3.889006,0.017708 7.908892,0.043881 12.037497,0.077500 4.128605,0.033619 8.365922,0.074685 12.685001,0.120000 4.319079,0.045315 8.719929,0.094879 13.177499,0.150000 4.457571,0.055121 8.971870,0.115799 13.515003,0.177500 4.543140,0.061701 9.115092,0.124425 13.690002,0.190000 4.574895,0.065575 9.152740,0.134001 13.705000,0.202500 4.552260,0.068499 9.078964,0.137070 13.554991,0.205000 4.476026,0.067930 8.901373,0.135220 13.247505,0.202500 4.346132,0.067280 8.613006,0.134551 12.777508,0.200000 4.164487,0.065449 8.226588,0.129077 12.159988,0.192500 3.933415,0.063423 7.738188,0.126643 11.392504,0.187500 3.654316,0.060857 7.158235,0.119352 10.490000,0.177500 3.331766,0.058148 6.491423,0.115948 9.460002,0.172500 2.968595,0.056552 5.746141,0.111855 8.314997,0.167500 2.568856,0.055645 4.929051,0.111631 7.064998,0.170000 2.135947,0.058369 4.047692,0.119121 5.727500,0.185000 1.679793,0.065879 3.127679,0.136885 4.340008,0.202500 1.212314,0.065615 2.189100,0.125839 2.902493,0.185000 0.713423,0.059161 1.163468,0.117260 1.372501,0.147500 0.209048,0.030240 0.177011,0.032620 0.145003,0.035000 Z";
    ms->dc = strlen(ms->d);

    //ms->d =

    //ms->qs = (char *)malloc(256 * sizeof(char));
    //strncpy(ms->qs, "M 47.306706,42.559409", 256);
    //ms->qsc = strlen(ms->qs);
    //ms->next.type = M;
    //ms->next.iend = 1;
    //ms->inext = 0;
    //ms->s = 0.0f;
    ms->vel = v;
    //ms->acc = STOP;
    update_next_command(ms);
    //read_command(&ms->inext, &ms->next, ms->dc, ms->d);
}

void superloop(machine_state *ms1, machine_state *ms2) {
    while ((ms1->curr.type != NOCOM) | (ms2->curr.type != NOCOM)) {
        if ((ms2->s >= ms2->curr.arclen)) {
            update_next_command(ms2);
            ms2->u = 0;
            ms2->ts = clock();
            //printf("%f,%f  %f,%f  %f,%f\n", creal(ms->curr.start), cimag(ms->curr.start), creal(ms->curr.ctrl1), cimag(ms->curr.ctrl1), creal(ms->curr.end), cimag(ms->curr.end));
        }
        if ((ms1->s >= ms1->curr.arclen)) {
            update_next_command(ms1);
            ms1->u = 0;
            ms1->ts = clock();
            //printf("%f,%f  %f,%f  %f,%f\n", creal(ms->curr.start), cimag(ms->curr.start), creal(ms->curr.ctrl1), cimag(ms->curr.ctrl1), creal(ms->curr.end), cimag(ms->curr.end));
        }

        high_lvl_ctrl(ms1);
        low_lvl_ctrl(ms1);
        visualize(ms1, BRIGHT_BLUE);

        high_lvl_ctrl(ms2);
        low_lvl_ctrl(ms2);
        visualize(ms2, BRIGHT_GREEN);

        //float dx = creal(ms1->q) - creal(ms2->q);
        //consider the direction of the derivative
        //ms2->vel -= dx;

        nanosleep(200 / 1000, 0);
    }
    //end of string
}