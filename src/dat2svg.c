//libraries
#include <complex.h>
#include <libgen.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//constants
#define NMAX 256
#define PI 3.141592654

void parseaf(int *n, float complex *q, char *fname) {
    FILE *fptr;

    float tmpx, tmpy;
    int i = 0;

    fptr = fopen(fname, "r");
    if (fptr == NULL) {
        exit(1);
    }
    fscanf(fptr, "%*[^\n]");

    while (fscanf(fptr, "%f %f", &tmpx, &tmpy) == 2) {
        q[i] = tmpx + I * tmpy;
        i++;
    }
    *n = i;
}

int leading_edge_idx(int n, float complex *q) {
    int i, ile;
    float xmin;
    xmin = crealf(q[0]);
    for (i = 0; i < n; i++) {
        if (crealf(q[i]) < xmin) {
            ile = i;
            xmin = crealf(q[i]);
        }
    }
    return ile;
}

void trisolve(float complex *x, float *b, float *a, float *c,
              float complex *r, int n) {
    int i;
    float tmp;
    for (i = 1; i < n; i++) {
        tmp = b[i] / a[i - 1];
        a[i] -= tmp * c[i - 1];
        r[i] -= tmp * r[i - 1];
    }
    x[n - 1] = r[n - 1] / a[n - 1];
    for (i = n - 2; i > -1; i--) {
        x[i] = (r[i] - c[i] * x[i + 1]) / a[i];
    }
}

void fitcubics(int *resc, float complex *res, float complex *q, int qc) {
    int i;
    float *a, *b, *c;
    float complex *r, *d;
    a = (float *)malloc(qc * sizeof(float));
    b = (float *)malloc(qc * sizeof(float));
    c = (float *)malloc(qc * sizeof(float));
    r = (float complex *)malloc(qc * sizeof(float complex));
    d = (float complex *)malloc(qc * sizeof(float complex));

    if ((a == NULL) | (b == NULL) | (c == NULL) | (r == NULL) | (d == NULL)) {
        exit(1);
    }

    for (i = 0; i < qc; i++) {
        b[i] = 1.0;
        c[i] = 1.0;
    }
    a[0] = 2.0;
    a[qc - 1] = 2.0;

    r[0] = 3 * (q[1] - q[0]);
    r[qc - 1] = 3 * (q[qc - 1] - q[qc - 2]);

    for (i = 1; i < qc - 1; i++) {
        a[i] = 4.0;
        r[i] = 3 * (q[i + 1] - q[i - 1]);
    }

    trisolve(d, b, a, c, r, qc);

    for (i = 0; i < qc - 1; i++) {
        res[i * 3] = q[i];
        res[i * 3 + 1] = (d[i] + 3 * q[i]) / 3;
        res[i * 3 + 2] = (3 * q[i + 1] - d[i + 1]) / 3;
    }
    res[(qc - 1) * 3] = q[qc - 1];

    *resc = (qc - 1) * 3 + 1;
    free(a);
    free(b);
    free(c);
    free(r);
    free(d);
}

void saveafsvg(char *fname, int zc, int nmid, float complex *z, float complex factor, float complex offset, bool relative, int height, int width) {
    int i, ilast;

    FILE *out;

    out = fopen(fname, "w+");

    if (out) {
        fprintf(out, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
        fprintf(out,
                "\
<svg\n\
   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n\
   xmlns:cc=\"http://creativecommons.org/ns#\"\n\
   xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n\
   xmlns:svg=\"http://www.w3.org/2000/svg\"\n\
   xmlns=\"http://www.w3.org/2000/svg\"\n\
   id=\"svg8\"\n\
   version=\"1.1\"\n\
   viewBox=\"0 0 %d %d\"\n\
   height=\"%dmm\"\n\
   width=\"%dmm\">\n\
   <defs\n\
     id=\"defs2\" />\n\
   <metadata\n\
     id=\"metadata5\">\n\
    <rdf:RDF>\n\
      <cc:Work\n\
         rdf:about=\"\">\n\
        <dc:format>image/svg+xml</dc:format>\n\
        <dc:type\n\
           rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n\
        <dc:title></dc:title>\n\
      </cc:Work>\n\
    </rdf:RDF>\n\
  </metadata>\n\
   ",
                width, height, height, width);

        fprintf(out, "<path\n");
        fprintf(out, "   id=\"path1\"\n");
        fprintf(out, "   d=\"");
        if (relative) {
            fprintf(out, "M");
            fprintf(out, " %f,%f", creal(z[0] * factor + offset), cimag(z[0] * factor + offset));
            fprintf(out, " c");
            ilast = 0;
            for (i = 1; i < nmid; i++) {
                ilast = ((i - 1) / 3) * 3;
                fprintf(out, " %f,%f", creal((z[i] - z[ilast]) * factor), cimag((z[i] - z[ilast]) * factor));
            }

            fprintf(out, " \"\n");

        } else {
            fprintf(out, "M");
            fprintf(out, " %f,%f", creal(z[0] * factor + offset), cimag(z[0] * factor + offset));
            fprintf(out, " C");
            for (i = 1; i < nmid; i++) {
                fprintf(out, " %f,%f", creal(z[i] * factor + offset), cimag(z[i] * factor + offset));
            }
            fprintf(out, " \"\n");
        }
        fprintf(out, "   style=\"fill:none;stroke:#32afff;stroke-width:0.26458332px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n");
        fprintf(out, "   />\n");

        fprintf(out, "   <path\n");
        fprintf(out, "   id=\"path2\"\n");
        fprintf(out, "   d=\"");
        if (relative) {
            fprintf(out, "M");
            fprintf(out, " %f,%f", creal(z[nmid - 1] * factor + offset), cimag(z[nmid - 1] * factor + offset));
            fprintf(out, " c");
            ilast = 0;
            for (i = nmid; i < zc; i++) {
                ilast = ((i - 1) / 3) * 3;
                fprintf(out, " %f,%f", creal((z[i] - z[ilast]) * factor), cimag((z[i] - z[ilast]) * factor));
            }

            fprintf(out, " \"\n");

        } else {
            fprintf(out, "M");
            fprintf(out, " %f,%f", creal(z[nmid - 1] * factor + offset), cimag(z[nmid - 1] * factor + offset));
            fprintf(out, " C");
            for (i = nmid; i < zc; i++) {
                fprintf(out, " %f,%f", creal(z[i] * factor + offset), cimag(z[i] * factor + offset));
            }
            fprintf(out, " \"\n");
        }
        fprintf(out, "   style=\"fill:none;stroke:#ef2929;stroke-width:0.26458332px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n");

        fprintf(out, "   />\n");

        fprintf(out, "</svg>");
        fclose(out);
    } else {
        exit(1);
    }
    return;
}

int main(int argc, char **argv) {
    int qc, zc, width, height;
    int midpoint, nmid;
    float complex q[NMAX];
    float complex z[(NMAX - 1) * 3 + 1];

    size_t len;

    float chord, angle, offset_x, offset_y;
    bool relative;
    char argument[128];
    char filename[128];

    //default args

    chord = 1.0f;
    angle = 0.0f;
    offset_x = 0.0f;
    offset_y = 0.0f;
    width = 210;
    height = 297;
    relative = false;

    for (int i = 0; i < argc; i++) {
        memcpy(argument, argv[i], strcspn(argv[i], " "));
        argument[strcspn(argv[i], " ")] = '\0';
        len = strlen(argument);

        if (strncmp(&argument[len - 4], ".dat", 4) == 0) {
            strncpy(filename, argument, len);
            strcpy(filename + len - 4, ".svg");
        }
        if (argument[0] == '-') {
            switch (argument[1]) {
                case 'c':
                    if (i + 1 < argc) sscanf(argv[i + 1], "%f", &chord);
                    break;
                case 'a':
                    if (i + 1 < argc) sscanf(argv[i + 1], "%f", &angle);
                    angle = (angle * PI) / 180;
                    break;
                case 'x':
                    if (i + 1 < argc) sscanf(argv[i + 1], "%f", &offset_x);
                    break;
                case 'y':
                    if (i + 1 < argc) sscanf(argv[i + 1], "%f", &offset_y);
                    break;
                case 'X':
                    if (i + 1 < argc) sscanf(argv[i + 1], "%d", &width);
                    break;
                case 'Y':
                    if (i + 1 < argc) sscanf(argv[i + 1], "%d", &height);
                    break;
                case 'r':
                    relative = true;
                    break;
                case 'h':
                    printf("dat2svg filename.dat\n -c (chord in mm)\n -a (angle of incidence in degrees)\n -x (x offset in mm)\n -y (y offset in mm)\n -X (width of sheet in mm)\n -Y (height of sheet in mm)\n -r to indicate that the svg path should use relative coordinates\n");
                    exit(0);
                default:
                    break;
            }
        }
    }

    if (argc > 1) {
        parseaf(&qc, q, argv[1]);

        midpoint = leading_edge_idx(qc, q);
        nmid = (midpoint)*3 + 1;

        fitcubics(&zc, z, q, qc);
        saveafsvg(filename, zc, nmid, z, CMPLX(chord * cos(-angle), chord * sin(-angle)), CMPLX(offset_x, offset_y), relative, height, width);
        //if (fp) {
        //    fprintf(fp, "M");
        //    fprintf(fp, " %f,%f", 100 * creal(z[i]), 100 * cimag(z[i]));
        //    fprintf(fp, " c");
        //    for (int i = 1; i < zc; i++) {
        //        fprintf(fp, " %f,%f", 100 * creal(z[i]), 100 * cimag(z[i]));
        //    }
        //    fprintf(fp, " z\n");
        //}
    }

    //generate svg file of airfoil
    return 0;
}