
#define NQUBICSMAX 8

typedef enum {
    NOCOM = -1,
    M = 0,
    m = 1,
    L = 2,
    l = 3,
    H = 4,
    h = 5,
    V = 6,
    v = 7,
    C = 8,
    c = 9,
    S = 10,
    s = 11,
    Q = 12,
    q = 13,
    T = 14,
    t = 15,
    A = 16,
    a = 17,
    Z = 18,
    z = 19,
} path_command;

typedef struct path_segment {
    path_command type;
    int iend;
    float complex ctrl[1 + 2 * NQUBICSMAX];
    float complex ders[3];
    float arclen;
} path_segment;

typedef struct machine_state {
    size_t dc;
    const char *d;
    int inext;

    path_segment curr;
    path_segment next;

    float s, vel;

    clock_t ts;

    float u;
    float complex q;
} machine_state;

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define sat(u, lb, ub) min(max(u, lb), ub)
#define sgn(a) ((a < 0) ? -1 : 1)

float pathlen(const char *d);

void initialize(machine_state *ms, float v, const char *d);
void superloop(machine_state *ms1, machine_state *ms2);