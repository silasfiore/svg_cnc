#include "graphic.h"

static struct {
    Display *dis;
    int screen;
    Window win;
    GC gc;
} graph;

void initgraph(int width, int height) {
    unsigned long c_border, c_background;

    c_border = RED;
    c_background = BLACK;

    graph.dis = XOpenDisplay((char *)0);

    graph.screen = DefaultScreen(graph.dis);
    graph.win = XCreateSimpleWindow(graph.dis, DefaultRootWindow(graph.dis), 0, 0, width, height, 5, c_border, c_background);
    XSetStandardProperties(graph.dis, graph.win, "", "", None, NULL, 0, NULL);
    XSelectInput(graph.dis, graph.win, ExposureMask);
    graph.gc = XCreateGC(graph.dis, graph.win, 0, 0);
    XSetBackground(graph.dis, graph.gc, WHITE);
    XSetForeground(graph.dis, graph.gc, BLACK);
    XClearWindow(graph.dis, graph.win);
    XMapRaised(graph.dis, graph.win);
}
void drawpixel(int x, int y, unsigned long color) {
    XSetForeground(graph.dis, graph.gc, color);
    XDrawPoint(graph.dis, graph.win, graph.gc, x, y);
    return;
}
void clear(void) {
    XClearWindow(graph.dis, graph.win);
    return;
}
void terminate(void) {
    XFreeGC(graph.dis, graph.gc);
    XDestroyWindow(graph.dis, graph.win);
    XCloseDisplay(graph.dis);
    exit(0);
    return;
}
