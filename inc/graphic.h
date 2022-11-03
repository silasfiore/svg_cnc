#include <X11/Xlib.h>
#include <X11/Xos.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <stdlib.h>

#define BLACK 0UL
#define RED 13369344UL
#define GREEN 5151238UL
#define YELLOW 12886016UL
#define BLUE 2675048655UL
#define MAGENTA 1349845115UL
#define CYAN 2550530202UL
#define WHITE 3620929743UL
#define BRIGHT_BLACK 1465188435UL
#define BRIGHT_RED 703529001UL
#define BRIGHT_GREEN 3800694836UL
#define BRIGHT_YELLOW 3925606479UL
#define BRIGHT_BLUE 2939289855UL
#define BRIGHT_MAGENTA 2142044328UL
#define BRIGHT_CYAN 3795058914UL
#define BRIGHT_WHITE 16777215UL

void initgraph(int width, int height);
void drawpixel(int x, int y, unsigned long color);
void clear(void);
void terminate(void);
