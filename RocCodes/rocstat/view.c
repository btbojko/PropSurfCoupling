#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#ifndef NOXPM
#  include <X11/xpm.h>
#endif

#include "amira.h"

#define FONT "-*-courier-medium-r-normal-*"

void view() {

	int i, j, k;

	Display *display = NULL;
	GC gc, gcs;
	XGCValues gcv;
	XSetWindowAttributes atts;
	Window window; Pixmap slice;

	/* open X display, if not already open */

	if(!display) {
		display = XOpenDisplay("");
		if(!display){
			fprintf(stderr, "Error: Can't open display %s.\n", getenv("DISPLAY"));
			exit(1);
		}
	}

	/* create main window */
	window = XCreateSimpleWindow(display, DefaultRootWindow(display),
			50, 50, mesh.size[0], mesh.size[1],   0, 8, 0);

	atts.backing_store = Always ;
	atts.backing_planes = 0xFFFFFFF;
	XChangeWindowAttributes(display, window, CWBackingStore|CWBackingPlanes, &atts);

	XStoreName(display, window, "rocstat");

	XClearWindow(display, window);
	XMapWindow(display, window);
	XFlush(display);

	gcv.function = GXcopy;
	gcv.plane_mask = 0xFFFFFFFF;
	gcv.subwindow_mode = ClipByChildren;
	gcv.clip_x_origin = gcv.clip_y_origin = 0;
	gcv.clip_mask = None;
	gcv.foreground = 1;
	gcv.background = 0;
	gcv.graphics_exposures = False;
	gcv.font = XLoadFont(display, FONT);

	gc = XCreateGC(display, window, GCForeground |GCBackground| GCPlaneMask |
			GCSubwindowMode | GCClipXOrigin |  GCClipYOrigin |
			GCClipMask | GCFunction | GCGraphicsExposures | GCFont, &gcv);

	/* create pixmap buffer for use with main window */
	slice = XCreatePixmap(display, window, mesh.size[0], mesh.size[1],
			DefaultDepth(display, DefaultScreen(display)));

	gcs   = XCreateGC(display, slice, GCForeground |GCBackground| GCPlaneMask |
			GCSubwindowMode | GCClipXOrigin |  GCClipYOrigin |
			GCClipMask | GCFunction | GCGraphicsExposures | GCFont, &gcv);

	/* show each slice on the screen */
	for (k = 0; k < mesh.size[2]; k++) {

		XSetForeground(display, gcs, 0xffffff);
		XFillRectangle(display, slice, gcs, 0, 0, mesh.size[0], mesh.size[1]);

		for (j = 0; j < mesh.size[1]; j++) {
			for (i = 0; i < mesh.size[0]; i++) {
				if ( (*mesh.voxel_content)(i, j, k) == 0){
					XSetForeground(display, gcs, 0x000000);
					XDrawPoint(display, slice, gcs, i, j);
				}
			}
		}

		XSetClipOrigin(display, gc, 0, 0);
		XCopyArea(display, slice, window, gc, 0, 0, mesh.size[0], mesh.size[1], 0, 0);
		XFlush(display);

		usleep(5000);
	}

	/* clean up before returning */
	XCloseDisplay(display); display = NULL;

	return;
}
