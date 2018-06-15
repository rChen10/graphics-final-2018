#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "math.h"
#include "gmath.h"

/*======== void scanline_convert() ==========
  Inputs: struct matrix *points
          int i
          screen s
          zbuffer zb
  Returns:

  Fills in polygon i by drawing consecutive horizontal (or vertical) lines.

  Color should be set differently for each polygon.
  ====================*/
void scanline_convert_flat( struct matrix *points, int i, screen s, zbuffer zb, color c) {

  int top, mid, bot, y;
  int distance0, distance1, distance2;
  double x0, x1, y0, y1, y2, dx0, dx1, z0, z1, dz0, dz1;
  int flip = 0;

  z0 = z1 = dz0 = dz1 = 0;

  y0 = points->m[1][i];
  y1 = points->m[1][i+1];
  y2 = points->m[1][i+2];

  /* color c; */
  /* c.red = (23 * (i/3))%255; */
  /* c.green = (109 * (i/3))%255; */
  /* c.blue = (c.blue+(227 * (i/3)))%255; */

  //find bot, mid, top
  if ( y0 <= y1 && y0 <= y2) {
    bot = i;
    if (y1 <= y2) {
      mid = i+1;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i+1;
    }
  }//end y0 bottom
  else if (y1 <= y0 && y1 <= y2) {
    bot = i+1;
    if (y0 <= y2) {
      mid = i;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i;
    }
  }//end y1 bottom
  else {
    bot = i+2;
    if (y0 <= y1) {
      mid = i;
      top = i+1;
    }
    else {
      mid = i+1;
      top = i;
    }
  }//end y2 bottom
  //printf("ybot: %0.2f, ymid: %0.2f, ytop: %0.2f\n", (points->m[1][bot]),(points->m[1][mid]), (points->m[1][top]));
  /* printf("bot: (%0.2f, %0.2f, %0.2f) mid: (%0.2f, %0.2f, %0.2f) top: (%0.2f, %0.2f, %0.2f)\n", */

  x0 = points->m[0][bot];
  x1 = points->m[0][bot];
  z0 = points->m[2][bot];
  z1 = points->m[2][bot];
  y = (int)(points->m[1][bot]);

  distance0 = (int)(points->m[1][top]) - y;
  distance1 = (int)(points->m[1][mid]) - y;
  distance2 = (int)(points->m[1][top]) - (int)(points->m[1][mid]);

  //printf("distance0: %d distance1: %d distance2: %d\n", distance0, distance1, distance2);
  dx0 = distance0 > 0 ? (points->m[0][top]-points->m[0][bot])/distance0 : 0;
  dx1 = distance1 > 0 ? (points->m[0][mid]-points->m[0][bot])/distance1 : 0;
  dz0 = distance0 > 0 ? (points->m[2][top]-points->m[2][bot])/distance0 : 0;
  dz1 = distance1 > 0 ? (points->m[2][mid]-points->m[2][bot])/distance1 : 0;

  while ( y <= (int)points->m[1][top] ) {
    //printf("\tx0: %0.2f x1: %0.2f y: %d\n", x0, x1, y);
    draw_line(x0, y, z0, x1, y, z1, s, zb, c);

    x0+= dx0;
    x1+= dx1;
    z0+= dz0;
    z1+= dz1;
    y++;

    if ( !flip && y >= (int)(points->m[1][mid]) ) {
      flip = 1;
      dx1 = distance2 > 0 ? (points->m[0][top]-points->m[0][mid])/distance2 : 0;
      dz1 = distance2 > 0 ? (points->m[2][top]-points->m[2][mid])/distance2 : 0;
      x1 = points->m[0][mid];
      z1 = points->m[2][mid];
    }//end flip code
  }//end scanline loop
}


void scanline_convert_gouraud( struct matrix *points, int i, screen s, zbuffer zb,
		struct vertex_hash *hash_table, double *view, color ambient, struct light_node *lights,
		double *areflect, double *dreflect, double *sreflect) {

  int top, mid, bot, y;
  int distance0, distance1, distance2;
  double x0, x1, y0, y1, y2, dx0, dx1, z0, z1, dz0, dz1;
  color c1, c2, c3;
  double *n1, *n2, *n3, *dc0, *dc1, *ci, *cf;
  int flip = 0;

  z0 = z1 = dz0 = dz1 = 0;

  y0 = points->m[1][i];
  y1 = points->m[1][i+1];
  y2 = points->m[1][i+2];

  dc0 = (double *)calloc(3, sizeof(double));
  dc1 = (double *)calloc(3, sizeof(double));
  ci = (double *)calloc(3, sizeof(double));
  cf = (double *)calloc(3, sizeof(double));
  
  //find bot, mid, top
  if ( y0 <= y1 && y0 <= y2) {
    bot = i;
    if (y1 <= y2) {
      mid = i+1;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i+1;
    }
  }//end y0 bottom
  else if (y1 <= y0 && y1 <= y2) {
    bot = i+1;
    if (y0 <= y2) {
      mid = i;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i;
    }
  }//end y1 bottom
  else {
    bot = i+2;
    if (y0 <= y1) {
      mid = i;
      top = i+1;
    }
    else {
      mid = i+1;
      top = i;
    }
  }//end y2 bottom
  //printf("ybot: %0.2f, ymid: %0.2f, ytop: %0.2f\n", (points->m[1][bot]),(points->m[1][mid]), (points->m[1][top]));
  /* printf("bot: (%0.2f, %0.2f, %0.2f) mid: (%0.2f, %0.2f, %0.2f) top: (%0.2f, %0.2f, %0.2f)\n", */
  n1 = calculate_hash_normal(points->m[0][top],
		  	     points->m[1][top],
			     points->m[2][top],
			     hash_table);

  n2 = calculate_hash_normal(points->m[0][mid],
		  	     points->m[1][mid],
			     points->m[2][mid],
			     hash_table);

  n3 = calculate_hash_normal(points->m[0][bot],
		  	     points->m[1][bot],
			     points->m[2][bot],
			     hash_table);

  c1 = get_lighting(n1, view, ambient, lights, areflect, dreflect, sreflect);
  c2 = get_lighting(n2, view, ambient, lights, areflect, dreflect, sreflect);
  c3 = get_lighting(n3, view, ambient, lights, areflect, dreflect, sreflect);


  x0 = points->m[0][bot];
  x1 = points->m[0][bot];
  z0 = points->m[2][bot];
  z1 = points->m[2][bot];
  ci[0] = c3.red; ci[1] = c3.green; ci[2] = c3.blue;
  cf[0] = c3.red; cf[1] = c3.green; cf[2] = c3.blue;
  y = (int)(points->m[1][bot]);

  distance0 = (int)(points->m[1][top]) - y;
  distance1 = (int)(points->m[1][mid]) - y;
  distance2 = (int)(points->m[1][top]) - (int)(points->m[1][mid]);

  //printf("distance0: %d distance1: %d distance2: %d\n", distance0, distance1, distance2);
  dx0 = distance0 > 0 ? (points->m[0][top]-points->m[0][bot])/distance0 : 0;
  dx1 = distance1 > 0 ? (points->m[0][mid]-points->m[0][bot])/distance1 : 0;
  dz0 = distance0 > 0 ? (points->m[2][top]-points->m[2][bot])/distance0 : 0;
  dz1 = distance1 > 0 ? (points->m[2][mid]-points->m[2][bot])/distance1 : 0;
 
  dc0[0] = distance0 > 0 ? (c1.red-c3.red)/distance0 : 0;
  dc0[1] = distance0 > 0 ? (c1.green-c3.green)/distance0 : 0;
  dc0[2] = distance0 > 0 ? (c1.blue-c3.blue)/distance0 : 0;
  dc1[0] = distance1 > 0 ? (c2.red-c3.red)/distance1 : 0;
  dc1[1] = distance1 > 0 ? (c2.green-c3.green)/distance1 : 0;
  dc1[2] = distance1 > 0 ? (c2.blue-c3.blue)/distance1 : 0;

  printf("%lf, %lf, %lf\n", dc0[0], dc0[1], dc0[2]);
  while ( y <= (int)points->m[1][top] ) {
    //printf("\tx0: %0.2f x1: %0.2f y: %d\n", x0, x1, y);
    draw_gouraud_line(x0, y, z0, x1, y, z1, s, zb, ci, cf);

    x0+= dx0;
    x1+= dx1;
    z0+= dz0;
    z1+= dz1;
    ci[0]+= dc0[0];
    ci[1]+= dc0[1];
    ci[2]+= dc0[2];
    cf[0]+= dc1[0];
    cf[1]+= dc1[1];
    cf[2]+= dc1[2];
    y++;

    if ( !flip && y >= (int)(points->m[1][mid]) ) {
      flip = 1;
      dx1 = distance2 > 0 ? (points->m[0][top]-points->m[0][mid])/distance2 : 0;
      dz1 = distance2 > 0 ? (points->m[2][top]-points->m[2][mid])/distance2 : 0;
      dc1[0] = distance2 > 0 ? (c1.red-c2.red)/distance2 : 0;
      dc1[1] = distance2 > 0 ? (c1.green-c2.green)/distance2 : 0;
      dc1[2] = distance2 > 0 ? (c1.blue-c2.blue)/distance2 : 0;
      x1 = points->m[0][mid];
      z1 = points->m[2][mid];
      cf[0] = c2.red; cf[1] = c2.green; cf[2] = c2.blue;
    }//end flip code
  }//end scanline loop
  //free heap memory
  free(n1);
  free(n2);
  free(n3);
  free(dc0);
  free(dc1);
  free(ci);
  free(cf);
}


void scanline_convert_phong( struct matrix *points, int i, screen s, zbuffer zb,
		struct vertex_hash *hash_table, double *view, color ambient, struct light_node *lights,
		double *areflect, double *dreflect, double *sreflect) {

  int top, mid, bot, y;
  int distance0, distance1, distance2;
  double x0, x1, y0, y1, y2, dx0, dx1, z0, z1, dz0, dz1;
  double *n1, *n2, *n3, *dn0, *dn1, *ni, *nf;
  int flip = 0;

  z0 = z1 = dz0 = dz1 = 0;

  y0 = points->m[1][i];
  y1 = points->m[1][i+1];
  y2 = points->m[1][i+2];

  dn0 = (double *)calloc(3, sizeof(double));
  dn1 = (double *)calloc(3, sizeof(double));
  ni = (double *)calloc(3, sizeof(double));
  nf = (double *)calloc(3, sizeof(double));

  //find bot, mid, top
  if ( y0 <= y1 && y0 <= y2) {
    bot = i;
    if (y1 <= y2) {
      mid = i+1;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i+1;
    }
  }//end y0 bottom
  else if (y1 <= y0 && y1 <= y2) {
    bot = i+1;
    if (y0 <= y2) {
      mid = i;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i;
    }
  }//end y1 bottom
  else {
    bot = i+2;
    if (y0 <= y1) {
      mid = i;
      top = i+1;
    }
    else {
      mid = i+1;
      top = i;
    }
  }//end y2 bottom
  //printf("ybot: %0.2f, ymid: %0.2f, ytop: %0.2f\n", (points->m[1][bot]),(points->m[1][mid]), (points->m[1][top]));
  /* printf("bot: (%0.2f, %0.2f, %0.2f) mid: (%0.2f, %0.2f, %0.2f) top: (%0.2f, %0.2f, %0.2f)\n", */
  n1 = calculate_hash_normal(points->m[0][top],
		  	     points->m[1][top],
			     points->m[2][top],
			     hash_table);

  n2 = calculate_hash_normal(points->m[0][mid],
		  	     points->m[1][mid],
			     points->m[2][mid],
			     hash_table);

  n3 = calculate_hash_normal(points->m[0][bot],
		  	     points->m[1][bot],
			     points->m[2][bot],
			     hash_table);

  x0 = points->m[0][bot];
  x1 = points->m[0][bot];
  z0 = points->m[2][bot];
  z1 = points->m[2][bot];
  y = (int)(points->m[1][bot]);
  ni[0] = n3[0]; ni[1] = n3[1]; ni[2] = n3[2];
  nf[0] = n3[0]; nf[1] = n3[1]; nf[2] = n3[2];

  distance0 = (int)(points->m[1][top]) - y;
  distance1 = (int)(points->m[1][mid]) - y;
  distance2 = (int)(points->m[1][top]) - (int)(points->m[1][mid]);

  //printf("distance0: %d distance1: %d distance2: %d\n", distance0, distance1, distance2);
  dx0 = distance0 > 0 ? (points->m[0][top]-points->m[0][bot])/distance0 : 0;
  dx1 = distance1 > 0 ? (points->m[0][mid]-points->m[0][bot])/distance1 : 0;
  dz0 = distance0 > 0 ? (points->m[2][top]-points->m[2][bot])/distance0 : 0;
  dz1 = distance1 > 0 ? (points->m[2][mid]-points->m[2][bot])/distance1 : 0;
  dn0[0] = distance0 > 0 ? (n1[0]-n3[0])/distance0 : 0;
  dn0[1] = distance0 > 0 ? (n1[1]-n3[1])/distance0 : 0;
  dn0[2] = distance0 > 0 ? (n1[2]-n3[2])/distance0 : 0;
  dn1[0] = distance1 > 0 ? (n2[0]-n3[0])/distance1 : 0;
  dn1[1] = distance1 > 0 ? (n2[1]-n3[1])/distance1 : 0;
  dn1[2] = distance1 > 0 ? (n2[2]-n3[2])/distance1 : 0;

  while ( y <= (int)points->m[1][top] ) {
    //printf("\tx0: %0.2f x1: %0.2f y: %d\n", x0, x1, y);
    draw_phong_line(x0, y, z0, x1, y, z1, s, zb, ni, nf, view, ambient, lights, areflect, dreflect, sreflect);

    x0+= dx0;
    x1+= dx1;
    z0+= dz0;
    z1+= dz1;
    ni[0]+=dn0[0];
    ni[1]+=dn0[1];
    ni[2]+=dn0[2];
    nf[0]+=dn1[0];
    nf[1]+=dn1[1];
    nf[2]+=dn1[2];
    y++;

    if ( !flip && y >= (int)(points->m[1][mid]) ) {
      flip = 1;
      dx1 = distance2 > 0 ? (points->m[0][top]-points->m[0][mid])/distance2 : 0;
      dz1 = distance2 > 0 ? (points->m[2][top]-points->m[2][mid])/distance2 : 0;
      dn1[0] = distance2 > 0 ? (n1[0]-n2[0])/distance2 : 0;
      dn1[1] = distance2 > 0 ? (n1[1]-n2[1])/distance2 : 0;
      dn1[2] = distance2 > 0 ? (n1[2]-n2[2])/distance2 : 0;
      x1 = points->m[0][mid];
      z1 = points->m[2][mid];
      nf[0] = n2[0];
      nf[1] = n2[1];
      nf[2] = n2[2];
    }//end flip code
  }//end scanline loop
  free(n1);
  free(n2);
  free(n3);
}
/*======== void add_polygon() ==========
  Inputs:   struct matrix *points
  struct vertex_hash *hash_table
  Returns:
  Adds the vertices in points to the hash_table matrix.
  Also adds the polygon's surface normal to the vertex's normal array.
  ====================*/
struct vertex_hash *add_matrix_hash(struct matrix *polygons, struct vertex_hash *hash_table){
  int point;
  for (point=0; point < polygons->lastcol-2; point+=3) {
      hash_table = add_vertex_hash_i( polygons->m[0][point],
        		 polygons->m[1][point],
        		 polygons->m[2][point],
			 hash_table);
      hash_table = add_vertex_hash_i( polygons->m[0][point+1],
        		 polygons->m[1][point+1],
        		 polygons->m[2][point+1],
        		 hash_table);
      hash_table = add_vertex_hash_i( polygons->m[0][point+2],
                 	 polygons->m[1][point+2],
                 	 polygons->m[2][point+2],
                 	 hash_table);
      double *surface_normal = calculate_normal(polygons, point);
      normalize(surface_normal);

      hash_table = add_normal_hash_i( polygons->m[0][point],
        		 polygons->m[1][point],
        		 polygons->m[2][point],
			 surface_normal,
			 hash_table);
      hash_table = add_normal_hash_i( polygons->m[0][point+1],
        		 polygons->m[1][point+1],
        		 polygons->m[2][point+1],
			 surface_normal,
        		 hash_table);
      hash_table = add_normal_hash_i( polygons->m[0][point+2],
                 	 polygons->m[1][point+2],
                 	 polygons->m[2][point+2],
			 surface_normal,
                 	 hash_table); 
  }
      return hash_table;
}

/*======== void add_polygon() ==========
  Inputs:   struct matrix *surfaces
  double x0
  double y0
  double z0
  double x1
  double y1
  double z1
  double x2
  double y2
  double z2
  Returns:
  Adds the vertices (x0, y0, z0), (x1, y1, z1)
  and (x2, y2, z2) to the polygon matrix. They
  define a single triangle surface.
  ====================*/
void add_polygon( struct matrix *polygons,
                  double x0, double y0, double z0,
                  double x1, double y1, double z1,
                  double x2, double y2, double z2 ) {

  add_point(polygons, x0, y0, z0);
  add_point(polygons, x1, y1, z1);
  add_point(polygons, x2, y2, z2);
}

/*======== void draw_polygons() ==========
  Inputs:   struct matrix *polygons
  screen s
  color c
  Returns:
  Goes through polygons 3 points at a time, drawing
  lines connecting each points to create bounding
  triangles
  ====================*/
void draw_polygons(struct matrix *polygons, screen s, zbuffer zb,
                   double *view, struct light_node *lights, color ambient,
                   double *areflect,
                   double *dreflect,
                   double *sreflect,
		   char *shading,
		   struct vertex_hash *hash_table) {
  if ( polygons->lastcol < 3 ) {
    printf("Need at least 3 points to draw a polygon!\n");
    return;
  }
	printf("%p", hash_table);

  int point;
  double *normal;
  struct vertex_hash *poly_hash = NULL;

  printf("%s\n", shading);
  //use surface normal if flat
  //vertex normal if phong or gouraudshading
  if(!strncmp(shading, "flat", 128)){
  }
  else{
          poly_hash = add_matrix_hash(polygons, hash_table);
  }
  
  for (point=0; point < polygons->lastcol-2; point+=3) {
    normal = calculate_normal(polygons, point);

    if ( dot_product(normal, view) > 0 ) {

      color c;

      if(!strncmp(shading, "gouraud", 128)){
	scanline_convert_gouraud( polygons, point, s, zb,
		poly_hash, view, ambient, lights,
		areflect, dreflect, sreflect);
      }
      else if(!strncmp(shading, "phong", 128)){
	scanline_convert_phong( polygons, point, s, zb,
		poly_hash, view, ambient, lights,
		areflect, dreflect, sreflect);
      }
      else{
        c = get_lighting(normal, view, ambient, lights, areflect, dreflect, sreflect);
        scanline_convert_flat(polygons, point, s, zb, c);
      }

 //     draw_line( polygons->m[0][point],
 //                polygons->m[1][point],
 //                polygons->m[2][point],
 //                polygons->m[0][point+1],
 //                polygons->m[1][point+1],
 //                polygons->m[2][point+1],
 //                s, zb, c);
 //     draw_line( polygons->m[0][point+2],
 //                polygons->m[1][point+2],
 //                polygons->m[2][point+2],
 //                polygons->m[0][point+1],
 //                polygons->m[1][point+1],
 //                polygons->m[2][point+1],
 //                s, zb, c);
 //     draw_line( polygons->m[0][point],
 //                polygons->m[1][point],
 //                polygons->m[2][point],
 //                polygons->m[0][point+2],
 //                polygons->m[1][point+2],
 //                polygons->m[2][point+2],
 //                s, zb, c);
    }
    //top level hashtable is already set to NULL
    //HASH_CLEAR(hh, poly_hash);
  }
}

/*======== void add_box() ==========
  Inputs:   struct matrix * edges
  double x
  double y
  double z
  double width
  double height
  double depth
  Returns:

  add the points for a rectagular prism whose
  upper-left corner is (x, y, z) with width,
  height and depth dimensions.
  ====================*/
void add_box( struct matrix * polygons,
              double x, double y, double z,
              double width, double height, double depth, struct vertex_hash *hash_table ) {

  double x1, y1, z1;
  x1 = x+width;
  y1 = y-height;
  z1 = z-depth;

  //front
  add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
  add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);

  //back
  add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
  add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);

  //right side
  add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
  add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
  //left side
  add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
  add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);

  //top
  add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
  add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
  //bottom
  add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
  add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);
}//end add_box

/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double cz
  double r
  double step
  Returns:

  adds all the points for a sphere with center
  (cx, cy, cz) and radius r.

  should call generate_sphere to create the
  necessary points
  ====================*/
void add_sphere( struct matrix * edges,
                 double cx, double cy, double cz,
                 double r, int step, struct vertex_hash *hash_table ) {

  struct matrix *points = generate_sphere(cx, cy, cz, r, step);

  int p0, p1, p2, p3, lat, longt;
  int latStop, longStop, latStart, longStart;
  latStart = 0;
  latStop = step;
  longStart = 0;
  longStop = step;

  step++;
  for ( lat = latStart; lat < latStop; lat++ ) {
    for ( longt = longStart; longt < longStop; longt++ ) {

      p0 = lat * (step) + longt;
      p1 = p0+1;
      p2 = (p1+step) % (step * (step-1));
      p3 = (p0+step) % (step * (step-1));

      //printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3);
      if (longt < step - 2)
        add_polygon( edges, points->m[0][p0],
                     points->m[1][p0],
                     points->m[2][p0],
                     points->m[0][p1],
                     points->m[1][p1],
                     points->m[2][p1],
                     points->m[0][p2],
                     points->m[1][p2],
                     points->m[2][p2]);
      if (longt > 0 )
        add_polygon( edges, points->m[0][p0],
                     points->m[1][p0],
                     points->m[2][p0],
                     points->m[0][p2],
                     points->m[1][p2],
                     points->m[2][p2],
                     points->m[0][p3],
                     points->m[1][p3],
                     points->m[2][p3]);
    }
  }
  free_matrix(points);
}

/*======== void generate_sphere() ==========
  Inputs:   double cx
  double cy
  double cz
  double r
  int step
  Returns: Generates all the points along the surface
  of a sphere with center (cx, cy, cz) and
  radius r.
  Returns a matrix of those points
  ====================*/
struct matrix * generate_sphere(double cx, double cy, double cz,
                                double r, int step ) {

  struct matrix *points = new_matrix(4, step * step);
  int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
  double x, y, z, rot, circ;

  rot_start = 0;
  rot_stop = step;
  circ_start = 0;
  circ_stop = step;

  for (rotation = rot_start; rotation < rot_stop; rotation++) {
    rot = (double)rotation / step;

    for(circle = circ_start; circle <= circ_stop; circle++){
      circ = (double)circle / step;

      x = r * cos(M_PI * circ) + cx;
      y = r * sin(M_PI * circ) *
        cos(2*M_PI * rot) + cy;
      z = r * sin(M_PI * circ) *
        sin(2*M_PI * rot) + cz;

      /* printf("rotation: %d\tcircle: %d\n", rotation, circle); */
      /* printf("rot: %lf\tcirc: %lf\n", rot, circ); */
      /* printf("sphere point: (%0.2f, %0.2f, %0.2f)\n\n", x, y, z); */
      add_point(points, x, y, z);
    }
  }

  return points;
}

/*======== void add_torus() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double cz
  double r1
  double r2
  double step
  Returns:

  adds all the points required to make a torus
  with center (cx, cy, cz) and radii r1 and r2.

  should call generate_torus to create the
  necessary points
  ====================*/
void add_torus( struct matrix * edges,
                double cx, double cy, double cz,
                double r1, double r2, int step, struct vertex_hash *hash_table ) {

  struct matrix *points = generate_torus(cx, cy, cz, r1, r2, step);

  int p0, p1, p2, p3, lat, longt;
  int latStop, longStop, latStart, longStart;
  latStart = 0;
  latStop = step;
  longStart = 0;
  longStop = step;

  //printf("points: %d\n", points->lastcol);

  for ( lat = latStart; lat < latStop; lat++ ) {
    for ( longt = longStart; longt < longStop; longt++ ) {
      p0 = lat * step + longt;
      if (longt == step - 1)
        p1 = p0 - longt;
      else
        p1 = p0 + 1;
      p2 = (p1 + step) % (step * step);
      p3 = (p0 + step) % (step * step);

      //printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3);
      add_polygon( edges, points->m[0][p0],
                   points->m[1][p0],
                   points->m[2][p0],
                   points->m[0][p3],
                   points->m[1][p3],
                   points->m[2][p3],
                   points->m[0][p2],
                   points->m[1][p2],
                   points->m[2][p2]);
      add_polygon( edges, points->m[0][p0],
                   points->m[1][p0],
                   points->m[2][p0],
                   points->m[0][p2],
                   points->m[1][p2],
                   points->m[2][p2],
                   points->m[0][p1],
                   points->m[1][p1],
                   points->m[2][p1]);
    }
  }
  free_matrix(points);
}
/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double cz
  double r
  int step
  Returns: Generates all the points along the surface
  of a torus with center (cx, cy, cz) and
  radii r1 and r2.
  Returns a matrix of those points
  ====================*/
struct matrix * generate_torus( double cx, double cy, double cz,
                                double r1, double r2, int step ) {

  struct matrix *points = new_matrix(4, step * step);
  int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
  double x, y, z, rot, circ;

  rot_start = 0;
  rot_stop = step;
  circ_start = 0;
  circ_stop = step;

  for (rotation = rot_start; rotation < rot_stop; rotation++) {
    rot = (double)rotation / step;

    for(circle = circ_start; circle < circ_stop; circle++){
      circ = (double)circle / step;

      x = cos(2*M_PI * rot) *
        (r1 * cos(2*M_PI * circ) + r2) + cx;
      y = r1 * sin(2*M_PI * circ) + cy;
      z = -1*sin(2*M_PI * rot) *
        (r1 * cos(2*M_PI * circ) + r2) + cz;

      //printf("rotation: %d\tcircle: %d\n", rotation, circle);
      //printf("torus point: (%0.2f, %0.2f, %0.2f)\n", x, y, z);
      add_point(points, x, y, z);
    }
  }
  return points;
}

/*======== void add_circle() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double r
  double step
  Returns:

  Adds the circle at (cx, cy) with radius r to edges
  ====================*/
void add_circle( struct matrix * edges,
                 double cx, double cy, double cz,
                 double r, int step ) {
  double x0, y0, x1, y1, t;
  int i;
  x0 = r + cx;
  y0 = cy;

  for (i=1; i<=step; i++) {
    t = (double)i/step;
    x1 = r * cos(2 * M_PI * t) + cx;
    y1 = r * sin(2 * M_PI * t) + cy;

    add_edge(edges, x0, y0, cz, x1, y1, cz);
    x0 = x1;
    y0 = y1;
  }
}


/*======== void add_curve() ==========
  Inputs:   struct matrix *points
  double x0
  double y0
  double x1
  double y1
  double x2
  double y2
  double x3
  double y3
  double step
  int type
  Returns:

  Adds the curve bounded by the 4 points passsed as parameters
  of type specified in type (see matrix.h for curve type constants)
  to the matrix points
  ====================*/
void add_curve( struct matrix *edges,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                int step, int type ) {

  double t, x, y;
  struct matrix *xcoefs;
  struct matrix *ycoefs;
  int i;

  xcoefs = generate_curve_coefs(x0, x1, x2, x3, type);
  ycoefs = generate_curve_coefs(y0, y1, y2, y3, type);

  /* print_matrix(xcoefs); */
  /* printf("\n"); */
  /* print_matrix(ycoefs); */

  for (i=1; i<=step; i++) {

    t = (double)i/step;
    x = xcoefs->m[0][0] *t*t*t + xcoefs->m[1][0] *t*t+
      xcoefs->m[2][0] *t + xcoefs->m[3][0];
    y = ycoefs->m[0][0] *t*t*t + ycoefs->m[1][0] *t*t+
      ycoefs->m[2][0] *t + ycoefs->m[3][0];

    add_edge(edges, x0, y0, 0, x, y, 0);
    x0 = x;
    y0 = y;
  }

  free_matrix(xcoefs);
  free_matrix(ycoefs);
}


/*======== void add_point() ==========
  Inputs:   struct matrix * points
  int x
  int y
  int z
  Returns:
  adds point (x, y, z) to points and increment points.lastcol
  if points is full, should call grow on points
  ====================*/
void add_point( struct matrix * points, double x, double y, double z) {

  if ( points->lastcol == points->cols )
    grow_matrix( points, points->lastcol + 100 );

  points->m[0][ points->lastcol ] = x;
  points->m[1][ points->lastcol ] = y;
  points->m[2][ points->lastcol ] = z;
  points->m[3][ points->lastcol ] = 1;
  points->lastcol++;
} //end add_point

/*======== void add_edge() ==========
  Inputs:   struct matrix * points
  int x0, int y0, int z0, int x1, int y1, int z1
  Returns:
  add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
  should use add_point
  ====================*/
void add_edge( struct matrix * points,
               double x0, double y0, double z0,
               double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
  Inputs:   struct matrix * points
  screen s
  color c
  Returns:
  Go through points 2 at a time and call draw_line to add that line
  to the screen
  ====================*/
void draw_lines( struct matrix * points, screen s, zbuffer zb, color c) {

  if ( points->lastcol < 2 ) {
    printf("Need at least 2 points to draw a line!\n");
    return;
  }
  int point;
  for (point=0; point < points->lastcol-1; point+=2)
    draw_line( points->m[0][point],
               points->m[1][point],
               points->m[2][point],
               points->m[0][point+1],
               points->m[1][point+1],
               points->m[2][point+1],
               s, zb, c);
}// end draw_lines




void draw_line(int x0, int y0, double z0,
               int x1, int y1, double z1,
               screen s, zbuffer zb, color c) {


  int x, y, d, A, B;
  int dy_east, dy_northeast, dx_east, dx_northeast, d_east, d_northeast;
  int loop_start, loop_end;
  double distance;
  double z, dz;

  //swap points if going right -> left
  int xt, yt;
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    z = z0;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    x1 = xt;
    y1 = yt;
    z1 = z;
  }

  x = x0;
  y = y0;
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);
  int wide = 0;
  int tall = 0;
  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ) { //octant 1/8
    wide = 1;
    loop_start = x;
    loop_end = x1;
    dx_east = dx_northeast = 1;
    dy_east = 0;
    d_east = A;
    distance = x1 - x;
    if ( A > 0 ) { //octant 1
      d = A + B/2;
      dy_northeast = 1;
      d_northeast = A + B;
    }
    else { //octant 8
      d = A - B/2;
      dy_northeast = -1;
      d_northeast = A - B;
    }
  }//end octant 1/8
  else { //octant 2/7
    tall = 1;
    dx_east = 0;
    dx_northeast = 1;
    distance = abs(y1 - y);
    if ( A > 0 ) {     //octant 2
      d = A/2 + B;
      dy_east = dy_northeast = 1;
      d_northeast = A + B;
      d_east = B;
      loop_start = y;
      loop_end = y1;
    }
    else {     //octant 7
      d = A/2 - B;
      dy_east = dy_northeast = -1;
      d_northeast = A - B;
      d_east = -1 * B;
      loop_start = y1;
      loop_end = y;
    }
  }

  z = z0;
  dz = (z1 - z0) / distance;
  //printf("\t(%d, %d) -> (%d, %d)\tdistance: %0.2f\tdz: %0.2f\tz: %0.2f\n", x0, y0, x1, y1, distance, dz, z);

  while ( loop_start < loop_end ) {

    plot( s, zb, c, x, y, z );
    if ( (wide && ((A > 0 && d > 0) ||
                   (A < 0 && d < 0)))
         ||
         (tall && ((A > 0 && d < 0 ) ||
                   (A < 0 && d > 0) ))) {
      y+= dy_northeast;
      d+= d_northeast;
      x+= dx_northeast;
    }
    else {
      x+= dx_east;
      y+= dy_east;
      d+= d_east;
    }
    z+= dz;
    loop_start++;
  } //end drawing loop
  plot( s, zb, c, x1, y1, z );
} //end draw_line

void draw_gouraud_line(int x0, int y0, double z0,
               int x1, int y1, double z1,
               screen s, zbuffer zb, double *ci, double *cf) {


  int x, y, d, A, B;
  int dy_east, dy_northeast, dx_east, dx_northeast, d_east, d_northeast;
  double *c, *dc;
  int loop_start, loop_end;
  double distance;
  double z, dz;

  //swap points if going right -> left
  int xt, yt;
  double *ct;
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    ct = ci;
    z = z0;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    ci = cf;
    x1 = xt;
    y1 = yt;
    z1 = z;
    cf = ct;
  }

  x = x0;
  y = y0;
  c = (double *)calloc(3, sizeof(double));
  c[0] = ci[0];
  c[1] = ci[1];
  c[2] = ci[2];
  dc = (double *)calloc(3, sizeof(double));
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);
  int wide = 0;
  int tall = 0;
  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ) { //octant 1/8
    wide = 1;
    loop_start = x;
    loop_end = x1;
    dx_east = dx_northeast = 1;
    dy_east = 0;
    d_east = A;
    distance = x1 - x;
    if ( A > 0 ) { //octant 1
      d = A + B/2;
      dy_northeast = 1;
      d_northeast = A + B;
    }
    else { //octant 8
      d = A - B/2;
      dy_northeast = -1;
      d_northeast = A - B;
    }
  }//end octant 1/8
  else { //octant 2/7
    tall = 1;
    dx_east = 0;
    dx_northeast = 1;
    distance = abs(y1 - y);
    if ( A > 0 ) {     //octant 2
      d = A/2 + B;
      dy_east = dy_northeast = 1;
      d_northeast = A + B;
      d_east = B;
      loop_start = y;
      loop_end = y1;
    }
    else {     //octant 7
      d = A/2 - B;
      dy_east = dy_northeast = -1;
      d_northeast = A - B;
      d_east = -1 * B;
      loop_start = y1;
      loop_end = y;
    }
  }

  z = z0;
  dz = (z1 - z0) / distance;
  dc[0] = (cf[0] - ci[0]) / distance;
  dc[1] = (cf[1] - ci[1]) / distance;
  dc[2] = (cf[2] - ci[2]) / distance;
  //printf("\t(%d, %d) -> (%d, %d)\tdistance: %0.2f\tdz: %0.2f\tz: %0.2f\n", x0, y0, x1, y1, distance, dz, z);

  color col;
  col.red = (int)c[0];
  col.green = (int)c[1];
  col.blue = (int)c[2];
  while ( loop_start < loop_end ) {

    plot( s, zb, col, x, y, z );
    if ( (wide && ((A > 0 && d > 0) ||
                   (A < 0 && d < 0)))
         ||
         (tall && ((A > 0 && d < 0 ) ||
                   (A < 0 && d > 0) ))) {
      y+= dy_northeast;
      d+= d_northeast;
      x+= dx_northeast;
    }
    else {
      x+= dx_east;
      y+= dy_east;
      d+= d_east;
    }
    z+= dz;
    c[0] = c[0] + dc[0] < 255 ? c[0] + dc[0] : 255;
    c[1] = c[1] + dc[1] < 255 ? c[1] + dc[1] : 255;
    c[2] = c[2] + dc[2] < 255 ? c[2] + dc[2] : 255;
    c[0] = c[0] + dc[0] > 0 ? c[0] + dc[0] : 0;
    c[1] = c[1] + dc[1] > 0 ? c[1] + dc[1] : 0;
    c[2] = c[2] + dc[2] > 0 ? c[2] + dc[2] : 0;
    col.red = (int)c[0];
    col.green = (int)c[1];
    col.blue = (int)c[2];
    loop_start++;
  } //end drawing loop
  plot( s, zb, col, x1, y1, z );
} //end draw_line

void draw_phong_line(int x0, int y0, double z0,
               int x1, int y1, double z1,
               screen s, zbuffer zb, 
	       double *ni, double *nf,
	       double *view, color ambient, struct light_node *lights,
	       double *areflect, double *dreflect, double *sreflect) {


  int x, y, d, A, B;
  int dy_east, dy_northeast, dx_east, dx_northeast, d_east, d_northeast;
  int loop_start, loop_end;
  double distance;
  double z, dz;
  double *n, *dn;

  //swap points if going right -> left
  int xt, yt;
  double *nt;
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    z = z0;
    nt = ni;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    ni = nf;
    x1 = xt;
    y1 = yt;
    z1 = z;
    nf = nt;
  }

  x = x0;
  y = y0;
  n = (double *)calloc(3, sizeof(double));
  n[0] = ni[0];
  n[1] = ni[1];
  n[2] = ni[2];
  dn = (double *)calloc(3, sizeof(double));
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);
  int wide = 0;
  int tall = 0;
  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ) { //octant 1/8
    wide = 1;
    loop_start = x;
    loop_end = x1;
    dx_east = dx_northeast = 1;
    dy_east = 0;
    d_east = A;
    distance = x1 - x;
    if ( A > 0 ) { //octant 1
      d = A + B/2;
      dy_northeast = 1;
      d_northeast = A + B;
    }
    else { //octant 8
      d = A - B/2;
      dy_northeast = -1;
      d_northeast = A - B;
    }
  }//end octant 1/8
  else { //octant 2/7
    tall = 1;
    dx_east = 0;
    dx_northeast = 1;
    distance = abs(y1 - y);
    if ( A > 0 ) {     //octant 2
      d = A/2 + B;
      dy_east = dy_northeast = 1;
      d_northeast = A + B;
      d_east = B;
      loop_start = y;
      loop_end = y1;
    }
    else {     //octant 7
      d = A/2 - B;
      dy_east = dy_northeast = -1;
      d_northeast = A - B;
      d_east = -1 * B;
      loop_start = y1;
      loop_end = y;
    }
  }

  z = z0;
  dz = (z1 - z0) / distance;

  dn[0] = (nf[0] - ni[0]) / distance;
  dn[1] = (nf[1] - ni[1]) / distance;
  dn[2] = (nf[2] - ni[2]) / distance;

  color c = get_lighting(n, view, ambient, lights, areflect, dreflect, sreflect);
  while ( loop_start < loop_end ) {

    plot( s, zb, c, x, y, z );
    if ( (wide && ((A > 0 && d > 0) ||
                   (A < 0 && d < 0)))
         ||
         (tall && ((A > 0 && d < 0 ) ||
                   (A < 0 && d > 0) ))) {
      y+= dy_northeast;
      d+= d_northeast;
      x+= dx_northeast;
    }
    else {
      x+= dx_east;
      y+= dy_east;
      d+= d_east;
    }
    z+= dz;
    n[0] += dn[0];
    n[1] += dn[1];
    n[2] += dn[2];
    normalize(n);
    c = get_lighting(n, view, ambient, lights, areflect, dreflect, sreflect);
    loop_start++;
  } //end drawing loop
  plot( s, zb, c, x1, y1, z );
} //end draw_line
