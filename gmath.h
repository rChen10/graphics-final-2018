#ifndef GMATH_H
#define GMATH_H

#include "matrix.h"
#include "ml6.h"

#define AMBIENT_C 0
#define DIFFUSE 1
#define SPECULAR 2
#define LOCATION 0
#define COLOR 1
#define RED 0
#define GREEN 1
#define BLUE 2
#define SPECULAR_EXP 4

//lighting functions
color get_lighting( double *normal, double *view, color alight, struct light_node *lights, double *areflect, double *dreflect, double *sreflect);
color calculate_ambient(color alight, double *areflect );
color calculate_diffuse(struct light_node light, double *dreflect, double *normal );
color calculate_specular(struct light_node light, double *sreflect, double *view, double *normal );
void limit_color( color * c );

//vector functions
void normalize( double *vector );
double dot_product( double *a, double *b );
double *calculate_normal(struct matrix *polygons, int i);
double *calculate_hash_normal(int x, int y, int z, struct vertex_hash *hash_table);

#endif
