/*========== ml6.h ==========

Header file for fucntions we will use in ml6

Sets the maximum XYES and YRES for images as well
as the maximum color value you want to use.

Creates the point structure in order to represent 
a pixel as a color triple
=========================*/
#ifndef ML6_H
#define ML6_H

#define XRES 500
#define YRES 500
#define MAX_COLOR 255

#ifndef UTHASH_SELF_SPEC
#define UTHASH_SELF_SPEC
//Defining uthash here so I don't have to deal with it.
#include "uthash.h"

typedef struct{
	int x;
	int y;
	int z;
} vertex_key_t;

struct vertex_hash{
	vertex_key_t id; //doubles as vertex coordinates
	double **normals; //array of normals
	int size; //size of normal array
	UT_hash_handle hh;
};

static struct vertex_hash *add_vertex_hash_p(int *new_id, struct vertex_hash *hash_table){
	struct vertex_hash *h;
	int i;
	struct vertex_hash temp;
	memset(&temp, 0, sizeof(struct vertex_hash));
	temp.id.x = new_id[0]; temp.id.y = new_id[1]; temp.id.z = new_id[2];
	HASH_FIND(hh, hash_table, &temp.id, sizeof(vertex_key_t), h);
	if( !h ){
		h = (struct vertex_hash *)calloc(1, sizeof(struct vertex_hash));
		h->id.x = new_id[0];
		h->id.y = new_id[1];
		h->id.z = new_id[2];
		HASH_ADD(hh, hash_table, id, sizeof(vertex_key_t), h);
		//printf("%d %d %d\n", h->id.x, h->id.y, h->id.z);
	}
	return hash_table;
}

static struct vertex_hash *add_vertex_hash_i(int x, int y, int z, struct vertex_hash *hash_table){
	int *id = (int *)calloc(3, sizeof(int));
	id[0] = x; id[1] = y; id[2] = z;
	struct vertex_hash *ret = add_vertex_hash_p(id, hash_table);
	free(id);
	return ret;
}
static struct vertex_hash *add_normal_hash_p(int *id, double *normal, struct vertex_hash *hash_table){
	struct vertex_hash *h;
	int i;
	double *new_normal = (double *)calloc(3, sizeof(double *));
	struct vertex_hash temp;
	memset(&temp, 0, sizeof(struct vertex_hash));
	temp.id.x = id[0]; temp.id.y = id[1]; temp.id.z = id[2];
	HASH_FIND(hh, hash_table, &temp.id, sizeof(vertex_key_t), h);

	if( h ){
		int duplicate = 0;
		for(i = 0; i < h->size; i++){
			if(h->normals[i][0] == normal[0]
			   && h->normals[i][1] == normal[1]
			   && h->normals[i][2] == normal[2])
			   duplicate = 1;
		}
		if(!duplicate){
			h->normals = realloc(h->normals, (h->size + 1) * sizeof(double *));
			for(i = 0; i < 3; i++){
				new_normal[i] = normal[i];
			}
			h->normals[h->size] = new_normal;
			h->size++;
		}
	}
	else{
		//printf("\nERROR: Vertex does not exist!\n");
	}
	return hash_table;
}

static struct vertex_hash *add_normal_hash_i(int x, int y, int z, double *normal, struct vertex_hash *hash_table){
	int *id = (int *)calloc(3, sizeof(int));
	id[0] = x; id[1] = y; id[2] = z;
	struct vertex_hash * ret = add_normal_hash_p(id, normal, hash_table);
	free(id);
	return ret;
}

#endif

/*
  Every point has an individual int for
  each color value
*/
struct point_t {

  int red;
  int green;
  int blue;
} point_t;

/*
  We can now use color as a data type representing a point.
  eg:
  color c;
  c.red = 0;
  c.green = 45;
  c.blue = 187;
*/
typedef struct point_t color;

//Ricky: Putting light node at the highest header so I don't have to deal with it.
struct light_node{
	char name[128];
	double location[3];
	color c;
	struct light_node *next;
};

/*
  Likewise, we can use screen as a data type representing
  an XRES x YRES array of colors.
  eg:
  screen s;
  s[0][0] = c;
*/
typedef struct point_t screen[XRES][YRES];

//z-buffer is a 2d array of doubles to store z values
typedef double zbuffer[XRES][YRES];
#endif
