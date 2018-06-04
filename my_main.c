/*========== my_main.c ==========

  This is the only file you need to modify in order
  to get a working mdl project (for now).

  my_main.c will serve as the interpreter for mdl.
  When an mdl script goes through a lexer and parser,
  the resulting operations will be in the array op[].

  Your job is to go through each entry in op and perform
  the required action from the list below:

  push: push a new origin matrix onto the origin stack
  pop: remove the top matrix on the origin stack

  move/scale/rotate: create a transformation matrix
                     based on the provided values, then
                     multiply the current top of the
                     origins stack by it.

  box/sphere/torus: create a solid object based on the
                    provided values. Store that in a
                    temporary matrix, multiply it by the
                    current top of the origins stack, then
                    call draw_polygons.

  line: create a line based on the provided values. Stores
        that in a temporary matrix, multiply it by the
        current top of the origins stack, then call draw_lines.

  save: call save_extension with the provided filename

  display: view the image live
  =========================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "parser.h"
#include "symtab.h"
#include "y.tab.h"

#include "matrix.h"
#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "stack.h"
#include "gmath.h"


/*======== void first_pass() ==========
  Inputs:
  Returns:

  Checks the op array for any animation commands
  (frames, basename, vary)

  Should set num_frames and basename if the frames
  or basename commands are present

  If vary is found, but frames is not, the entire
  program should exit.

  If frames is found, but basename is not, set name
  to some default value, and print out a message
  with the name being used.
  ====================*/
void first_pass() {
  //in order to use name and num_frames throughout
  //they must be extern variables
  extern int num_frames;
  extern char name[128];
  extern struct light_node *lights;
  extern color ambient;
  extern char shading[128];

  int i;
  struct light_node *last_light = 0;

  for (i=0;i<lastop;i++) {
    //printf("%d: ",i);
    switch (op[i].opcode)
      {
	case FRAMES:
		printf("frames\n");
		if(op[i].op.frames.num_frames > 0)
			num_frames = op[i].op.frames.num_frames;
		else{
			printf("Invalid frame count\n");
			exit(1);
		}
		break;
  	case BASENAME:
		printf("basename\n");
		strncpy(name, op[i].op.basename.p->name, 128);
		break;
	case LIGHT:
		printf(""); //need space here, apparently
		struct light_node *tmp = (struct light_node *)calloc(1, sizeof(struct light_node));
		strncpy(tmp->name, op[i].op.light.p->name, 128);
		tmp->c.red = op[i].op.light.c[0]; tmp->c.green = op[i].op.light.c[1]; tmp->c.blue = op[i].op.light.c[2];
		
		if(!last_light)
			last_light->next = tmp;
		else
			lights = tmp;
		last_light = tmp;	
		break;
	case AMBIENT:
		ambient.red = op[i].op.ambient.c[0]; ambient.green = op[i].op.light.c[1]; ambient.blue = op[i].op.light.c[2];
		break;
	case SHADING:
		strncpy(shading, op[i].op.shading.p->name, 128);
		break;
     }
  }

  //set default light if none found
  if(!lights){
  lights = (struct light_node *)calloc(1, sizeof(struct light_node));

  strncpy(lights->name, "default", 128);

  lights->location[0] = 0.5;
  lights->location[1] = 0.75;
  lights->location[2] = 1;

  printf("call\n");
  lights->c.red = 0;
  lights->c.green = 155;
  lights->c.blue = 255;
  }
}

/*======== struct vary_node ** second_pass() ==========
  Inputs:
  Returns: An array of vary_node linked lists

  In order to set the knobs for animation, we need to keep
  a seaprate value for each knob for each frame. We can do
  this by using an array of linked lists. Each array index
  will correspond to a frame (eg. knobs[0] would be the first
  frame, knobs[2] would be the 3rd frame and so on).

  Each index should contain a linked list of vary_nodes, each
  node contains a knob name, a value, and a pointer to the
  next node.

  Go through the opcode array, and when you find vary, go
  from knobs[0] to knobs[frames-1] and add (or modify) the
  vary_node corresponding to the given knob with the
  appropirate value.
  ====================*/
struct vary_node ** second_pass() {
	struct vary_node **knobs = (struct vary_node **)calloc( num_frames, sizeof(struct vary_node *) );
	int i;
	int j;

  	for (i=0;i<lastop;i++) {
    	switch (op[i].opcode)
      		{
		case VARY:
			for(j = op[i].op.vary.start_frame; j <= op[i].op.vary.end_frame; j++){
				struct vary_node *last_node = 0;
				int is_not_last = 1;
				
				if(knobs[j] != 0){
					last_node = knobs[j];
					while(is_not_last){
						if(last_node->next == 0){
							is_not_last = 0;
						}
						else{
							last_node = last_node->next;
						}
					}
				}	

				double frame_range = op[i].op.vary.end_frame - op[i].op.vary.start_frame;
				double value_range = op[i].op.vary.end_val - op[i].op.vary.start_val;
				char *knob_name = op[i].op.vary.p->name;
				
				struct vary_node *vn = (struct vary_node *)calloc(1, sizeof(struct vary_node));

				strncpy(vn->name, knob_name, 128);
				vn->next = 0;
				vn->value = op[i].op.vary.start_val +
					( (value_range / frame_range) * (j - op[i].op.vary.start_frame) );
				if (last_node){
					last_node->next = vn;
				}
				else
					knobs[j] = vn;
			}
		break;
		}
	}
  	return knobs;
}

/*======== void print_knobs() ==========
Inputs:
Returns:

Goes through symtab and display all the knobs and their
currnt values
====================*/
void print_knobs() {
  int i;

  printf( "ID\tNAME\t\tTYPE\t\tVALUE\n" );
  for ( i=0; i < lastsym; i++ ) {

    if ( symtab[i].type == SYM_VALUE ) {
      printf( "%d\t%s\t\t", i, symtab[i].name );

      printf( "SYM_VALUE\t");
      printf( "%6.2f\n", symtab[i].s.value);
    }
  }
}

/*======== void my_main() ==========
  Inputs:
  Returns:

  This is the main engine of the interpreter, it should
  handle most of the commadns in mdl.

  If frames is not present in the source (and therefore
  num_frames is 1, then process_knobs should be called.

  If frames is present, the enitre op array must be
  applied frames time. At the end of each frame iteration
  save the current screen to a file named the
  provided basename plus a numeric string such that the
  files will be listed in order, then clear the screen and
  reset any other data structures that need it.

  Important note: you cannot just name your files in
  regular sequence, like pic0, pic1, pic2, pic3... if that
  is done, then pic1, pic10, pic11... will come before pic2
  and so on. In order to keep things clear, add leading 0s
  to the numeric portion of the name. If you use sprintf,
  you can use "%0xd" for this purpose. It will add at most
  x 0s in front of a number, if needed, so if used correctly,
  and x = 4, you would get numbers like 0001, 0002, 0011,
  0487
  ====================*/
void my_main() {

  int i;
  int j;
  struct matrix *tmp;
  struct stack *systems;
  screen t;
  zbuffer zb;
  color g;
  g.red = 0;
  g.green = 0;
  g.blue = 0;
  double step_3d = 20;
  double theta;
  double knob_value, xval, yval, zval;

  //Lighting values here for easy access
  double light[2][3];
  double view[3];
  double areflect[3];
  double dreflect[3];
  double sreflect[3];

  ambient.red = 50;
  ambient.green = 50;
  ambient.blue = 50;

  view[0] = 0;
  view[1] = 0;
  view[2] = 1;

  areflect[RED] = 0.1;
  areflect[GREEN] = 0.1;
  areflect[BLUE] = 0.1;

  dreflect[RED] = 0.5;
  dreflect[GREEN] = 0.5;
  dreflect[BLUE] = 0.5;

  sreflect[RED] = 0.5;
  sreflect[GREEN] = 0.5;
  sreflect[BLUE] = 0.5;

  systems = new_stack();
  tmp = new_matrix(4, 1000);
  clear_screen( t );
  clear_zbuffer(zb);

  printf("first pass\n");
  strncpy(shading, "flat", 128);
  int k;
  num_frames = 1;
  first_pass();
  struct vary_node **knob_frames = second_pass();
  printf("third\n");
  for(j=0;j<num_frames;j++){
	print_knobs();
	free(systems);
	free(tmp);
  	tmp = new_matrix(4, 1000);
  	clear_screen( t );
  	clear_zbuffer(zb);
	systems = new_stack();
  	for ( k=0; k < lastsym; k++ ) {
    		if ( symtab[k].type == SYM_VALUE ) {
		int has_name = 1;
		struct vary_node *current_knob = knob_frames[j];
      		while (strncmp(current_knob->name, symtab[k].name, 128)){
			printf("name: %s\n", current_knob->name);
			if(current_knob->next == 0){
				has_name = 0;
				break;
			}
			current_knob = current_knob->next;
			printf("p: %p\n", current_knob->next);
		}
		if (has_name)
      			symtab[k].s.value = current_knob->value;
    		}
  	}
  for (i=0;i<lastop;i++) {
    //printf("%d: ",i);
    switch (op[i].opcode)
      {

      case SPHERE:
        /* printf("Sphere: %6.2f %6.2f %6.2f r=%6.2f", */
        /* 	 op[i].op.sphere.d[0],op[i].op.sphere.d[1], */
        /* 	 op[i].op.sphere.d[2], */
        /* 	 op[i].op.sphere.r); */
        if (op[i].op.sphere.constants != NULL)
          {
            //printf("\tconstants: %s",op[i].op.sphere.constants->name);
          }
        if (op[i].op.sphere.cs != NULL)
          {
            //printf("\tcs: %s",op[i].op.sphere.cs->name);
          }
        add_sphere(tmp, op[i].op.sphere.d[0],
                   op[i].op.sphere.d[1],
                   op[i].op.sphere.d[2],
                   op[i].op.sphere.r, step_3d);
        matrix_mult( peek(systems), tmp );
        draw_polygons(tmp, t, zb, view, lights, ambient,
                      areflect, dreflect, sreflect);
        tmp->lastcol = 0;
        break;
      case TORUS:
        /* printf("Torus: %6.2f %6.2f %6.2f r0=%6.2f r1=%6.2f", */
        /* 	 op[i].op.torus.d[0],op[i].op.torus.d[1], */
        /* 	 op[i].op.torus.d[2], */
        /* 	 op[i].op.torus.r0,op[i].op.torus.r1); */
        if (op[i].op.torus.constants != NULL)
          {
            //printf("\tconstants: %s",op[i].op.torus.constants->name);
          }
        if (op[i].op.torus.cs != NULL)
          {
            //printf("\tcs: %s",op[i].op.torus.cs->name);
          }
        add_torus(tmp,
                  op[i].op.torus.d[0],
                  op[i].op.torus.d[1],
                  op[i].op.torus.d[2],
                  op[i].op.torus.r0,op[i].op.torus.r1, step_3d);
        matrix_mult( peek(systems), tmp );
        draw_polygons(tmp, t, zb, view, lights, ambient,
                      areflect, dreflect, sreflect);
        tmp->lastcol = 0;
        break;
      case BOX:
        /* printf("Box: d0: %6.2f %6.2f %6.2f d1: %6.2f %6.2f %6.2f", */
        /* 	 op[i].op.box.d0[0],op[i].op.box.d0[1], */
        /* 	 op[i].op.box.d0[2], */
        /* 	 op[i].op.box.d1[0],op[i].op.box.d1[1], */
        /* 	 op[i].op.box.d1[2]); */
        if (op[i].op.box.constants != NULL)
          {
            //printf("\tconstants: %s",op[i].op.box.constants->name);
          }
        if (op[i].op.box.cs != NULL)
          {
            //printf("\tcs: %s",op[i].op.box.cs->name);
          }
        add_box(tmp,
                op[i].op.box.d0[0],op[i].op.box.d0[1],
                op[i].op.box.d0[2],
                op[i].op.box.d1[0],op[i].op.box.d1[1],
                op[i].op.box.d1[2]);
        matrix_mult( peek(systems), tmp );
        draw_polygons(tmp, t, zb, view, lights, ambient,
                      areflect, dreflect, sreflect);
        tmp->lastcol = 0;
        break;
      case LINE:
        /* printf("Line: from: %6.2f %6.2f %6.2f to: %6.2f %6.2f %6.2f",*/
        /* 	 op[i].op.line.p0[0],op[i].op.line.p0[1], */
        /* 	 op[i].op.line.p0[1], */
        /* 	 op[i].op.line.p1[0],op[i].op.line.p1[1], */
        /* 	 op[i].op.line.p1[1]); */
        if (op[i].op.line.constants != NULL)
          {
            //printf("\n\tConstants: %s",op[i].op.line.constants->name);
          }
        if (op[i].op.line.cs0 != NULL)
          {
            //printf("\n\tCS0: %s",op[i].op.line.cs0->name);
          }
        if (op[i].op.line.cs1 != NULL)
          {
            //printf("\n\tCS1: %s",op[i].op.line.cs1->name);
          }
        add_edge(tmp,
                 op[i].op.line.p0[0],op[i].op.line.p0[1],
                 op[i].op.line.p0[2],
                 op[i].op.line.p1[0],op[i].op.line.p1[1],
                 op[i].op.line.p1[2]);
        matrix_mult( peek(systems), tmp );
        draw_lines(tmp, t, zb, g);
        tmp->lastcol = 0;
        break;
      case MOVE:
        xval = op[i].op.move.d[0];
        yval = op[i].op.move.d[1];
        zval = op[i].op.move.d[2];
        printf("Move: %6.2f %6.2f %6.2f",
               xval, yval, zval);
        if (op[i].op.move.p != NULL)
          {
            printf("\tknob: %s",op[i].op.move.p->name);
    	    if ( op[i].op.move.p->type == SYM_VALUE ) {
		xval *= op[i].op.move.p->s.value;
		yval *= op[i].op.move.p->s.value;
		zval *= op[i].op.move.p->s.value;
	    }
          }
        tmp = make_translate( xval, yval, zval );
        matrix_mult(peek(systems), tmp);
        copy_matrix(tmp, peek(systems));
        tmp->lastcol = 0;
        break;
      case SCALE:
        xval = op[i].op.scale.d[0];
        yval = op[i].op.scale.d[1];
        zval = op[i].op.scale.d[2];
        printf("Scale: %6.2f %6.2f %6.2f",
               xval, yval, zval);
        if (op[i].op.scale.p != NULL)
          {
            printf("\tknob: %s",op[i].op.scale.p->name);
    	    if ( op[i].op.scale.p->type == SYM_VALUE ) {
		xval *= op[i].op.scale.p->s.value;
		yval *= op[i].op.scale.p->s.value;
		zval *= op[i].op.scale.p->s.value;
	    }
          }
        tmp = make_scale( xval, yval, zval );
        matrix_mult(peek(systems), tmp);
        copy_matrix(tmp, peek(systems));
        tmp->lastcol = 0;
        break;
      case ROTATE:
        xval = op[i].op.rotate.axis;
        theta = op[i].op.rotate.degrees;
        printf("Rotate: axis: %6.2f degrees: %6.2f",
               xval, theta);
        if (op[i].op.rotate.p != NULL)
          {
            printf("\tknob: %s",op[i].op.rotate.p->name);
    	    if ( op[i].op.rotate.p->type == SYM_VALUE ) {
		theta *= op[i].op.rotate.p->s.value;
	    }
          }
        theta*= (M_PI / 180);
        if (op[i].op.rotate.axis == 0 )
          tmp = make_rotX( theta );
        else if (op[i].op.rotate.axis == 1 )
          tmp = make_rotY( theta );
        else
          tmp = make_rotZ( theta );

        matrix_mult(peek(systems), tmp);
        copy_matrix(tmp, peek(systems));
        tmp->lastcol = 0;
        break;
      case PUSH:
        //printf("Push");
        push(systems);
        break;
      case POP:
        //printf("Pop");
        pop(systems);
        break;
      case SAVE:
        //printf("Save: %s",op[i].op.save.p->name);
        save_extension(t, op[i].op.save.p->name);
        break;
      case DISPLAY:
        //printf("Display");
        display(t);
        break;
      } //end opcode switch
    printf("\n");
  }//end operation loop
  //if(num_frames > 1){
  //  mkdir(name, 0777);
  //  char file_name[128];
  //  sprintf(file_name, "./%s/%s%03d.png", name, name, j);
  //  printf("%s\n", file_name);
  //  save_extension(t, file_name);
  //}
    
  }
  if(num_frames > 1)
    make_animation(name);

}
