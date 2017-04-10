/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

#include "utils.h"

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
int MAX_DEPTH;

void buildScene(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

 ///////////////////////////////////////
 // TO DO: For Assignment 3 you have to use
 //        the simple scene provided
 //        here, but for Assignment 4 you
 //        *MUST* define your own scene.
 //        Part of your mark will depend
 //        on how nice a scene you
 //        create. Use the simple scene
 //        provided as a sample of how to
 //        define and position objects.
 ///////////////////////////////////////

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 // Let's add a plane
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 // o=newPlane(.05,.75,.05,0.5,.55,.8,.75,1,1,2);  // Note the plane is highly-reflective (rs=rg=.75) so we
 //                        // should see some reflections if all is done properly.
 //                        // Colour is close to cyan, and currently the plane is
 //                        // completely opaque (alpha=1). The refraction index is
 //                        // meaningless since alpha=1
 // Scale(o,6,6,1);                // Do a few transforms...
 // RotateZ(o,PI/1.20);
 // RotateX(o,PI/2.25);
 // Translate(o,0,-3,10);
 // invert(&o->T[0][0],&o->Tinv[0][0]);        // Very important! compute
 //                        // and store the inverse
 //                        // transform for this object!
 // insertObject(o,&object_list);          // Insert into object list

 // // Let's add a couple spheres

 // o=newSphere(1,0,0,0,1,.25,.25,1,1,50);        // Scene signature
 // // o=newSphere(.05,.95,.0,.35,1,.25,.25,1,1,50);     // Ambient and diffuse 
 // // o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,50);    // Full Phong
 // Scale(o,.75,.5,1.5);
 // RotateY(o,PI/2);
 // Translate(o,-1.45,1.1,3.5);
 // invert(&o->T[0][0],&o->Tinv[0][0]);
 // insertObject(o,&object_list);

 // o=newSphere(1,0,0,0,.75,.95,.55,1,1,50);        // Scene signature
 // // o=newSphere(.05,.95,.0,.75,.75,.95,.55,1,1,50);     // Ambient and diffuse 
 // // o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,50);    // Full Phong
 // Scale(o,.5,2.0,1.0);
 // RotateZ(o,PI/1.5);
 // Translate(o,1.75,1.25,5.0);
 // invert(&o->T[0][0],&o->Tinv[0][0]);
 // insertObject(o,&object_list);

 // // Insert a single point light source.
 // p.px=0;
 // p.py=15.5;
 // p.pz=-5.5;
 // p.pw=1;
 // l=newPLS(&p,.95,.95,.95);
 // insertPLS(l,&light_list);

 	// ASSIGNMENT 4: BUILD OWN SCENE

 	// Background wall
    o=newPlane(0.3, 0.8, 0.25, 0.1, 1.0, 1.0, 1.0, 1, 1, 4);
    Scale(o, 13, 5, 1);
    RotateY(o, PI/6);
    Translate(o, 2.55, 0.1, 8);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, "pac.ppm");
    insertObject(o, &object_list);

    // Floor
    o=newPlane(0.3, 0.8, 0.3, 0.2, 1.0, 1.0, 1.0, 1, 1, 4);
    Scale(o, 5, 13, 1);
    RotateZ(o, PI/3);
    RotateX(o, PI/2);
    Translate(o, -5.0, -4.0, 8.0);
    invert(&o->T[0][0], &o->Tinv[0][0]);
    loadTexture(o, "blue.ppm");
    insertObject(o, &object_list);

    // Pacman figure
    o=newSphere(0.25, 0.3, 0.4, 0.2, 0.8, 0.8, 0.8, 0.1, 0.5, 6);
    Scale(o,0.65,0.65,0.65);
    Translate(o,-2.3,0.3,4.0);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    RotateY(o, -PI/3);
    loadTexture(o, "face.ppm");
    insertObject(o,&object_list);

    // Pacman food
    o=newSphere(0.25, 0.3, 0.4, 0.2, 0.8, 0.8, 0.8, 0.1, 1.0, 6);
    Scale(o,0.3,0.3,0.3);
    Translate(o,0.0,-0.5,3.5);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    loadTexture(o, "white.ppm");
    insertObject(o,&object_list);

    // Pacman food
    o=newSphere(0.25, 0.3, 0.4, 0.2, 0.8, 0.8, 0.8, 0.1, 1.5, 6);
    Scale(o,0.3,0.3,0.3);
    Translate(o,1.5,-1,3.2);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    loadTexture(o, "white.ppm");
    insertObject(o,&object_list);

	// Back lane barrier
    o=newCylinder(.05,.45,.45,.01,0,0,1,1,1,6);
    Scale(o,0.3,12,0.3);
    RotateZ(o,PI/2.1);
    RotateY(o,PI/7);
    Translate(o,-5,-2,8);
    invert(&o->T[0][0],&o->Tinv[0][0]);            
    insertObject(o,&object_list);   

    // Front lane barrier
    o=newCylinder(.05,.45,.45,.01,0,0,1,1,1,6);
    Scale(o,0.3,1.5,0.3);
    RotateZ(o,PI/2.1);
    RotateY(o,PI/7);
    Translate(o,-2.7,-1.9,2.9);
    invert(&o->T[0][0],&o->Tinv[0][0]);            
    insertObject(o,&object_list);   

    // Light source
    addAreaLight(1., 1., 0.0, 1.0, 0.0,0.5, 10.0,-5.0, 5, 5,1, 1, 1, &object_list, &light_list);

 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //

 struct colourRGB tmp_col;  // Accumulator for colour components
 double R,G,B;          // Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj->texImg==NULL)     // Not textured, use object colour
 {
  R=obj->col.R;
  G=obj->col.G;
  B=obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }

 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////

 struct pointLS *light;

 for (light=light_list; light != NULL; light = light->next) {

  struct point3D* p0 = &light->p0;
  struct colourRGB* colour = &light->col; 

  // Shadow ray
  struct point3D *shadow_d = newPoint(p0->px - p->px, p0->py - p->py, p0->pz - p->pz);
  struct ray3D *shadow_ray = newRay(p, shadow_d);

  // Intersection
  double lambda = -1;
  struct object3D *h_obj;
  struct point3D p_i;
  struct point3D n_i;
  double a = 0;
  double b = 0;
  findFirstHit(shadow_ray, &lambda, obj, &h_obj, &p_i, &n_i, &a, &b);
  
  // Check for intersection
  if (0 < lambda && lambda < 1){ // Intersects
    // Ambient
    tmp_col.R += 0;
    tmp_col.G += 0;
    tmp_col.B += 0;   
  } else { // No intersection

    // View point ray
    struct point3D * V = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
    V->pw = 0;
    normalize(V);

    // Normal ray
    struct point3D * N = newPoint(n->px, n->py, n->pz);
    N->pw = 0;
    normalize(N);

    // Incident light ray
    struct point3D * L = newPoint(p0->px-p->px, p0->py-p->py, p0->pz-p->pz);
    L->pw = 0;
    normalize(L);

    // Reflected light ray
    double LN = dot(L, N);
    struct point3D * Rl = newPoint(2 * LN * N->px - L->px, 2 * LN * N->py - L->py, 2 * LN * N->pz - L->pz);
    Rl->pw = 0;
    normalize(Rl); 
    
    // Phong
    double max1 = max(0, dot(N, L));
    double max2 = max(0, pow(dot(V, Rl), obj->shinyness));
    double phong = obj->alb.ra + (obj->alb.rd * max1) + (obj->alb.rs * max2); // Calculate phong with albedos a, d and s

    tmp_col.R += R * (phong * colour->R);
    tmp_col.G += G * (phong * colour->G);
    tmp_col.B += B * (phong * colour->B);
   
    free(L);
    free(N);
    free(Rl);
    free(V);
  }
 }

col->R = tmp_col.R;
col->G = tmp_col.G;
col->B = tmp_col.B;

 if (depth < MAX_DEPTH){

  struct point3D *reverse =newPoint(0, 0, 0);
  reverse->pw = 0;

  struct point3D* ray_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
  ray_d->pw = 0;
  normalize(ray_d);

  // Reflective ray
  normalize(n);
  double rdn = dot(ray_d, n);
  reverse->px = 2 * rdn * n->px - ray_d->px;
  reverse->py = 2 * rdn * n->py - ray_d->py;
  reverse->pz = 2 * rdn * n->pz - ray_d->pz;
  reverse->pw = 0;
  normalize(reverse);

  // ASSIGNMENT 4: GLOSSY REFLECTION
  struct point3D * u = cross(reverse, n);
  normalize(u);

  struct point3D * v = cross(reverse, u);
  normalize(v);

  double reflect_R = 0;
  double reflect_G = 0;
  double reflect_B = 0;

  int num_reflect = 0;
  while (num_reflect < 10) {
    struct point3D* reflect = newPoint(reverse->px, reverse->py, reverse->pz);
    reflect->pw = reverse->pw;

    double theta =  2 * PI * ((double)rand() / (double)RAND_MAX * 0.05);
    double phi = 2 * PI * ((double)rand() / (double)RAND_MAX * 0.05);

    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);

    // Into world coordinates
    reflect->px =  x * u->px + y * v->px + z * reflect->px;
    reflect->py =  x * u->py + y * v->py + z * reflect->py;
    reflect->pz =  x * u->pz + y * v->pz + z * reflect->pz;
    reflect->pw =  0;
    normalize(reflect);

    // Angle
    if (0 <= dot(n, reflect)) {
      struct ray3D *reflect_p = newRay(p, reflect);

      struct colourRGB sub_color;
      sub_color.R = 0;
      sub_color.G = 0;
      sub_color.B = 0;

      rayTrace(reflect_p, depth + 1, &sub_color, obj);
      reflect_R += sub_color.R;
      reflect_G += sub_color.G;
      reflect_B += sub_color.B;

      num_reflect += 1;

      free(reflect_p);

    }

    free(reflect);
  }

  col->R += (reflect_R / num_reflect);
  col->G += (reflect_G / num_reflect);
  col->B += (reflect_B / num_reflect);
  
  free(u);
  free(v);
  free(reverse);
  free(ray_d);
}

return;


}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////

  struct object3D *current_obj = object_list;

  while (current_obj!=NULL) {
    if (current_obj != Os){

      double intersect_lamba;
      struct point3D intersect_p;
      struct point3D intersect_n;
      current_obj->intersect(current_obj, ray, &intersect_lamba, &intersect_p, &intersect_n, a, b);

      if (0 <= intersect_lamba && (*lambda == -1 || intersect_lamba < *lambda)) {
        *obj = current_obj;
        *lambda = intersect_lamba;
        *p = intersect_p;
        *n = intersect_n;
      }
    }

    current_obj = current_obj->next;
  
  }

  return;

}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //

 double lambda;     // Lambda at intersection
 double a,b;        // Texture coordinates
 struct object3D *obj;  // Pointer to object at intersection
 struct point3D p;  // Intersection point
 struct point3D n;  // Normal at intersection
 struct colourRGB I;    // Colour returned by shading function

 if (depth>MAX_DEPTH)   // Max recursion depth reached. Return invalid colour.
 {
  col->R=-1;
  col->G=-1;
  col->B=-1;
  return;
 }

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////

  // Closest intersection
  lambda = -1;
  findFirstHit(ray, &lambda, NULL, &obj, &p, &n, &a, &b);

  if (lambda < 0) {
    return;
  } 

  if(obj == NULL){
    return;
  }

  if (obj != Os) {
    rtShade(obj, &p, &n, ray, depth, a, b, &I);

    if (Os != NULL) {

      col->R += I.R * Os->alb.rg;
      col->G += I.G * Os->alb.rg;
      col->B += I.B * Os->alb.rg;

    } else {

      col->R += I.R;
      col->G += I.G;
      col->B += I.B;

    } 
  }

  return;
  
}

int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;  // Will hold the raytraced image
 struct view *cam;  // Camera and view for this scene
 int sx;        // Size of the raytraced image
 int antialiasing;  // Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];    // Name of the output file for the raytraced .ppm image
 struct point3D e;      // Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;         // Increase along u and v directions for pixel coordinates
 struct point3D pc,d;       // Point structures to keep the coordinates of a pixel and
                // the direction or a ray
 struct ray3D *ray;     // Structure to keep the ray from e to a pixel
 struct colourRGB col, sub_col;      // Return colour for raytraced pixels
 struct colourRGB background;   // Background colour
 int i,j;           // Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 3, you can use
 //        the simple scene already provided. But
 //        for Assignment 4 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();      // Create a scene. This defines all the
            // objects in the world of the raytracer

 //////////////////////////////////////////
 // TO DO: For Assignment 3 you can use the setup
 //        already provided here. For Assignment 4
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-1;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
 g.px=0;
 g.py=0;
 g.pz=1;
 g.pw=0;

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=0;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);      // du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);     // here we use wl, wt, and wsize. du=dv since the image is
                // and dv is negative since y increases downward in pixel
                // coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 fprintf(stderr,"Rendering row: ");
 for (j=0;j<sx;j++)     // For each of the pixels in the image
 {
  fprintf(stderr,"%d/%d, ",j,sx);
  for (i=0;i<sx;i++)
  {

      if (antialiasing) {
        struct colourRGB colour;
        colour.R = 0.0;
        colour.G = 0.0;
        colour.B = 0.0;
        
        int k, l;
        for(k=0; k<3; k++) {
          for (l=0; l<3; l++) {
            struct point3D pc, d; 

            // Pixel coordinate
            pc.px = cam->wl+(i+(k*0.25) + (float)rand() / (float)RAND_MAX * 0.25)*du;
            pc.py = cam->wt+(j+(l*0.25) + (float)rand() / (float)RAND_MAX * 0.25)*dv;
            pc.pz = cam->f;
            pc.pw = 1;

            // Ray direction
            d.px = pc.px - cam->e.px;
            d.py = pc.py - cam->e.py;
            d.pz = pc.pz - cam->e.pz;
            d.pw = 0;

            // To world space
            matVecMult(cam->C2W, &pc);
            matVecMult(cam->C2W, &d);
            col.R = 0.0;
            col.G = 0.0;
            col.B = 0.0;

            ray = newRay(&pc, &d);
            rayTrace(ray, 0, &col, NULL);

            colour.R += col.R;
            colour.G += col.G;
            colour.B += col.B;
          }  
        }

        colour.R  =  colour.R /( 3 * 3);
        colour.G  =  colour.G /( 3 * 3);
        colour.B  =  colour.B /( 3 * 3);

        if (colour.R < 0) {

          rgbIm[(j*sx + i)*3] = background.R * 255;
          rgbIm[(j*sx + i)*3 + 1] = background.G * 255;
          rgbIm[(j*sx + i)*3 + 2] = background.B * 255;

        } else {

          if (1 < colour.R){
            rgbIm[(j*sx + i)*3] = 255;
          } else {
            rgbIm[(j*sx + i)*3] = colour.R * 255;
          }

          if (1 < colour.G){
            rgbIm[(j*sx + i)*3+1] = 255;
          } else {
            rgbIm[(j*sx + i)*3+1] = colour.G * 255;
          }

          if (1 < colour.B){
            rgbIm[(j*sx + i)*3+2] = 255;
          } else {
            rgbIm[(j*sx + i)*3+2] = colour.B * 255;
          }

        }

      } else {

        // Pixel coordinate
        pc = *newPoint(cam->wl+(i+0.5)*du, cam->wt+(j+0.5)*dv, cam->f);

        // Ray direction
        d = *newPoint(pc.px-cam->e.px, pc.py-cam->e.py, pc.pz-cam->e.pz);
        d.pw = 0.0;

        // To world space
        matVecMult(cam->C2W, &pc);
        matVecMult(cam->C2W, &d);

        col.R = 0.0;
        col.G = 0.0;
        col.B = 0.0;

        ray = newRay(&pc, &d);
        rayTrace(ray, 0, &col, NULL);

        if (col.R < 0) {

          rgbIm[(j*sx + i)*3] = background.R * 255;
          rgbIm[(j*sx + i)*3 + 1] = background.G * 255;
          rgbIm[(j*sx + i)*3 + 2] = background.B * 255;

        } else {

          if (1 < col.R){
            rgbIm[(j*sx + i)*3] = 255;
          } else {
            rgbIm[(j*sx + i)*3] = col.R * 255;
          }

          if (1 < col.G){
            rgbIm[(j*sx + i)*3+1] = 255;
          } else {
            rgbIm[(j*sx + i)*3+1] = col.G * 255;
          }

          if (1 < col.B){
            rgbIm[(j*sx + i)*3+2] = 255;
          } else {
            rgbIm[(j*sx + i)*3+2] = col.B * 255;
          }

        }

      }

  } // end for i
 } // end for j

 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list);       // Object and light lists
 deleteImage(im);               // Rendered image
 free(cam);                 // camera view
 exit(0);
}
