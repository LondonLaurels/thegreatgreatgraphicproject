#include "raytracer.h"
#include "scene_types.h"
#include "ray.h"
#include "image.h"
#include "kdtree.h"
#include <stdio.h>

#define MAX_DEPTH 1


/// acne_eps is a small constant used to prevent acne when computing intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {
  float t;
  Geometry geom= obj->geom;
  vec3 n =geom.plane.normal;//normal
  float dist=geom.plane.dist;
  if(!(dot(ray->dir,n))) return false;//no solution
  t=(dot(ray->orig,n)+dist)/dot(ray->dir,n)*(-1);

  //test in tmin<t<tmax
  if(t<ray->tmin||t>ray->tmax) return false;

  if (intersection!=NULL){
  //ok, put in the structure
  intersection-> position= (ray->orig + t*ray->dir);
  intersection-> normal= n;
  intersection-> mat= &(obj->mat); 
  ray->tmax= t;
  }
  return true;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {
  //calculate intersection
  float t;
  Geometry geom = obj->geom;
  vec3 dist= ray->orig-geom.sphere.center;
  float b=2*dot(ray->dir,dist);
  float delta= (b*b)-4*(dot(dist,dist)-(obj->geom.sphere.radius*obj->geom.sphere.radius));
  if (delta<0){
    return false;

  //if intersection: calculate t
  }else if(!delta){//only one solution
    t= -b*0.5;

    //test if tmin< t <tmax
    if (!(t>ray->tmin && t<ray->tmax)) return false;

  }else{
    float root= pow(delta,0.5);
    float t1= (-b+root)/2;
    float t2= (-b-root)/2;
    if (t1<t2) t=t1;
    else t=t2;

    //test if tmin< t <tmax
    if (!(t>ray->tmin && t<ray->tmax)){
      if(t1<t2) t=t2; else t=t1;
      if (!(t>ray->tmin && t<ray->tmax)) return false;
    }
  }
  if(intersection!=NULL){
    //if true: MAJ Intersection && ray
    intersection-> position= (ray->orig + t*ray->dir);
    intersection-> normal= normalize(intersection-> position-geom.sphere.center);
    intersection-> mat= &(obj->mat); 
    ray->tmax= t;
  }
  return true;
}

bool intersectTriangle (Ray *ray, Intersection *intersection, Object *triangle) {
  point3 p1 = triangle->geom.triangle.p1;
  point3 p2 = triangle->geom.triangle.p2;
  point3 p3 = triangle->geom.triangle.p3;
  
  vec3 n = cross((p2 - p1), (p3 - p1));
  n = normalize(n);
  float cos_theta = dot(n, ray->dir);
  
  //ray parallel
  if (cos_theta == 0) return false;
  
  float D = dot(-n, p1);

  float t = -(dot(n, ray->orig) + D) / cos_theta;
  
  if (t < 0 || t > ray->tmax || t < ray->tmin) return false;
  
  vec3 interPoint = rayAt(*ray, t);
  
  vec3 edge0 = p2 - p1;
  vec3 edge1 = p3 - p2;
  vec3 edge2 = p1 - p3;
  
  vec3 c0 = interPoint - p1;
  vec3 c1 = interPoint - p2;
  vec3 c2 = interPoint - p3;
  
  if (dot(n, cross(edge0, c0)) > 0 &&
      dot(n, cross(edge1, c1)) > 0 &&
      dot(n, cross(edge2, c2)) > 0){
    if(intersection!=NULL){
      intersection->normal = n;
      intersection->position = interPoint;
      intersection->mat = &triangle->mat;
      ray->tmax = t;
    }
    return true;
  }
  
  return false;
}


bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;

  for (Object *obj:scene->objects){
    switch (obj->geom.type){
      case SPHERE:
        hasIntersection |= intersectSphere(ray,intersection,obj);
        break;  
      case PLANE:
        hasIntersection |= intersectPlane(ray, intersection, obj);
        break;
      case TRIANGLE:
        hasIntersection |= intersectTriangle(ray, intersection, obj);      
        break;
      default:
        perror("Unknow object\n");   
    }
  }
  return hasIntersection;;
}

/* --------------------------------------------------------------------------- */
/*
 *  The following functions are coded from Cook-Torrance bsdf model description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba renderer)
 */

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha) {
  float ret;
  float cos2= (NdotH*NdotH);
  float tan2=(1-cos2)/cos2;
  float alpha2= alpha*alpha;
  ret= expf(-tan2/alpha2)/(M_PI*alpha2*cos2*cos2);
  return ret;

}

// Fresnel term computation. Implantation of the exact computation. we can use the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {
  //light reflected
  float cos2= LdotH*LdotH;
  float sin2= powf(extIOR/intIOR,2)*(1-cos2); 
  if(sin2>1){
      return 1;
  }
  float cost= sqrtf(1-sin2);
  float a= extIOR*LdotH;
  float b= intIOR*cost;
  float c= extIOR*cost;
  float d= intIOR*LdotH;
  float rs= powf(a-b,2)/powf(a+b,2);
  float rp= powf(c-d,2)/powf(c+d,2);
  return 0.5f*(rs+rp);

}


// Shadowing and masking function. Linked with the NDF. Here, Smith function, suitable for Beckmann NDF
float RDM_chiplus(float c) {
  return (c > 0.f) ? 1.f : 0.f;
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {
  float tan=sqrtf(1-powf(DdotN,2))/DdotN;
  float b= 1/(alpha*tan);
  float k= DdotH/DdotN;
  if(k>0 && b<1.6){
    return (((3.535*b)+(2.181*powf(b,2)))/(1+2.276*b+2.577*powf(b,2)));
  }
  if(k>0) return 1;

  return 0;
}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN, float alpha){
  float ret= RDM_G1(LdotH,LdotN,alpha)*RDM_G1(VdotH,VdotN,alpha);
  return ret;
}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {
  //specular term
  float d= RDM_Beckmann(NdotH, m->roughness);
  float f= RDM_Fresnel(LdotH, 1, m->IOR);
  float g= RDM_Smith(LdotH,LdotN,VdotH,VdotN, m-> roughness);
  return m->specularColor*((d*f*g)/(4*LdotN*VdotN)); 
}

// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {
  return m->diffuseColor/(float)M_PI;
}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {

  color3 bsdf= RDM_bsdf_d(m)+ RDM_bsdf_s(LdotH,NdotH,VdotH,LdotN,VdotN,m);
  return bsdf;

}

/* --------------------------------------------------------------------------- */

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat ){
  color3 ret = color3(0.f);
  float cos= dot(l,n);
  if(cos<0)return ret;
  //ret= (mat->diffuseColor/((float) M_PI))*cos*lc;
  //BSDF
  vec3 h= (v+l)/length(v+l);
  float lh= dot(l,h);
  float nh=dot(n,h);
  float vh=dot(v,h);
  float vn=dot(v,n);
  ret= lc*RDM_bsdf(lh,nh,vh,cos,vn,mat)*cos;
  return ret;
}

//! if tree is not null, use intersectKdTree to compute the intersection instead of intersect scene
color3 trace_ray(Scene * scene, Ray *ray, KdTree *tree) {  
  color3 direct_color = color3(0.f);
  color3 reflect_color = color3(0.f);
  color3 final_color = scene->skyColor;
  Intersection intersection;

  // version 4
  if (intersectScene(scene, ray, &intersection))
  {
    vec3 v = -(ray->dir);
    vec3 n = intersection.normal;

    // Trace ray to each light to compute shades 
    for (Light *light :scene->lights)
    {
      //create shade ray
      Ray shade_r;
      point3 l =light->position;
      vec3 p = intersection.position;
      vec3 lp = normalize(l-p);
      point3 o = p+acne_eps*lp;
      rayInit(&shade_r, o, lp, 0, distance(l,p));
      if (!(intersectScene(scene, &shade_r, NULL)))
      {
        color3 lc = light->color;
        direct_color += shade(n,v,lp,lc,intersection.mat);
      }
    }

    // Compute reflect ray
    if (ray->depth < MAX_DEPTH)
    {
      Ray reflect_ray;
      vec3 reflect_dir;
      point3 o_r = intersection.position + acne_eps*n;

      reflect_dir = reflect(ray->dir, n);
      rayInit(&reflect_ray, o_r, reflect_dir, 0, ray->tmax, ray->depth+1);

      // Compute reflect color
      reflect_color = trace_ray(scene, &reflect_ray, tree);

      float LH = dot(reflect_ray.dir, n);
      float fresnel_term = RDM_Fresnel(LH, 1, intersection.mat->IOR);

      reflect_color = reflect_color * fresnel_term;
    }

    final_color = direct_color + reflect_color;
  }

  return final_color;
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing and kdtree initializaion
  float aspect = 1.f/scene->cam.aspect;
    
  KdTree *tree =  NULL;


  //! \todo initialize KdTree

  float delta_y = 1.0f / (img->height * 0.5f); //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step 
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) * aspect * scene->cam.ydir;

  float delta_x = 1.0f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x = (0.5f - img->width * 0.5f) / (img->width * 0.5f) *scene->cam.xdir;
    
  for(size_t j=0; j<img->height; j++) {
    if(j!=0) printf("\033[A\r");
    float progress = (float)j/img->height*100.f;
    printf("progress\t[");
    int cpt = 0;
    for(cpt = 0; cpt<progress; cpt+=5) printf(".");
    for(       ; cpt<100; cpt+=5) printf(" ");
    printf("]\n");
#pragma omp parallel for
    for(size_t i=0; i<img->width; i++) {
      color3 *ptr = getPixelPtr(img, i,j);
      vec3 const_arg = scene->cam.center + ray_delta_x + ray_delta_y + float(i)*dx +float(j)*dy;
      for (int y = 0 ; y < 2 ; ++y)
      {
        for (int x = 0 ; x < 2 ; ++x)
        {
          vec3 ray_dir = const_arg+ float(x)*(dx/2.0f) +  float(y)*(dy/2.0f) ;
          Ray rx;
          rayInit(&rx, scene->cam.position, normalize(ray_dir));
          
          *ptr += trace_ray(scene, &rx, tree);
        }
      }

      *ptr = *ptr / 4.0f;
    }
  }
}