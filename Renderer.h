/**
@file Renderer.h
@author JOL
*/
#pragma once
#ifndef _RENDERER_H_
#define _RENDERER_H_

#include "Color.h"
#include "Image2D.h"
#include "Ray.h"

/// Namespace RayTracer
namespace rt {

  struct Background {
    virtual Color backgroundColor( const Ray& ray ) = 0;
  };
  
  struct MyBackground : public Background {
    Color backgroundColor( const Ray& ray ) {
      //Set les couleurs de fond
      Color white = Color( 1.0, 1.0, 1.0 );
      Color black = Color( 0.0, 0.0, 0.0 );
      Color blue = Color( 0.0, 0.0, 1.0 );

      Real z = ray.direction[2];
    
      if(z >= 0 && z <= 0.5)
        return z * 2 * blue + (1 - z * 2) * white;
      else if(z > 0.5)
        return (z - 0.5) * 2 * black + (1 - (z - 0.5) * 2) * blue;
      Color result = Color( 0.0, 0.0, 0.0 );
      Real x = -0.5f * ray.direction[ 0 ] / ray.direction[ 2 ];
      Real y = -0.5f * ray.direction[ 1 ] / ray.direction[ 2 ];
      Real d = sqrt( x*x + y*y );
      Real t = std::min( d, 30.0f ) / 30.0f;
      x -= floor( x );
      y -= floor( y );
      if ( ( ( x >= 0.5f ) && ( y >= 0.5f ) ) || ( ( x < 0.5f ) && ( y < 0.5f ) ) )
        result += (1.0f - t)*Color( 0.2f, 0.2f, 0.2f ) + t * Color( 1.0f, 1.0f, 1.0f );
      else
        result += (1.0f - t)*Color( 0.4f, 0.4f, 0.4f ) + t * Color( 1.0f, 1.0f, 1.0f );
      return result;
    }
  };

  inline void progressBar( std::ostream& output,
                           const double currentValue, const double maximumValue)
  {
    static const int PROGRESSBARWIDTH = 60;
    static int myProgressBarRotation = 0;
    static int myProgressBarCurrent = 0;
    // how wide you want the progress meter to be
    double fraction = currentValue /maximumValue;
    
    // part of the progressmeter that's already "full"
    int dotz = static_cast<int>(floor(fraction * PROGRESSBARWIDTH));
    if (dotz > PROGRESSBARWIDTH) dotz = PROGRESSBARWIDTH;
    
    // if the fullness hasn't changed skip display
    if (dotz == myProgressBarCurrent) return;
    myProgressBarCurrent = dotz;
    myProgressBarRotation++;
    
    // create the "meter"
    int ii=0;
    output << "[";
    // part  that's full already
    for ( ; ii < dotz;ii++) output<< "#";
    // remaining part (spaces)
    for ( ; ii < PROGRESSBARWIDTH;ii++) output<< " ";
    static const char* rotation_string = "|\\-/";
    myProgressBarRotation %= 4;
    output << "] " << rotation_string[myProgressBarRotation]
           << " " << (int)(fraction*100)<<"/100\r";
    output.flush();
  }
  
  /// This structure takes care of rendering a scene.
  struct Renderer {

    /// The scene to render
    Scene* ptrScene;
    /// The origin of the camera in space.
    Point3 myOrigin;
    /// (myOrigin, myOrigin+myDirUL) forms a ray going through the upper-left
    /// corner pixel of the viewport, i.e. pixel (0,0)
    Vector3 myDirUL;
    /// (myOrigin, myOrigin+myDirUR) forms a ray going through the upper-right
    /// corner pixel of the viewport, i.e. pixel (width,0)
    Vector3 myDirUR;
    /// (myOrigin, myOrigin+myDirLL) forms a ray going through the lower-left
    /// corner pixel of the viewport, i.e. pixel (0,height)
    Vector3 myDirLL;
    /// (myOrigin, myOrigin+myDirLR) forms a ray going through the lower-right
    /// corner pixel of the viewport, i.e. pixel (width,height)
    Vector3 myDirLR;
    /// Background color or scene
    Background* ptrBackground;

    int myWidth;
    int myHeight;

    Renderer() : ptrScene( 0 ) {
      ptrBackground = new MyBackground();
    }
    Renderer( Scene& scene ) : ptrScene( &scene ) {
      ptrBackground = new MyBackground();
    }

    Color background( const Ray& ray )
    {
      Color result = Color( 0.0, 0.0, 0.0 );
      for ( Light* light : ptrScene->myLights )
        {
          Real cos_a = light->direction( ray.origin ).dot( ray.direction );
          if ( cos_a > 0.99f )
            {
              Real a = acos( cos_a ) * 360.0 / M_PI / 8.0;
              a = std::max( 1.0f - a, 0.0f );
              result += light->color( ray.origin ) * a * a;
            }
        }
      if ( ptrBackground != 0 ) result += ptrBackground->backgroundColor( ray );
      return result;
    }

    void setScene( rt::Scene& aScene ) { ptrScene = &aScene; }
    
    void setViewBox( Point3 origin, 
                     Vector3 dirUL, Vector3 dirUR, Vector3 dirLL, Vector3 dirLR )
    {
      myOrigin = origin;
      myDirUL = dirUL;
      myDirUR = dirUR;
      myDirLL = dirLL;
      myDirLR = dirLR;
    }

    void setResolution( int width, int height )
    {
      myWidth  = width;
      myHeight = height;
    }


    /// The main rendering routine
    void render( Image2D<Color>& image, int max_depth )
    {
      std::cout << "Rendering into image ... might take a while." << std::endl;
      image = Image2D<Color>( myWidth, myHeight );
      for ( int y = 0; y < myHeight; ++y ) 
        {
          Real    ty   = (Real) y / (Real)(myHeight-1);
          progressBar( std::cout, ty, 1.0 );
          Vector3 dirL = (1.0f - ty) * myDirUL + ty * myDirLL;
          Vector3 dirR = (1.0f - ty) * myDirUR + ty * myDirLR;
          dirL        /= dirL.norm();
          dirR        /= dirR.norm();
          for ( int x = 0; x < myWidth; ++x ) 
            {
              Real    tx   = (Real) x / (Real)(myWidth-1);
              Vector3 dir  = (1.0f - tx) * dirL + tx * dirR;
              Ray eye_ray  = Ray( myOrigin, dir, max_depth );
              Color result = trace( eye_ray );
              image.at( x, y ) = result.clamp();
            }
        }
      std::cout << "Done." << std::endl;
    }

    Color shadow( const Ray& ray, Color light_color ) {
      Color result = light_color;
      Real epsilon = 0.0001f;
      Ray r = ray;
      GraphicalObject* obj;
      Point3 p;
      while(result.max() > 0.003f) {
        r.origin += r.direction * epsilon;
        if(ptrScene->rayIntersection(r, obj, p) < 0) {
          Material mat = obj->getMaterial(p);
          result = result * mat.coef_refraction * mat.diffuse;
        } else {
          break;
        }
      }
      return result;
    }

    Ray refractionRay( const Ray& aRay, const Point3& p, Vector3 N, const Material& m ) {
      Real tmp;
      Real r = m.in_refractive_index / m.out_refractive_index;
      Real c  = (-1.0f) * N.dot(aRay.direction);

      //When the ray is inside the object and go out
      if(aRay.direction.dot(N) <= 0 ) {
          r = 1.0f /r;
      }
      if(c>0)
          tmp = r*c - sqrt(1 - ( (r*r) * (1 - (c*c) )));
      else {
          tmp = r * c + sqrt(1 - ((r * r) * (1 - (c * c))));
      }

      Vector3 vRefrac = Vector3(r*aRay.direction + tmp * N);

      //Total reflexion
      if( 1 - ( (r*r) * (1 - (c*c) )) < 0) {
          vRefrac = aRay.direction - (2 * (aRay.direction.dot(N)) * N);
      }

      Ray rRefrac = Ray(p + vRefrac * 0.01f,vRefrac,aRay.depth -1);

      return rRefrac;
    }

    /// The rendering routine for one ray.
    /// @return the color for the given ray.
    Color trace( const Ray& ray )
    {
      assert( ptrScene != 0 );
      Color result = Color( 0.0, 0.0, 0.0 );
      GraphicalObject* obj_i = 0; // pointer to intersected object
      Point3           p_i;       // point of intersection

      // Look for intersection in this direction.
      Real ri = ptrScene->rayIntersection( ray, obj_i, p_i );
      // Nothing was intersected
      if ( ri >= 0.0f ) return background(ray);
      Material mat = obj_i->getMaterial(p_i);

      // Reflexion
      if(ray.depth > 0 && mat.coef_reflexion != 0) {
        Vector3 n_refl = obj_i->getNormal(p_i);
        Vector3 v_refl  = ray.direction - (2 * (ray.direction.dot(n_refl)) * n_refl);
        Ray ray_refl = Ray(p_i + v_refl  * 0.001f,v_refl ,ray.depth -1);
        Color c_refl = trace(ray_refl);
        result += c_refl * mat.specular * mat.coef_reflexion;
      }

      // Refraction
      if(ray.depth > 0 && mat.coef_refraction != 0) {
        Ray ray_refrac = refractionRay(ray,p_i,obj_i->getNormal(p_i),mat);
        Color c_refrac = trace(ray_refrac);
        result += c_refrac * mat.diffuse  * mat.coef_refraction;
      }
      
      result += illumination(ray, obj_i, p_i);
      return result;
    }

    /// Calcule l'illumination de l'objet \a obj au point \a p, sachant que l'observateur est le rayon \a ray.
    Color illumination( const Ray& ray, GraphicalObject* obj, Point3 p ){
      Material mat = obj->getMaterial(p);
      //set the base color with the ambient color
      Color result = Color(0,0,0);
      //Pour chaque source de lumiière
      for(Light* light : ptrScene->myLights) {
        Ray pointToLight = Ray(p,light->direction(p));
        Color colorShadow = shadow(pointToLight,light->color(p));
        Vector3 d = light->direction(p);
        Vector3 normLight = obj->getNormal(p);
        //On calcule la direction mirroir
        Vector3 mirrorD = ray.direction - (2 * (ray.direction.dot(normLight)) * normLight);
        Real cos_D = d.dot(normLight)/(d.norm()*normLight.norm());
        if (cos_D > 0.0) result += mat.diffuse * cos_D * light->color(p);
        Real cos_S = mirrorD.dot(d)/(mirrorD.norm()*d.norm());
        if(cos_S > 0) result += mat.specular * pow(cos_S, mat.shinyness) * light->color(p);
        result = result * colorShadow;
      }
      result += mat.ambient;
      return result;
    }

    /// Calcule le vecteur réfléchi à W selon la normale N.
    //Vector3 reflect( const Vector3& W, Vector3 N ) 
  };

  
} // namespace rt

#endif // #define _RENDERER_H_
