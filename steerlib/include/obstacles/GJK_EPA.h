/*!
*
* \author VaHiD AzIzI
*
*/


#ifndef GJK_EPA_H_
#define GJK_EPA_H_


#include "util/Geometry.h"
#include <utility>
#include <cmath>	

#include <algorithm>
#include <vector>
#include <tuple>


namespace SteerLib
{

    class STEERLIB_API GJK_EPA
    {
        public:
			GJK_EPA();

            /*
             *
             *  DO NOT CHANGE THE FUNCTION DEFINITION FOR intersect()
             *  -----------------------------------------------------
             *  Return type: 
             *  If there is a collision detected by the GJK algorithm the return value of intersect is TRUE, else FALSE
             *  -----------------------------------------------------
             *  Arguments:
             *  float& return_penetration_depth
             *  If there is a collision, this variable needs to be updated to the penetration depth calculated by the EPA algorithm
             *  TODO: YOU NEED TO POPULATE THIS VARIABLE
             *  
             *  Util::Vector& return_penetration_vector
             *  If there is a collision, this variable needs to be updated to the penetration vector calculated by the EPA algorithm
             *  TODO: YOU NEED TO POPULATE THIS VARIABLE
             *  
             *  NOTE: The above arguments are passed by reference, thereby their value change is reflected in the calling function.
             *  THIS IS CALLED PASS BY REFERENCE.
             *  This is way for functions to "return" multiple values
             *  ie. 
             *  funct()
             *  {
             *      return (boolean, val)
             *  }
             *  
             *  can be implemented as:
             *  
             *  boolean funct(ref& value)
             *  {
             *      value = val;
             *  }
             *  
             *  const std::vector<Util::Vector>& _shapeA
             *  The first input shape passed by the collisionAI. Every point p(x,y,z) is a "Util::Vector p(x,y,z);"
             *  std::vector<point> gives the list of points in the polygon as "std::vector<Util::Vector>" 
             *  NOTE: ALL THE POINTS LIE ON THE XZ PLANE, ie. y=0 for all p. So all Polygon = { p1(x1,0,z2), p2(x2,0,z2), p3(x3,0,z3) ...}
             *  where Polygon is the std::vector and p1..pn are Util::Vector of (x,0,z)
             *  
             *  const std::vector<Util::Vector>& _shapeB
             *  The second input shape passed by the collisionAI. 
             *  ----------------------------------------------------
             *  Algorithms:
             *  Implement the GJK and EPA Algorithms in 2D since for all points (x,0,z), y=0.
             *  ----------------------------------------------------
             *  Pseudocode:
             *  GJK(Polygon A, Polygon B)
             *  {
             *      if (A collides with B)
             *          return (Simplex, true)
             *      else
             *          return (NULL, false)
             *  }
             *  
             *  EPA(Polygon A, Polygon B, Simplex S)
             *  {
             *      Calculate penetration
             *      return (penetration_depth, penetration_vector)
             *  }
             *  
             *  intersect(Polygon A, Polygon B)
             *  {
             *      (Simplex, is_colliding) = GJK(A, B)
             *      if ( is_colliding == true)
             *      {
             *          (penetration_depth, penetration_vector) = EPA(A, B, Simplex)
             *          return (true, penetration_depth, penetration_vector)
             *      }
             *      else
             *      {
             *          return (false, 0, NULL)    
             *      }    
             *  }
             *  ----------------------------------------------------
             *  Implementation:
             *  Write your own helper functions, as long as you do not change the function definition of intersect.
             *  ONLY MODIFY THIS CLASS.
             *  DO NOT LEAVE ANY DEBUG PRINTS IN YOUR CODE IN THIS CLASS.
             *  ----------------------------------------------------
             *  Testing:
             *  In order to test your code, execute the command, inside /build/bin
             *  ./steersim -testcase polygons_test -ai collisionAI
             *  
             *  The output you see, should be:
             *  
             *  $ ./steersim -testcase polygons_test -ai collisionAI
             *   
             *  loaded module collisionAI
             *  loaded module testCasePlayer
             *  Initializing...
             *  Preprocessing...
             *   Collision detected between polygon No.0 and No.1 with a penetration depth of 2 and penetration vector of (-0,0,-1)
             *  Simulation is running...
             *  ----------------------------------------------------
             *  Additional Notes:
             *  Type of vertices are Vector which you can check Vector class in util/Geometry.h to get familier with overloaded operators
             *  intersection will be called for all polygons, for all unique pairs, ie. if we have 3 polygons it is called 3 times
             *  one between polygon No. 1 and No. 2
             *  one between polygon No. 1 and No. 3
             *  one between polygon No. 2 and No. 3
             *  For grading your results will be compared by our results with testcase polygons1, so please run your code with
             *  original version of polygon1
             *  DO NOT MODIFY polygon1.xml
             */
            static bool intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB);

			static std::vector<Util::Vector> A_min_B(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB);
			/*static std::vector<std::pair<float,Util::Vector>> shortest_projection(std::vector<Util::Vector>& W, Util::Vector& v);*/
			static bool triangle_contain_point(std::vector<Util::Vector>& W, Util::Vector& v);
			static bool triangle_contain_point_2(std::vector<Util::Vector>& W, Util::Vector& v);
			static float calculateS(Util::Vector A, Util::Vector B, Util::Vector C);
			static int sign(float val);
			static float determinant(Util::Vector& a, Util::Vector& b, Util::Vector& c);
			static std::tuple<Util::Vector,float,Util::Vector> shortest_projection2(std::vector<Util::Vector>& W, Util::Vector& v);
			static float dist(Util::Vector& a, Util::Vector& b);
			static std::pair<Util::Vector, float> ortho_projection_point(Util::Vector& a, Util::Vector& b, Util::Vector& c);
			static std::tuple<Util::Vector, float, Util::Vector, Util::Vector> ortho_projection_point2(Util::Vector& a, Util::Vector& b, Util::Vector& c);
			static std::pair<Util::Vector, float> ortho_projection_point3(Util::Vector& a, Util::Vector& b, Util::Vector& c);
			static std::vector<std::pair<float,Util::Vector>> projection(std::vector<Util::Vector>& W, Util::Vector& v);
	
			static std::vector<Util::Vector> QuickHull(std::vector<Util::Vector> W);
			static void lower(Util::Vector A, Util::Vector B, std::vector<Util::Vector> L, std::vector<Util::Vector>& CovxHull);
			static void upper(Util::Vector A, Util::Vector B, std::vector<Util::Vector> U, std::vector<Util::Vector>& CovxHull);

			static void earclipping(std::vector<std::vector<Util::Vector>>& triangle_list, std::vector<Util::Vector> pts);
			static std::vector<std::vector<Util::Vector>> triangulation(std::vector<Util::Vector> pts);
			
			static bool GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex);
			static std::pair<float,Util::Vector> EPA(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex);
			static bool intersect_for_convex(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB);
        private:
			
			

    }; // class GJK_EPA

} // namespace SteerLib


#endif /* GJK_EPA_H_ */
