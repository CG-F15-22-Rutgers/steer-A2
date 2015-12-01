/*!
*
* \author VaHiD AzIzI
*
*/


#include "obstacles/GJK_EPA.h"
#include <iostream>
using namespace std;
SteerLib::GJK_EPA::GJK_EPA()
{		
}


bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{
	std::vector<Util::Vector> A_minus_B_unpruned = A_min_B(_shapeA, _shapeB);
	std::vector<Util::Vector> A_minus_B = QuickHull(A_minus_B_unpruned);
	//cout << "All the Points" << endl;
	//for (int i = 0; i < A_minus_B.size(); i++)
	//{
	//	cout << A_minus_B[i].x << " " << A_minus_B[i].z<<endl;
	//}

	float error0 = 0.0001f;
	bool close_enough = false;
	Util::Vector origin;
	Util::Vector v = A_minus_B[rand() % (A_minus_B.size()-1)];
	simplex.push_back(v);

	while(!close_enough && v != origin) 
	{
		
		auto temp = projection(A_minus_B, v); 
		Util::Vector close_to_origin = std::get<1>(temp[0]);// get point in A_minus_B that has the shortest projection on v
		

		if(simplex.size() == 1)
		{
			
			if (close_to_origin == v)
				return false;
			else
			{				
				simplex.push_back(close_to_origin);
				//std::cout << "simplex 1st pt" <<" " << simplex[0].x << " " << simplex[0].z << std::endl;
				v = std::get<0>(ortho_projection_point(v, close_to_origin, origin));
			}
			//cout << "size = 1" << endl;
		}
		else if(simplex.size() == 2)
		{

			//int i = 0;
			//while(std::find(simplex.begin(), simplex.end(), std::get<1>(temp[i])) != simplex.end())
			//	i++;
			//simplex.push_back(std::get<1>(temp[i]));
			//std::tuple<Util::Vector, float, Util::Vector> temp2;
			//temp2 = shortest_projection2(simplex, origin);
			//v = std::get<0>(temp2);

			if (close_to_origin == simplex[0] || close_to_origin == simplex[1])
				return false;
			else
			{
				simplex.push_back(close_to_origin);
				//std::cout << "simplex 2nd pt" << " " << simplex[1].x << " " << simplex[1].z << std::endl;
			}

			std::tuple<Util::Vector, float, Util::Vector> temp2;
			temp2 = shortest_projection2(simplex, origin);
			v = std::get<0>(temp2);
			//cout << "size = 2" << endl;
		}
		else if (simplex.size() > 2)
		{
			//std::cout << "simplex:" << std::endl;
			//for (auto i = simplex.begin(); i != simplex.end(); i++)
			//{
			//	std::cout << (*i).x << "  "  << (*i).z << std::endl;
			//}

			float u = close_to_origin*v/v.length();
			
			close_enough = std::abs((std::abs(u)-v.length())) < error0 ;
			if (triangle_contain_point_2(simplex, origin))//check if convex hull of W contain origin
			{
				//cout << "oky" << endl;
				return true;
			}
			if(!close_enough)
			{
				//std::cout << "close_to_origin:" << close_to_origin.x << " " << close_to_origin.z << std::endl;
				//std::cout << "close_to_origin projection length:" << u << std::endl;
				//std::cout << "vector v:" << v.x << " " << v.z << " " << v.length() << std::endl;

				//
				//std::cout << "all points:" << std::endl;
				//for (auto i = A_minus_B.begin(); i != A_minus_B.end(); i++)
				//{
				//std::cout << (*i).x << "   "  << (*i).z<<std::endl;
				//}

				std::tuple<Util::Vector, float, Util::Vector> temp2 = shortest_projection2(simplex, origin);

				simplex.erase(std::remove(simplex.begin(), simplex.end(), std::get<2>(temp2)), simplex.end()); // update simplex as the closest triangle to origin
				simplex.push_back(close_to_origin);
				temp2 = shortest_projection2(simplex, origin);// update v as the origin's shortest projection point on all W's edges that is also inside W
				v = std::get<0>(temp2);
				//std::cout << "loop " << std::endl;
			}
			//cout << "size = 3" << endl;
		}
	}
	return  false;
}

std::vector<std::pair<float,Util::Vector>> SteerLib::GJK_EPA::projection(std::vector<Util::Vector>& W, Util::Vector& v)
{
	std::vector<std::pair<float,Util::Vector>> temp;
	for (auto i= W.begin(); i !=W.end(); i++)
	{
		temp.push_back(std::make_pair((*i)*v,*i));
	}
	std::sort(temp.begin(), temp.end(), [](auto &a, auto &b) {return std::get<0>(a) < std::get<0>(b); });
	return temp;
}


std::tuple<Util::Vector,float,Util::Vector> SteerLib::GJK_EPA::shortest_projection2(std::vector<Util::Vector>& W, Util::Vector& v)
{
	Util::Vector a,b,c, origin;
	a = W[0];
	b = W[1];
	c = W[2];
	std::vector<std::tuple<Util::Vector,float,Util::Vector>> temp; // 1st element is the projection point, 2nd is orthogonal distance, 3rd is the point to remove
	temp.push_back(std::make_tuple(ortho_projection_point(a, b, origin).first, ortho_projection_point(a, b, origin).second,c));
	temp.push_back(std::make_tuple(ortho_projection_point(b, c, origin).first, ortho_projection_point(b, c, origin).second,a));
	temp.push_back(std::make_tuple(ortho_projection_point(c, a, origin).first, ortho_projection_point(c, a, origin).second,b));
	temp.push_back(std::make_tuple(a,dist(a, origin), b));
	temp.push_back(std::make_tuple(b,dist(b, origin), c));
	temp.push_back(std::make_tuple(c,dist(c, origin), a));
	std::sort(temp.begin(), temp.end(), [](auto &a, auto &b) {return std::get<1>(a) < std::get<1>(b); });
	
	int n=0;
	/*for (int i = 0; i < temp.size(); i++)
	{
		cout << get<0>(temp[i]).x << " " << get<0>(temp[i]).z << " length:" << get<1>(temp[i]) << endl;
	}*/

	while (!triangle_contain_point_2(W, std::get<0>(temp[n])))
	{
		//std::cout << "point not in the triangle:"<<std::endl;
		//std::cout << get<0>(temp[n]).x << " " << get<0>(temp[n]).z << endl;
		n++;
		//getchar();
	}


	return temp[n];
	
}

bool SteerLib::GJK_EPA::triangle_contain_point(std::vector<Util::Vector>& W, Util::Vector& v)
{
	float err0 = 0.00000000001f;
	Util::Vector a, b, c;
	a = W[0];
	b = W[1];
	c = W[2];
	if ((sign(determinant(a, b, c)) == sign(determinant(a, b, v)) && sign(determinant(a, b, v)) == sign(determinant(b, c, v)) && sign(determinant(b, c, v)) == sign(determinant(c, a, v)) && sign(determinant(c, a, v)) == sign(determinant(a, b, c))) || std::abs(determinant(c, a, v)) < err0 || std::abs(determinant(a, b, v)) < err0 || std::abs(determinant(b, c, v)) < err0)
		return true;
	else
		return false;
}

bool SteerLib::GJK_EPA::triangle_contain_point_2(std::vector<Util::Vector>& W, Util::Vector& v)
{
	float y1 = W[0].z;
	float x1 = W[0].x;
	float y2 = W[1].z;
	float x2 = W[1].x;
	float y3 = W[2].z;
	float x3 = W[2].x;
	float x = v.x;
	float y = v.z;
	/*float d = ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
	float a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / d;
	float b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / d;
	float c = 1 - a - b;

	return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;*/

	float S = calculateS(W[0]-W[1],W[1]- W[2],W[2] - W[0]);
	float S1 = calculateS(W[0] - W[1], W[0] - v, W[1] - v);
	float S2 = calculateS(W[0] - W[2], W[0] - v, W[2] - v);
	float S3 = calculateS(W[1] - W[2], W[2] - v, W[1] - v);

	float error = 0.01;
	//cout << abs(S - S1 - S2 - S3) << endl;
	
	//cout << S << " " << S1 << " " << S2 << " " << S3 << endl;
	if (abs(S - S1 - S2 - S3) <= 0.01)
		return true;
	else return false;


}

float SteerLib::GJK_EPA::calculateS(Util::Vector A, Util::Vector B, Util::Vector C)
{
	float a = A.length();
	float b = B.length();
	float c = C.length();
	float s = 0.5*(a + b + c);
	float check = s*(s - a)*(s - b)*(s - c);
	if (check < 0)
		return 0;
	float area = sqrt(s*(s - a)*(s - b)*(s - c));
	return area;
}

int SteerLib::GJK_EPA::sign(float val)
{	
	float val0 = 0.0f;
	return ((val < val0) - (val0 < val));
}

float SteerLib::GJK_EPA::determinant(Util::Vector& a, Util::Vector& b, Util::Vector& c)
{
	float det = a.x*b.z + a.z*c.x + b.x*c.z - a.x*c.z - a.z*b.x - b.z*c.x;
	return det;
	//return negative if a,b,c form a right turn, and vice versa, 0 if linear
}

float SteerLib::GJK_EPA::dist(Util::Vector& a, Util::Vector& b)
{
	return(pow(a.x-b.x, 2)+pow(a.z-b.z,2));
}

std::pair<Util::Vector, float> SteerLib::GJK_EPA::ortho_projection_point(Util::Vector& a, Util::Vector& b, Util::Vector& c)
{
	//cout << "calculate ortho_projection_point" << endl;
	//cout << a.x << " " << a.z << " " << b.x << " " << b.z << endl;
	float t = - ((a.x-c.x)*(b.x-a.x) + (a.z-c.z)*(b.z-a.z))/((b.x-a.x)*(b.x-a.x) + (b.z-a.z)*(b.z-a.z));
	//cout << t << endl;
	Util::Vector u;
	u.x = t*b.x + (1.0f - t)*a.x;
	u.y = 0.0f;
	u.z= t*b.z + (1.0f - t)*a.z;
	//cout << u.x <<" "<< u.z<<endl;
	float len_square = pow(u.x, 2) + pow(u.z, 2);
	return std::make_pair(u, len_square);
}

std::tuple<Util::Vector, float, Util::Vector, Util::Vector> SteerLib::GJK_EPA::ortho_projection_point2(Util::Vector& a, Util::Vector& b, Util::Vector& c)
{
	float t = - ((a.x-c.x)*(b.x-a.x) + (a.z-c.z)*(b.z-a.z))/((b.x-a.x)*(b.x-a.x) + (b.z-a.z)*(b.z-a.z));
	Util::Vector u;
	u.x = t*b.x + (1.0f - t)*a.x-c.x;
	u.y = 0.0f;
	u.z = t*b.z + (1.0f - t)*a.z-c.z;
	float len_square = pow(u.x, 2) + pow(u.z, 2);
	return std::make_tuple(u, len_square, a, b);
}

std::pair<Util::Vector, float> SteerLib::GJK_EPA::ortho_projection_point3(Util::Vector& a, Util::Vector& b, Util::Vector& c)
{
	//cout << "calculate ortho_projection_point" << endl;
	//cout << a.x << " " << a.z << " " << b.x << " " << b.z << endl;
	float t = -((a.x - c.x)*(b.x - a.x) + (a.z - c.z)*(b.z - a.z)) / ((b.x - a.x)*(b.x - a.x) + (b.z - a.z)*(b.z - a.z));
	//cout << t << endl;
	Util::Vector u;
	u.x = t*b.x + (1.0f - t)*a.x;
	u.y = 0.0f;
	u.z = t*b.z + (1.0f - t)*a.z;
	//cout << u.x <<" "<< u.z<<endl;
	float len_square = pow(u.x-c.x, 2) + pow(u.z-c.z, 2);
	return std::make_pair(c, len_square);
}

std::vector<Util::Vector> SteerLib::GJK_EPA::A_min_B(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	std::vector<Util::Vector> A_B;
	Util::Vector a_b;
	for(auto i = _shapeA.begin(); i != _shapeA.end(); i++ )
		for(auto j = _shapeB.begin(); j != _shapeB.end(); j++)
		{
			
			a_b = (*i) - (*j);
			
			A_B.push_back(a_b);
		}
	return A_B;
}

std::vector<Util::Vector> SteerLib::GJK_EPA::QuickHull(std::vector<Util::Vector> W)
{
	std::vector<Util::Vector> CovxHull;
	std::vector<Util::Vector> Up, Low;
	float err = 0.00001f;
	// find a line seperate all pts into upper part and lower part
	std::sort(W.begin(), W.end(), [](auto &a, auto &b) {return a.x < b.x; });
	Util::Vector A = W[0];
	Util::Vector B = W[W.size() - 1];
	for (auto i = W.begin(); i != W.end(); i++)
	{
		if (determinant(A, B, (*i)) > 0 && std::abs(determinant(A, B, (*i))) > err )
			Up.push_back(*i);
		else if (determinant(A, B, (*i)) < 0 && std::abs(determinant(A, B, (*i))) > err)
			Low.push_back(*i);
	}
	//std::cout << "pt A" << "  " << A.x << "  " << A.z << std::endl;
	//std::cout << "pt B" << "  " << B.x << "  " << B.z << std::endl;
	//std::cout << "all pts in Low" << std::endl;
	//for (auto i = Low.begin(); i != Low.end(); i++)
	//{
	//	std::cout << (*i).x << " "  << (*i).z << std::endl;
	//}
	//std::cout << "all pts in Up" << std::endl;
	//for (auto i = Up.begin(); i != Up.end(); i++)
	//{
	//	std::cout << (*i).x << " " << (*i).z << std::endl;
	//}
	//std::cout << "line found" << std::endl;
	lower(A, B, Low, CovxHull);
	upper(A, B, Up, CovxHull);
	return CovxHull;
}


void SteerLib::GJK_EPA::lower(Util::Vector A, Util::Vector B, std::vector<Util::Vector> L, std::vector<Util::Vector>& CovxHull)
{
	if (L.size() == 0)
	{
		CovxHull.push_back(A);
	//	std::cout << "low finished" << std::endl;
	}
	else
	{
		float err = 0.000001f;
		std::vector<Util::Vector> left, right;
		Util::Vector lowest_pt;
		//find the lowest point
		std::vector<std::pair<Util::Vector, float>> temp;

		//std::cout << "pt A" << "  " <<   A.x << "  " << A.z << std::endl;
		//std::cout << "pt B" << "  " <<   B.x << "  " << B.z << std::endl;
		//std::cout << "all pts in L" << std::endl;
		for (auto i = L.begin(); i != L.end(); i++)
		{
			temp.push_back(ortho_projection_point3(A, B, (*i)));
		//	std::cout << (*i).x << " "  << (*i).z << std::endl;
		}
		std::sort(temp.begin(), temp.end(), [](auto &a, auto &b) {return std::get<1>(a) > std::get<1>(b); });
		lowest_pt = std::get<0>(temp[0]);

		//std::cout << "lowest pt" << " " << lowest_pt.x << " " <<  lowest_pt.z << std::endl;
		// divide and conquer
		for (auto i = L.begin(); i != L.end(); i++)
		{
			if (determinant(A, lowest_pt, (*i)) < 0 && std::abs(determinant(A, lowest_pt, (*i))) > err)
				left.push_back(*i);
			if (determinant(lowest_pt, B, (*i)) < 0 && std::abs(determinant(lowest_pt, B, (*i))) > err)
				right.push_back(*i);
		}
		//std::cout << "loop2" << std::endl;
		lower(A, lowest_pt, left, CovxHull);
		lower(lowest_pt, B, right, CovxHull);

	}
}

void SteerLib::GJK_EPA::upper(Util::Vector A, Util::Vector B, std::vector<Util::Vector> U, std::vector<Util::Vector>& CovxHull)
{
	if (U.size() == 0)
	{
		CovxHull.push_back(B);
		//std::cout << "up finished" << std::endl;
	}
	else
	{
		float err = 0.000001f;
		std::vector<Util::Vector> left, right;
		Util::Vector highest_pt;
		//find the highest point
		std::vector<std::pair<Util::Vector, float>> temp;
		for (auto i = U.begin(); i != U.end(); i++)
		{
			temp.push_back(ortho_projection_point3(A, B, (*i)));
		}
		std::sort(temp.begin(), temp.end(), [](auto &a, auto &b) {return std::get<1>(a) > std::get<1>(b); });
		highest_pt = std::get<0>(temp[0]);

		// divide and conquer
		for (auto i = U.begin(); i != U.end(); i++)
		{
			if (determinant(A, highest_pt, (*i)) > 0 && std::abs(determinant(A, highest_pt, (*i))) > err)
				left.push_back(*i);
			if (determinant(highest_pt, B, (*i)) > 0 && std::abs(determinant(highest_pt, B, (*i))) > err)
				right.push_back(*i);
		}
	//	std::cout << "loop1" << std::endl;
		upper(highest_pt, B, right, CovxHull);
		upper(A, highest_pt, left, CovxHull);
		
	}
}


std::pair<float,Util::Vector>  SteerLib::GJK_EPA::EPA(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{
	Util::Vector origin, a, b, c;
	std::vector<Util::Vector> A_minus_B_unpruned = A_min_B(_shapeA, _shapeB);
	std::vector<Util::Vector> A_minus_B = QuickHull(A_minus_B_unpruned);
	a = simplex[0];
	b = simplex[1];
	c = simplex[2];
//	cout << "all node:" << endl;
	for (int i = 0; i < A_minus_B.size(); i++)
	{
//		cout << A_minus_B[i].x << " " << A_minus_B[i].z << endl;
	}

	/*cout << "shapeA:" << endl;
	for (int i = 0; i < _shapeA.size(); i++)
	{
		cout << _shapeA[i].x << " " << _shapeA[i].z << endl;
	}
	cout << "shapeB:" << endl;
	for (int i = 0; i < _shapeB.size(); i++)
	{
		cout << _shapeB[i].x << " " << _shapeB[i].z << endl;
	}*/

	std::vector<std::tuple<Util::Vector, float, Util::Vector, Util::Vector>> queue;
	std::vector<std::tuple<Util::Vector, float, Util::Vector, Util::Vector>>::iterator min;

	queue.push_back(ortho_projection_point2(a, b, origin));
	queue.push_back(ortho_projection_point2(b, c, origin));
	queue.push_back(ortho_projection_point2(c, a, origin));

	//for (auto i = queue.begin(); i != queue.end(); i++)
	//{
	//	std::cout << std::get<1>(*i) << std::endl;
	//}
	//
	
	bool converge = false;
	//int loop = 0;
	while (1)
	{
		min = std::min_element(queue.begin(), queue.end(), [](auto &a, auto &b) {return std::get<1>(a) < std::get<1>(b); });
		//min = queue.begin();
		//for (auto it = queue.begin()+1; it != queue.end(); it++)
		//{
		//	if (get<1>(*min) > get<1>(*it))
		//	{
		//		min = it;
		//	}
		//}
		int i = min - queue.begin();

		Util::Vector left_node, right_node, dir;
		float len1 = std::get<1>(*min);
		dir = std::get<0>(*min);


		left_node = std::get<2>(*min);
		right_node = std::get<3>(*min);
		
		
		auto temp = projection(A_minus_B, dir);
	
		Util::Vector expand_node = std::get<1>(temp[A_minus_B.size()-1]);
		//cout << "queue node" << endl;
		//for (int i = 0; i < queue.size(); i++)
		//{
		//	cout << get<0>(queue[i]).x << " " << get<0>(queue[i]).z << " " << get<1>(queue[i]) << endl;
		//}


		//cout << "all node:" << endl;
		//for (int i = 0; i < A_minus_B.size(); i++)
		//{
		//	cout << A_minus_B[i].x << " " << A_minus_B[i].z << endl;
		//}
		//cout << "expand node: " << expand_node.x<<" "<<expand_node.y<< endl;

		//std::cout << "expand node" << std::endl;
		//std::cout << expand_node.x << "   " << expand_node.z << std::endl;

		//std::cout << "right node" << std::endl;
		//std::cout << right_node.x << "   "  << right_node.z << std::endl;

		//std::cout << "left node" << std::endl;
	 //   std::cout << left_node.x << "    " << left_node.z << std::endl;
	

		if (expand_node == left_node || expand_node == right_node)
		{	
			//std::cout << "EPA vector" << std::endl;
			//std::cout << dir.x << "    " << dir.z << std::endl;
			//std::cout << "depth" << std::endl;
			//std::cout << std::sqrt(len1) << std::endl;

			
			return std::pair<float, Util::Vector>(std::sqrt(len1), dir);
	//		std::cout << "len1 is" <<'  ' << std::sqrt(len1) << std::endl;
		}
		else
		{	
	//		float len2 = pow(std::get<0>(temp[0]), 2);
	//		std::cout <<"len2 is"<<  ' ' << std::sqrt(len2) << std::endl;
	//		converge = ((len2 - len1) < 0.0000001f);
	//		std::cout << len2 - len1 << std::endl;
	/*		if(converge)
			{
			
				return std::pair<float,Util::Vector>( std::sqrt(len2), dir);
			}
			else
			*/
				
			queue.push_back(ortho_projection_point2(left_node, expand_node, origin));
			queue.push_back(ortho_projection_point2(right_node, expand_node, origin));
			queue.erase(queue.begin()+i);
			
		}
	}
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	bool is_colliding;
	std::vector<Util::Vector> simplex;
	//check if A and B are convex
	if (QuickHull(_shapeA).size() != _shapeA.size() || QuickHull(_shapeB).size() != _shapeB.size())
	{
		//std::cout << "not convex" << std::endl;
		std::vector<std::vector<Util::Vector>> triangle_list_A, triangle_list_B;
		triangle_list_A = triangulation(_shapeA);
		triangle_list_B = triangulation(_shapeB);
		int size_A = triangle_list_A.size();
		int size_B = triangle_list_B.size();
		std::cout << "concave polygon detected, A is decomposed into " << size_A << " parts and B is decomposed into " << size_B << " parts" <<  std::endl;
		//std::cout << "decomp finished" << std::endl;
		//std::cout << "A original" << std::endl;
		//for (auto i = _shapeA.begin(); i != _shapeA.end(); i++)
		//{
		//	std::cout << (*i).x << " " << (*i).z << std::endl;
		//}
		//std::cout << "A list" << std::endl;
		//for (auto i = triangle_list_A.begin(); i != triangle_list_A.end(); i++)
		//{
		//	std::cout << "one triangle" << std::endl;
		//	std::cout << (*i)[0].x << " " << (*i)[0].z << std::endl;
		//	std::cout << (*i)[1].x << " " << (*i)[1].z << std::endl;
		//	std::cout << (*i)[2].x << " " << (*i)[2].z << std::endl;
		//}
		//std::cout << "B list" << std::endl;
		//for (auto i = triangle_list_B.begin(); i != triangle_list_B.end(); i++)
		//{
		//	std::cout << "one triangle" << std::endl;
		//	std::cout << (*i)[0].x << " " << (*i)[0].z << std::endl;
		//	std::cout << (*i)[1].x << " " << (*i)[1].z << std::endl;
		//	std::cout << (*i)[2].x << " " << (*i)[2].z << std::endl;
		//}
		bool flag = false;
		for (auto i = triangle_list_A.begin(); i != triangle_list_A.end(); i++)
			for (auto j = triangle_list_B.begin(); j != triangle_list_B.end(); j++)
			{
				simplex.clear();
				is_colliding = GJK((*i), (*j), simplex);
				//std::cout << "GJK finished" << std::endl;
				if (is_colliding == true)
				{
					//std::cout << "triangle from A" << std::endl;
					//std::cout << (*i)[0].x << " " << (*i)[0].z << std::endl;
					//std::cout << (*i)[1].x << " " << (*i)[1].z << std::endl;
					//std::cout << (*i)[2].x << " " << (*i)[2].z << std::endl;
					//std::cout << "triangle from B" << std::endl;
					//std::cout << (*j)[0].x << " " << (*j)[0].z << std::endl;
					//std::cout << (*j)[1].x << " " << (*j)[1].z << std::endl;
					//std::cout << (*j)[2].x << " " << (*j)[2].z << std::endl;
					return_penetration_depth = std::get<0>(EPA((*i), (*j), simplex));
					return_penetration_vector = std::get<1>(EPA((*i), (*j), simplex));
					int num_i = i - triangle_list_A.begin() + 1;
					int num_j = j - triangle_list_B.begin() + 1;
					std::cout << "penetration_depth is" << return_penetration_depth << "and penetration_vector is (" << return_penetration_vector.x << ", 0.0f, " << return_penetration_vector.z <<") for part " << num_i << "in A and part "<< num_j << "in B" <<std::endl;
					flag = true;
				}
				if (flag == true)
				{
					return_penetration_depth = 0.0f;
					return_penetration_vector = Util::Vector(0.0f, 0.0f, 0.0f);
					return true;
				}
			}
		//std::cout << "no collision" << std::endl;
		return false;
	}
	else
		intersect_for_convex(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB);
 }	


bool SteerLib::GJK_EPA::intersect_for_convex(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	bool is_colliding;
	std::vector<Util::Vector> simplex;
	is_colliding = GJK(_shapeA, _shapeB, simplex);
	if (is_colliding == true)
	{
		return_penetration_depth = std::get<0>(EPA(_shapeA, _shapeB, simplex));
		return_penetration_vector = std::get<1>(EPA(_shapeA, _shapeB, simplex));
		//		std::cout << return_penetration_depth << "   " << return_penetration_vector.x << "  " << return_penetration_vector.z << std::endl;
		return true;
	}
	else
		return false; // There is no collision
}

void SteerLib::GJK_EPA::earclipping(std::vector<std::vector<Util::Vector>>& triangle_list, std::vector<Util::Vector> pts)
{
	std::vector<Util::Vector> ear, rest;
	rest = pts; 
	//std::cout << "pts" << std::endl;
	//for (auto i = pts.begin(); i != pts.end(); i++)
	//{
	//	std::cout << (*i).x << " " << (*i).z << std::endl;
	//}
	if (pts.size() == 3)
	{
		triangle_list.push_back(pts);
		//std::cout << "size == 3" << std::endl;
	}
	else
	{
		// the ear to clip
		std::vector<Util::Vector>::iterator min;
		min = std::min_element(pts.begin(), pts.end(), [](auto &a, auto &b) {return a.x < b.x; });
		int n = min - pts.begin();
		ear.resize(3);
		if (n > 0)
		{
			ear[0] = pts[n - 1];
			ear[1] = pts[n];
			ear[2] = pts[n + 1];
		}
		else
		{
			ear[0] = pts[pts.size()-1];
			ear[1] = pts[0];
			ear[2] = pts[1];
		}
		//std::cout << "int n" << n << std::endl;
		//for (auto i = ear.begin(); i != ear.end(); i++)
		//{
		//	std::cout << "ear" << " " << (*i).x << " " << (*i).z << std::endl;
		//}

		//find the point that has the largest projection in the direction towards the ear
		Util::Vector dir, pt;
		dir = std::get<0>(ortho_projection_point2(ear[0], ear[2], ear[1])); //the direction outwards the ear
		std::vector<std::pair<float, Util::Vector>> temp = projection(pts, dir);
		int i = 0;
		while (std::get<1>(temp[i]) == ear[0] || std::get<1>(temp[i]) == ear[1] || std::get<1>(temp[i]) == ear[2])
		{
			i++;
		}
		pt = std::get<1>(temp[i]);
		//std::cout << "pt" << " " << pt.x << " " << pt.z << std::endl;
		if (!triangle_contain_point_2(ear, pt))
		{
			//std::cout << "not contain" << std::endl;
			triangle_list.push_back(ear);
			rest.erase(rest.begin() + n);
			//std::cout << "rest" << std::endl;
			//for (auto i = rest.begin(); i != rest.end(); i++)
			//{
			//	std::cout << (*i).x << " " << (*i).z << std::endl;
			//}
			earclipping(triangle_list, rest);
		}
		else
		{
			//std::cout << "contain" << std::endl;
			std::vector<Util::Vector> rest1 = rest, rest2 = rest;
			std::vector<Util::Vector>::iterator it;
			it = std::find(pts.begin(), pts.end(), pt);
			int m = it - pts.begin();
			if (m > n)
			{
				rest1.erase(rest1.begin() + n + 1, rest1.begin() + m);
				rest2.erase(rest2.begin() + m + 1, rest2.end());
				rest2.erase(rest2.begin(), rest2.begin() + n);
			}
			else if (n > m)
			{
				rest1.erase(rest1.begin() + m + 1, rest1.begin() + n);
				rest2.erase(rest2.begin() + n + 1, rest2.end());
				rest2.erase(rest2.begin(), rest2.begin() + m);
			}
			//std::cout << "rest1" << std::endl;
			//for (auto i = rest1.begin(); i != rest1.end(); i++)
			//{
			//	std::cout << (*i).x << " " << (*i).z << std::endl;
			//}
			//std::cout << "rest2" << std::endl;
			//for (auto i = rest2.begin(); i != rest2.end(); i++)
			//{
			//	std::cout << (*i).x << " " << (*i).z << std::endl;
			//}
			earclipping(triangle_list, rest1);
			earclipping(triangle_list, rest2);
		}
	}
}

std::vector<std::vector<Util::Vector>> SteerLib::GJK_EPA::triangulation(const std::vector<Util::Vector> pts)
{
	std::vector<std::vector<Util::Vector>> triangle_list;
	earclipping(triangle_list, pts);
	return triangle_list;
}
