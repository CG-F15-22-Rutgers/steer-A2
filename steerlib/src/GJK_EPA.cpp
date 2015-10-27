/*!
*
* \author VaHiD AzIzI
*
*/


#include "obstacles/GJK_EPA.h"


SteerLib::GJK_EPA::GJK_EPA()
{		
}


bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{
	std::vector<Util::Vector> A_minus_B = A_min_B(_shapeA, _shapeB);
	
	float error0 = 0.00001f;
	bool close_enough = false;
	Util::Vector origin;
	Util::Vector v = A_minus_B[rand() % (A_minus_B.size()-1)];
	simplex.push_back(v);
	
	while((!close_enough && v != origin) || simplex.size() < 3)
	{
		auto temp = shortest_projection(A_minus_B, v); // get point in A_minus_B that has the shortest projection on v
		
		if(simplex.size() == 1)
		{
			if (std::get<1>(temp[0]) == v)
			{
				simplex.push_back(std::get<1>(temp[1]));
				v = std::get<1>(temp[1]);
			}
			else
			{					
				simplex.push_back(std::get<1>(temp[0]));
				v = std::get<1>(temp[0]);
			}
		}
		else if(simplex.size() == 2)
		{
			int i = 0;
			while(std::find(simplex.begin(), simplex.end(), std::get<1>(temp[i])) != simplex.end())
				i++;
			simplex.push_back(std::get<1>(temp[i]));
			std::tuple<Util::Vector, float, Util::Vector> temp2;
			temp2 = shortest_projection2(simplex, origin);
			v = std::get<0>(temp2);
		}
		else if (simplex.size() > 2)
		{
			float u = std::get<1>(temp[0])*v/v.length();
			close_enough = std::abs((std::abs(u)-v.length())) < error0 ;
			if(triangle_contain_point(simplex, origin))//check if convex hull of W contain origin
				return true;
			if(!close_enough)
			{
				auto temp2 = shortest_projection2(simplex, origin);
				simplex.erase(std::remove(simplex.begin(), simplex.end(), std::get<2>(temp2)), simplex.end()); // upadte simplex as the closest triangle to origin
				simplex.push_back(std::get<1>(temp[0]));
				temp2 = shortest_projection2(simplex, origin);// update v as the origin's shortest projection point on all W's edges that is also inside W
				v = std::get<0>(temp2);
			}
		}
	}
	return  false;
}

std::vector<std::pair<float,Util::Vector>> SteerLib::GJK_EPA::shortest_projection(std::vector<Util::Vector>& W, Util::Vector& v)
{
	std::vector<std::pair<float,Util::Vector>> temp;
	for (auto i= W.begin(); i !=W.end(); i++)
	{
		temp.push_back(std::make_pair(std::abs((*i)*v),*i));
	}
	std::sort(temp.begin(), temp.end(), [](auto &a, auto &b) {return std::get<0>(a) < std::get<0>(b); });
	return temp;
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
	while (!triangle_contain_point(W, std::get<0>(temp[n])))
		n++;
	return temp[n];
	
}

bool SteerLib::GJK_EPA::triangle_contain_point(std::vector<Util::Vector>& W, Util::Vector& v)
{
	Util::Vector a,b,c;
	a = W[0];
	b = W[1];
	c = W[2];
	if ((sign(determinant(a,b,c)) ==  sign(determinant(a,b,v)) && sign(determinant(a,b,v)) == sign(determinant(b,c,v)) && sign(determinant(b,c,v)) == sign(determinant(c,a,v)) && sign(determinant(c,a,v)) == sign(determinant(a,b,c))) || sign(determinant(c,a,v)) == 0 || sign(determinant(a,b,v)) == 0 || sign(determinant(b,c,v)) == 0 )
		return true;
	else 
		return false;
}

int SteerLib::GJK_EPA::sign(float val)
{	
	float val0 = 0.0f;
	return ((val < val0) - (val0 < val));
}

float SteerLib::GJK_EPA::determinant(Util::Vector& a, Util::Vector& b, Util::Vector& c)
{
	return (a.x*b.z+a.z*c.x+b.x*c.z-a.x*c.z-a.z*b.x-b.z*c.x);
}

float SteerLib::GJK_EPA::dist(Util::Vector& a, Util::Vector& b)
{
	return(pow(a.x-b.x, 2)+pow(a.z-b.z,2));
}

std::pair<Util::Vector, float> SteerLib::GJK_EPA::ortho_projection_point(Util::Vector& a, Util::Vector& b, Util::Vector& c)
{
	float t = - ((a.x-c.x)*(b.x-a.x) + (a.z-c.z)*(b.z-a.z))/((b.x-a.x)*(b.x-a.x) + (b.z-a.z)*(b.z-a.z));
	Util::Vector u;
	u.x = (t*b.x + (1.0f - t))*a.x;
	u.y = 0.0f;
	u.z= t*b.z + (1.0f - t)*a.z;
	float len_square = pow(t*(b.x-a.x)+a.x-c.x, 2) + pow(t*(b.z-a.z)+a.z-c.z, 2);
	return std::make_pair(u, len_square);
}

std::tuple<Util::Vector, float, Util::Vector, Util::Vector> SteerLib::GJK_EPA::ortho_projection_point2(Util::Vector& a, Util::Vector& b, Util::Vector& c)
{
	float t = - ((a.x-c.x)*(b.x-a.x) + (a.z-c.z)*(b.z-a.z))/((b.x-a.x)*(b.x-a.x) + (b.z-a.z)*(b.z-a.z));
	Util::Vector u;
	u.x = (t*b.x + (1.0f - t))*a.x-c.x;
	u.y = 0.0f;
	u.z = t*b.z + (1.0f - t)*a.z-c.z;
	float len_square = pow(t*(b.x-a.x)+a.x-c.x, 2) + pow(t*(b.z-a.z)+a.z-c.z, 2);
	return std::make_tuple(u, len_square, a, b);
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


std::pair<float,Util::Vector>  SteerLib::GJK_EPA::EPA(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{
	Util::Vector origin, a, b, c;
	std::vector<Util::Vector> A_minus_B = A_min_B(_shapeA, _shapeB);
	a = simplex[0];
	b = simplex[1];
	c = simplex[2];
	std::vector<std::tuple<Util::Vector, float, Util::Vector, Util::Vector>> queue;
	std::vector<std::tuple<Util::Vector, float, Util::Vector, Util::Vector>>::iterator min;

	queue.push_back( ortho_projection_point2(a, b, origin));
	queue.push_back( ortho_projection_point2(b, c, origin));
	queue.push_back( ortho_projection_point2(c, a, origin));
	
	bool converge = false;
	while (!converge)
	{
		min = std::min_element(queue.begin(), queue.end(), [](auto &a, auto &b) {return std::get<1>(a) < std::get<1>(b); });
		Util::Vector left_node, right_node, dir;
		float len1 = std::get<1>(*min);
		dir= std::get<0>(*min);
		left_node = std::get<2>(*min);
		right_node = std::get<3>(*min);
		
		auto temp = projection(A_minus_B, dir);
		Util::Vector expand_node = std::get<1>(temp[0]);
		float len2 = pow(std::get<0>(temp[0])/dir.length(), 2);
		converge = (len2 - len1 < 0.00001f);
		if(converge)
		{
			return std::pair<float,Util::Vector>( std::sqrt(len2), dir);
		}
		else
		{
			queue.push_back( ortho_projection_point2(left_node, expand_node, origin));
			queue.push_back( ortho_projection_point2(right_node, expand_node, origin));
			queue.erase(min);
		}
	}
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	bool is_colliding;
	std::vector<Util::Vector> simplex;
	is_colliding = GJK(_shapeA, _shapeB, simplex);
	if(is_colliding == true)
	{
		std::make_pair(return_penetration_depth, return_penetration_vector) = EPA(_shapeA, _shapeB, simplex);
		return true;
	}
	else	
		return false; // There is no collision
}	





