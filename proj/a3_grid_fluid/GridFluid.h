#ifndef __GridFluid_h__
#define __GridFluid_h__
#include "Common.h"
#include "Grid.h"
#include "Particles.h"

//////////////////////////////////////////////////////////////////////////
////Grid fluid simulator
template<int d> class GridFluid
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;
public:
	Grid<d> grid;
	Array<VectorD> u;		////velocity on grid nodes
	Array<real> div_u;		////velocity divergence on grid nodes (right hand side of the Poisson equation)
	Array<real> p;			////pressure
	Array<real> vor;		////vorticity

	int node_num=0;
	VectorD src_pos=VectorD::Ones()*(real).5;
	VectorD src_vel=VectorD::Unit(0)*(real)1.5;
	real src_radius=(real).1;


	virtual void Initialize()
	{
		int n=64;
		VectorDi cell_counts=VectorDi::Ones()*n;
		cell_counts[1]/=2;
		real dx=(real)2./(real)n;
		VectorD domain_min=VectorD::Zero();
		grid.Initialize(cell_counts,dx,domain_min);
		node_num=grid.node_counts.prod();

		u.resize(node_num,VectorD::Unit(0)*(real).01);
		div_u.resize(node_num,(real)0);
		p.resize(node_num,(real)0);
		vor.resize(node_num,(real)0);
	}
/////////////////////////////////////////////////////////////////////////////////////

	virtual void Advection(real dt)
	{
		////Advection using the semi-Lagrangian method
		Array<VectorD> u_copy=u;
		for(int i=0;i<node_num;i++){
			u[i]=VectorD::Zero();

			// Find midpoint
			VectorD result  = Pos(i) - u[i]*(dt/2);
			
			// Find the value at this point using interpolation function
			VectorD inter = Interpolate(u_copy,result);

			// Use the midpoint to find the starting point
			result = Pos(i) -inter*(dt/2);

			// Find the value at this point using interpolation function
			u[i] = Interpolate(u_copy, result);
		}
	}

	virtual void Projection()
	{
		real dx=grid.dx;
		real dx2=grid.dx*grid.dx;

		////Projection step 1: calculate the velocity divergence on each node
		////Read this sample code to learn how to access data with the node index and coordinate
		std::fill(div_u.begin(),div_u.end(),(real)0);
		for(int i=0;i<node_num;i++){
			if(Bnd(i))continue;		////ignore the nodes on the boundary
			VectorDi node=Coord(i);
			div_u[i]=(real)0;

			for(int j=0;j<d;j++){
				VectorD u_1=u[Idx(node-VectorDi::Unit(j))];
				VectorD u_2=u[Idx(node+VectorDi::Unit(j))];
				div_u[i]+=(u_2[j]-u_1[j])/(2*dx);}
		}

		////Projection step 2: solve the Poisson's equation -lap p= div u
		////using the Gauss-Seidel iterations
		std::fill(p.begin(),p.end(),(real)0);

		// 40 iterations
		for(int iter=0;iter<40;iter++){
			//for every coordinate
			for(int i=0;i<node_num;i++){
				if(Bnd(i))continue;		////ignore the nodes on the boundary
				VectorDi node=Coord(i);


				// Pressure value, start of the equation -( div * u)
				p[i] = -1 * div_u[i];

				// Laplacian Operator - Code doesnt rely on dimension 
				// for each dimension a, find a-1 and a+1 values, add together and multiply by 1/(dx^2) 
				for(int j=0;j<d;j++){
					// Vector has +-1 in the jth direction (dimension)
					real u_1=p[Idx(node-VectorDi::Unit(j))];
					real u_2=p[Idx(node+VectorDi::Unit(j))];
					p[i]+=((u_2+u_1)/(pow(dx,2)));}

				// Finish Laplacian by multiplying by (2 * dimension / dx^2) 
				p[i] /= (2*d)/(pow(dx,2));



				/*Your implementation ends*/
			}
		}

		////Projection step 3: correct velocity with the pressure gradient
		for(int i=0;i<node_num;i++){
			if(Bnd(i))continue;		////ignore boundary nodes
			VectorDi node=Coord(i);
			VectorD grad_p=VectorD::Zero();

			/*Your implementation starts*/
			for(int j=0;j<d;j++){
				real p_1=p[Idx(node-VectorDi::Unit(j))];
				real p_2=p[Idx(node+VectorDi::Unit(j))];
				grad_p[j] = (p_2 - p_1) * (1/(2*dx));
			}
			u[i] = u[i] - grad_p;
			/*Your implementation ends*/
		}
	}

	void Source()
	{
		for(int i=0;i<node_num;i++){
			VectorD pos=grid.Node(i);
			if((pos-src_pos).norm()<src_radius){u[i]=src_vel;}
		}
	}

	////this function works for 2D only
	void Vorticity_Confinement(const real dt)
	{
		real dx=grid.dx;

		////Vorticity confinement step 1: update vorticity
		std::fill(vor.begin(),vor.end(),(real)0);
		for(int i=0;i<node_num;i++){
			if(Bnd(i))continue;		////ignore boundary nodes
			VectorDi node=Coord(i);
			vor[i]=(real)0;

			//Curl - code changes for # of dimensions

			// x+1
			VectorD u_1=u[Idx(node+VectorDi::Unit(0))];
			// x-1
			VectorD u_2=u[Idx(node-VectorDi::Unit(0))];
			// y+1
			VectorD u_3=u[Idx(node+VectorDi::Unit(1))];
			// y-1
			VectorD u_4=u[Idx(node-VectorDi::Unit(1))];

			//Calculate Vorticity with curl approximation
			vor[i] = ((u_1[1] - u_2[1] - u_3[0] + u_4[0])/(2*dx));
		}

		////Vorticity confinement step 2: update N = (grad(|vor|)) / |grad(|vor|)|
		Array<VectorD> N(node_num,VectorD::Zero());
		for(int i=0;i<node_num;i++){
			if(Bnd(i))continue;		////ignore boundary nodes
			VectorDi node=Coord(i);
			N[i]=VectorD::Zero();
			
			// Calculate divergence of vorticity (same as projection step 1 code) 
			for(int j=0;j<d;j++){
				real vor_1=vor[Idx(node-VectorDi::Unit(j))];
				real vor_2=vor[Idx(node+VectorDi::Unit(j))];
				N[i][j] = (vor_2 - vor_1);
			}

			//Normalize 
			N[i].normalize();

		}

		////Vorticity confinement step 3: calculate confinement force and use it to update velocity
		real vor_conf_coef=(real)4;
		for(int i=0;i<node_num;i++){
			if(Bnd(i))continue;		////ignore boundary nodes
			VectorD f=vor_conf_coef*dx*Cross(N[i],vor[i]);
			u[i]+=f*dt;	////we don't have mass by assuming density=1
		}
	}



/////////////////////////////////////////////////////////////////////////////////////

	void Advance(const real dt)
	{
		Source();
		Advection(dt);
		Vorticity_Confinement(dt);
		Projection();
	}

	//////////////////////////////////////////////////////////////////////////
	////Helper functions
	////Please make sure to read these helper functions to understand how to access grid data
	//////////////////////////////////////////////////////////////////////////
	////return the node index given its coordinate
	int Idx(const Vector2i& node_coord) const
	{return grid.Node_Index(node_coord);}

	////return the coordinate given its index
	VectorDi Coord(const int node_index) const
	{return grid.Node_Coord(node_index);}

	////return the node position given its index
	VectorD Pos(const int node_index) const
	{return grid.Node(node_index);}

	////check if a node is on the boundary of the grid
	////given its coordinate or index
	bool Bnd(const Vector2i& node_coord) const
	{
		for(int i=0;i<d;i++){
			if(node_coord[i]==0||node_coord[i]==grid.node_counts[i]-1)
				return true;
		}
		return false;
	}
	bool Bnd(const int node_index) const
	{return Bnd(Coord(node_index));}

	Vector2 Cross(const Vector2& v,const real w) const
	{return Vector2(v[1]*w,-v[0]*w);}

	VectorD Interpolate(const Array<VectorD>& u,VectorD& pos)
	{
		////clamp pos to ensure it is always inside the grid
		real epsilon=grid.dx*(real)1e-3;
		for(int i=0;i<d;i++){
			if(pos[i]<=grid.domain_min[i])pos[i]=grid.domain_min[i]+epsilon;
			else if(pos[i]>=grid.domain_max[i])pos[i]=grid.domain_max[i]-epsilon;}

		////calculate the index, fraction, and interpolated values from the array
		VectorD cell_with_frac=(pos-grid.domain_min)/grid.dx;
		VectorDi cell=cell_with_frac.template cast<int>();
		VectorD frac=cell_with_frac-cell.template cast<real>();
		return Interpolate_Helper(cell,frac,u);
	}
	//////////////////////////////////////////////////////////////////////////

protected:
	////2D bi-linear interpolation
	Vector2 Interpolate_Helper(const Vector2i& cell,const Vector2& frac,const Array<Vector2>& u)
	{
		return ((real)1-frac[0])*((real)1-frac[1])*u[grid.Node_Index(cell)]
			+frac[0]*((real)1-frac[1])*u[grid.Node_Index(Vector2i(cell[0]+1,cell[1]))]
			+((real)1-frac[0])*frac[1]*u[grid.Node_Index(Vector2i(cell[0],cell[1]+1))]
			+frac[0]*frac[1]*u[grid.Node_Index(Vector2i(cell[0]+1,cell[1]+1))];
	}

	////Particles, for visualization only
public:
	Particles<d> particles;

	void Initialize_Visualization_Particles()
	{
		particles.Resize(1);
		particles.X(0)=VectorD::Zero();
	}

	void Update_Visualization_Particles(const real dt)
	{
		for(int i=0;i<particles.Size();i++){
			VectorD v=Interpolate(u,particles.X(i));
			particles.X(i)+=v*dt;
		}
	}
};

#endif