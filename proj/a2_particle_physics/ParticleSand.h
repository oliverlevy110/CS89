//#####################################################################
// Particle Sand (DEM)
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################

#ifndef __ParticleSand_h__
#define __ParticleSand_h__
#include "Common.h"
#include "Particles.h"
#include "ImplicitGeometry.h"

template<int d> class ParticleSand
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using MatrixD=Matrix<real,d>;
public:
	Particles<d> particles;
	real ks=(real)2e2;		////stiffness for the collision force
	real kd=(real).5e1;		////damping for the collision force
	VectorD g=VectorD::Unit(1)*(real)-1.;	////gravity

	Array<ImplicitGeometry<d>* > env_objects;	////list of implicit geometries describing the environment, by default it has one element, a circle with its normals pointing inward (Bowl)
	Array<Vector2i> particle_particle_collision_pairs;
	Array<Vector2i> particle_environment_collision_pairs;
	
	virtual void Advance(const real dt)
	{
		////Clear forces on particles
		for(int i=0;i<particles.Size();i++){
			particles.F(i)=VectorD::Zero();}

		////Accumulate body forces
		for(int i=0;i<particles.Size();i++){
			particles.F(i)+=particles.M(i)*g;}

		Particle_Environment_Collision_Detection();
		Particle_Environment_Collision_Response();
		Particle_Particle_Collision_Detection();
		Particle_Particle_Collision_Response();

		for(int i=0;i<particles.Size();i++){
			particles.V(i)+=particles.F(i)/particles.M(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): detect collision between particles and env_objects (implicit surface) and record the detected collisions in particle_environment_collision_pairs
	////env_objects is a list of implicit geometries, by default there is only one geometry (the bowl) in the list
	////Each element in particle_environment_collision_pairs is a Vector2i, with the first index for the particle and the second index for the env_objects
	virtual void Particle_Environment_Collision_Detection()
	{
		particle_environment_collision_pairs.clear();
		/* Your implementation start */
		/* Your implementation end */

		for(int i=0;i<particles.Size();i++){
			for(int j=0;j<env_objects.size();j++){
				if(env_objects[j]->Phi(particles.X(i)) - particles.R(i) < 0.0){
					Vector2i result (i,j);
					particle_environment_collision_pairs.push_back(result);
				}
			}
		}
	}
		
	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): compute the penalty-based collision force for the particles that are colliding with the env_objects
	////The collision response force consists of a spring force and a damping force
	virtual void Particle_Environment_Collision_Response()
	{
		for(int pair_i=0;pair_i<particle_environment_collision_pairs.size();pair_i++){
			int i=particle_environment_collision_pairs[pair_i][0];	////particle index
			int j=particle_environment_collision_pairs[pair_i][1];	////env_objects index
			VectorD collision_force=VectorD::Zero();

			VectorD env_vel=VectorD::Zero();

			double signed_dist = env_objects[j]->Phi(particles.X(i));
			VectorD surface_normal = env_objects[j]->Normal(particles.X(i));

			VectorD spring_force = ks*(signed_dist - particles.R(i))*(-1.0*surface_normal);
			VectorD dampening_force = kd*((env_vel - particles.V(i)).dot(-1.0*surface_normal))*(-1.0*surface_normal);
			
			collision_force = spring_force + dampening_force;

			particles.F(i)+=collision_force;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): find all the pairs of particles that are colliding each other and record the detected pairs in particle_particle_collision_pairs
	////Each element in particle_particle_collision_pairs is a Vector2i specifying the indices of the two colliding particles
	virtual void Particle_Particle_Collision_Detection()
	{
		particle_particle_collision_pairs.clear();
		/* Your implementation start */
		/* Your implementation end */
		for(int i=0;i<particles.Size();i++){
			for(int j=0;j<particles.Size();j++){
				if(i != j && ((particles.X(i)-particles.X(j)).norm())-(particles.R(i)+particles.R(j)) < 0.0){
					Vector2i result (i,j);
					particle_particle_collision_pairs.push_back(result);
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): compute penalty-based collision forces for pairs of colliding particles in particle_particle_collision_pairs and add the forces to particle.F 
	////The collision response force for each pair consists of a spring force and a damping force
	virtual void Particle_Particle_Collision_Response()
	{
		for(int pair_i=0;pair_i<particle_particle_collision_pairs.size();pair_i++){
			int i=particle_particle_collision_pairs[pair_i][0];	////the first particle index in the pair
			int j=particle_particle_collision_pairs[pair_i][1];	////the second particle index in the pair
			VectorD collision_force=VectorD::Zero();

			/* Your implementation start */
			/* Your implementation end */
			VectorD X_i=particles.X(i);
			VectorD X_j=particles.X(j);

			VectorD V_i=particles.V(i);
			VectorD V_j=particles.V(j);

			float l_o = particles.R(i) + particles.R(j);


			VectorD n = (X_j - X_i).normalized();

			VectorD fs_ij = (ks)*((X_j - X_i).norm()-(l_o)) * n;
			VectorD fd_ij = (kd)*((V_j - V_i).dot(n))*(n);
			VectorD f_ij = fs_ij + fd_ij;
			
			particles.F(i) = particles.F(i) + f_ij;
			particles.F(j)= particles.F(j) + (-1*f_ij);

		}

	}
};

#endif
