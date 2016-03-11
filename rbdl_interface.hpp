#ifndef RBDL_INTERFACE_HPP_
#define RBDL_INTERFACE_HPP_
#include <iostream>
#include <rbdl/rbdl.h>
#include <rbdl/addons/urdfreader/urdfreader.h>
#include <Eigen/Dense>

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

namespace shared {
class RBDLInterface {
public:
	
	RBDLInterface();
	
	bool load_model(const std::string &model_file);
	
	bool forward_dynamics(std::vector<double> &q, 
			              std::vector<double> &qdot, 
			              std::vector<double> &tau,
			              Eigen::VectorXd &qddot);
	
	bool contact_impulse(std::vector<double> &q, 
                         std::vector<double> &qdot, 
                         std::vector<double> &tau,
                         std::string &body_name,
                         std::vector<double> &body_point,
                         std::vector<double> &world_normal,						              
                         Eigen::VectorXd &res);
	
	bool forward_dynamics_constraints(std::vector<double> &q, 
	                                  std::vector<double> &qdot, 
	                                  std::vector<double> &tau,
				                      std::string &body_name,
						              std::vector<double> &body_point,
						              std::vector<double> &world_normal,						              
	                                  Eigen::VectorXd &qddot);
	
	void setViscous(std::vector<double> &viscous);
	
	void setGravity(double &gravity);
	
	void setPositionConstraints(std::vector<double> &lowerPositionConstraints,
			                    std::vector<double> &upperPositionConstraints);
	
	void getEEJacobian(std::vector<double> &q,
			           Eigen::MatrixXd &res);
	
private:	
	RigidBodyDynamics::Model* model_;
	
	VectorNd q_;
    VectorNd qdot_;
	VectorNd tau_;
	
	std::vector<double> viscous_;
	
	std::vector<double> lower_position_constraints_;
	std::vector<double> upper_position_constraints_;
	
};

}

#endif