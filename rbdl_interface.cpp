#include "rbdl_interface.hpp"

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

using std::cout;
using std::endl;

namespace shared {
RBDLInterface::RBDLInterface():
	model_(new Model),
	q_(),
	qdot_(),
	tau_(),
	viscous_(),
	lower_position_constraints_(),
	upper_position_constraints_(){
	
}

bool RBDLInterface::load_model(const std::string &model_file) {
	const char *c = model_file.c_str();    
	if (!Addons::URDFReadFromFile (c, model_, false)) {
		std::cerr << "Error loading model " << model_file << std::endl;
		abort();
	}
	cout << "loaded model with dof " << model_->dof_count << endl;
		
	for (std::map<std::string, unsigned int>::iterator it=model_->mBodyNameMap.begin(); it!=model_->mBodyNameMap.end(); ++it) {
	    std::cout << it->first << " => " << it->second << '\n';
	}	
	
	q_ = Eigen::VectorXd::Zero(model_->dof_count);
	qdot_ = Eigen::VectorXd::Zero(model_->dof_count);
	tau_ = Eigen::VectorXd::Zero(model_->dof_count);
	
	viscous_.clear();
	for (size_t i = 0; i < model_->dof_count; i++) {
		viscous_.push_back(0.0);
		lower_position_constraints_.push_back(-100000.0);
		upper_position_constraints_.push_back(100000.0);		
	}
}

void RBDLInterface::setViscous(std::vector<double> &viscous) {
	viscous_.clear();
	for (size_t i = 0; i < viscous.size(); i++) {
		viscous_[i] = viscous[i];
	}
}

void RBDLInterface::setGravity(double &gravity) {
	model_->gravity = Vector3d(0.0, 0.0, -gravity);
}

void RBDLInterface::setPositionConstraints(std::vector<double> &lowerPositionConstraints,
			                               std::vector<double> &upperPositionConstraints) {
	for (size_t i = 0; i < lower_position_constraints_.size(); i++) {
		lower_position_constraints_[i] = lowerPositionConstraints[i];
		upper_position_constraints_[i] = upperPositionConstraints[i];
	}
}

void RBDLInterface::getEEJacobian(std::vector<double> &q,
			                      Eigen::MatrixXd &res) {
	q_ = VectorNd::Zero(model_->dof_count);
	for (size_t i = 0; i < q.size(); i++) {
		q_[i] = q[i];
	}
	unsigned int id = 0;
	for (std::map<std::string, unsigned int>::iterator it=model_->mBodyNameMap.begin(); it!=model_->mBodyNameMap.end(); ++it) {
		if (it->first == "end_effector") {
			id = it->second;
		}		
	}
	
	CalcBodySpatialJacobian(*model_, q_, id, res);
	res.row(0).swap(res.row(3));
	res.row(1).swap(res.row(4));
	res.row(2).swap(res.row(5));
	
}

bool RBDLInterface::forward_dynamics(std::vector<double> &q, 
			                         std::vector<double> &qdot, 
			                         std::vector<double> &tau,
			                         Eigen::VectorXd &qddot) {
	for (size_t i = 0; i < q.size(); i++) {
		tau[i] -= viscous_[i] * qdot[i];
	}
	
	q_ = VectorNd::Zero(model_->dof_count);
	qdot_ = VectorNd::Zero(model_->dof_count);
	tau_ = VectorNd::Zero(model_->dof_count);
	for (size_t i = 0; i < q.size(); i++) {
		q_[i] = q[i];
		qdot_[i] = qdot[i];
		tau_[i] = tau[i];
	}
	
	ForwardDynamics(*model_, q_, qdot_, tau_, qddot);
	/**for (size_t i = 0; i < q.size(); i++) {
		if (i == 0) { 
		    cout << "lower_position_constraints_[i] " << lower_position_constraints_[i] << endl;
		    cout << "q(i) " << q[i] << endl;
		}
		
		if (q[i] < lower_position_constraints_[i]) {
			cout << "hello " << i << endl;
		}
		if (q[i] < lower_position_constraints_[i] + 0.0001 && qddot(i) < 0) {
			cout << "set zero!!!" << endl;
			qddot(i) = 0;
		}
		else if (q[i] > upper_position_constraints_[i] - 0.0001 && qddot(i) > 0) {
			cout << "set zero1!!!" << endl;
			qddot(i) = 0;
		}
	}*/
	return true;
}

}