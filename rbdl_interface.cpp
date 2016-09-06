#include "rbdl_interface.hpp"

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

using std::cout;
using std::endl;

namespace frapu
{
RBDLInterface::RBDLInterface():
    model_(new Model),
    q_(),
    qdot_(),
    tau_(),
    viscous_(),
    lower_position_constraints_(),
    upper_position_constraints_()
{

}

bool RBDLInterface::load_model(const std::string& model_file)
{
    const char* c = model_file.c_str();
    if (!shared::URDFReadFromFile(c, model_, false)) {
        std::cerr << "Error loading model " << model_file << std::endl;
        abort();
    }

    for (std::map<std::string, unsigned int>::iterator it = model_->mBodyNameMap.begin(); it != model_->mBodyNameMap.end(); ++it) {
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

void RBDLInterface::setViscous(std::vector<double>& viscous)
{
    viscous_.clear();
    for (size_t i = 0; i < viscous.size(); i++) {
        viscous_.push_back(viscous[i]);
    }
}

void RBDLInterface::setGravity(double& gravity)
{
    model_->gravity = Vector3d(0.0, 0.0, -gravity);
}

void RBDLInterface::setPositionConstraints(std::vector<double>& lowerPositionConstraints,
        std::vector<double>& upperPositionConstraints)
{
    for (size_t i = 0; i < lower_position_constraints_.size(); i++) {
        lower_position_constraints_[i] = lowerPositionConstraints[i];
        upper_position_constraints_[i] = upperPositionConstraints[i];
    }
}

void RBDLInterface::getEEJacobian(std::vector<double>& q,
                                  Eigen::MatrixXd& res)
{
    q_ = VectorNd::Zero(model_->dof_count);
    for (size_t i = 0; i < q.size(); i++) {
        q_[i] = q[i];
    }
    unsigned int id = 0;
    for (std::map<std::string, unsigned int>::iterator it = model_->mBodyNameMap.begin(); it != model_->mBodyNameMap.end(); ++it) {
        if (it->first == "end_effector") {
            id = it->second;
        }
    }

    CalcBodySpatialJacobian(*model_, q_, id, res);
    res.row(0).swap(res.row(3));
    res.row(1).swap(res.row(4));
    res.row(2).swap(res.row(5));

}

bool RBDLInterface::contact_impulse(std::vector<double>& q,
                                    std::vector<double>& qdot,
                                    std::vector<double>& tau,
                                    std::string& body_name,
                                    std::vector<double>& body_point,
                                    std::vector<double>& world_normal,
                                    Eigen::VectorXd& res)
{
    unsigned int body_id = model_->mBodyNameMap[body_name];
    Vector3d bpoint(body_point[0], body_point[1], body_point[2]);
    Vector3d wnormal(world_normal[0], world_normal[1], world_normal[2]);
    ConstraintSet constraint_set;
    constraint_set.AddConstraint(body_id,
                                 bpoint,
                                 Vector3d(0, 0, 1));

    constraint_set.Bind(*model_);
    q_ = VectorNd::Zero(model_->dof_count);
    qdot_ = VectorNd::Zero(model_->dof_count);
    for (size_t i = 0; i < q.size(); i++) {
        q_[i] = q[i];
        qdot_[i] = qdot[i];
    }

    ComputeContactImpulsesDirect(*model_, q_, qdot_, constraint_set, res);

    return true;
}

bool RBDLInterface::forward_dynamics_constraints(std::vector<double>& q,
        std::vector<double>& qdot,
        std::vector<double>& tau,
        std::string& body_name,
        std::vector<double>& body_point,
        std::vector<double>& world_normal,
        Eigen::VectorXd& qddot)
{
    unsigned int body_id = model_->mBodyNameMap[body_name];
    Vector3d bpoint(body_point[0], body_point[1], body_point[2]);
    Vector3d wnormal(world_normal[0], world_normal[1], world_normal[2]);
    ConstraintSet constraint_set;
    constraint_set.AddConstraint(body_id,
                                 bpoint,
                                 wnormal);
    constraint_set.Bind(*model_);
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

    ForwardDynamicsContactsDirect(*model_, q_, qdot_, tau_, constraint_set, qddot);

    return true;



}

bool RBDLInterface::forward_dynamics(std::vector<double>& q,
                                     std::vector<double>& qdot,
                                     std::vector<double>& tau,
                                     Eigen::VectorXd& qddot)
{
    q_ = VectorNd::Zero(model_->dof_count);
    qdot_ = VectorNd::Zero(model_->dof_count);
    tau_ = VectorNd::Zero(model_->dof_count);
    for (size_t i = 0; i < q.size(); i++) {
        q_[i] = q[i];
        qdot_[i] = qdot[i];
        tau_[i] = tau[i] - viscous_[i] * qdot[i];;
    }

    ForwardDynamics(*model_, q_, qdot_, tau_, qddot);

    return true;
}

}
