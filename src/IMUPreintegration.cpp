#include "IMUPreintegration.h"

IMUPreintegrator::IMUPreintegrator()
{
  //delta measurements ,position velocity rotation(matrix)
  _delta_P.setZero();
  _delta_V.setZero();
  _delta_R.setIdentity();
  
  _J_P_Biasa.setZero();
  _J_P_Biasg.setZero();
  _J_V_Biasg.setZero();
  _J_V_Biasa.setZero();
  _J_R_Biasg.setZero();
  
  _cov_P_V_Phi.setZero();
  
  _delta_time = 0;
} 

IMUPreintegrator::IMUPreintegrator(const IMUPreintegrator& pre):
	_delta_P(pre._delta_P),_delta_V(pre._delta_V),_delta_R(pre._delta_R),
	_J_P_Biasg(pre._J_P_Biasg),_J_P_Biasa(pre._J_P_Biasa),_J_V_Biasg(_J_V_Biasg),_J_V_Biasa(pre._J_P_Biasa),_J_R_Biasg(pre._J_R_Biasg),
	_cov_P_V_Phi(pre._cov_P_V_Phi),_delta_time(pre._delta_time)
{  
}

void IMUPreintegrator::reset()
{
  //delta measurements, position /velocity/rotation(matrix)
  _delta_P.setZero();			//P_k+1 = P_k + V_k*dt + R_k+1*a_k*dt*dt/2
  _delta_V.setZero();			//V_k+1 = V_k + R_k*a_k*dt
  _delta_R.setIdentity();		//R_k+1 = R_k*exp(w_k*dt)
  
  //jacobian of delta measurements w.r.t bias of gyro/acc
  _J_P_Biasg.setZero();
  _J_P_Biasa.setZero();
  _J_V_Biasa.setZero();
  _J_V_Biasg.setZero();
  _J_R_Biasg.setZero();
  
  //noise covariance propagation of delta measurements
  _cov_P_V_Phi.setZero();
  
  _delta_time = 0;
  
}

//incrementally update 1.delta measurements 2.jacobians 3.covariance
//acc = acc_measurement -bias_a, last measurements not the current
//omega = gyro_measurement - bias_g
void IMUPreintegrator::update(const Vector3d& omega, const Vector3d& acc, const double& dt)
{
  double dt2 = dt*dt;
  
  Matrix3d dR = Expmap(omega*dt);
  Matrix3d Jr = JacobianR(omega*dt);
  
  // noise covariance propagation of delta measurements
  // err_k+1 = A*err_k + B*err_gyro + C*err_acc
  Eigen::Matrix3d I3x3 = Eigen::Matrix3d::Identity();
  Matrix9d A = Matrix9d::Identity();
  A.block<3,3>(6,6) = dR.transpose();
  A.block<3,3>(3,6) = -_delta_R*skew(acc)*dt2;
  A.block<3,3>(0,6) = -0.5*_delta_R*skew(acc)*dt2;
  A.blcok<3,3>(0,3) = I3x3*dt;
  Eigen::Matrix<double,9,3> Bg = Eigen::Matrix<double,9,3>::Zero();
  Bg.block<3,3>(6,0) =  Jr*dt;
  Eigen::Matrix<double,9,3> ca = Eigen::Matrix<double,9,3>::Zero();
  ca.block<3,3>(3,0) = _delta_R*dt;
  ca.block<3,3>(0,0) = 0.5*_delta_R*dt2;
  
  _cov_P_V_Phi = A*_cov_P_V_Phi*A.transpose()+
		 Bg*IMUData::getGyrMeasCov()*Bg.transpose()+
		 ca*IMUData::getAccMeasCov()*ca.transpose();
  //jacobian of delta measurements w.r.t bias of gyr/acc
  //update P first Then V,then R,because P's update need V&R's previous state
  _J_P_Biasg += _J_V_Biasg*dt - 0.5*_delta_R*skew(acc)*_J_P_Biasg*dt2;
  _J_P_Biasa += _J_V_Biasa*dt - 0.5*_delta_R*dt2;
  _J_V_Biasa += -_delta_R*dt;
  _J_V_Biasg += -_delta_R*skew(acc)*_J_R_Biasg*dt;
  _J_R_Biasg = dR.transpose()*_J_R_Biasg-Jr*dt;
  
  // delta measurements position/velocity/rotation
  // update P first then V then R;
  _delta_P += _delta_V*dt + 0.5*_delta_R*acc*dt2;
  _delta_V += _delta_R*acc*dt;
  _delta_R = normalizeRotationM(_delta_R*dR);
  
  _delta_time+=dt;
  
}
