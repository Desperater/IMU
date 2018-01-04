#include "imudata.h"

extern std::pair< C::iterator, C::iterator > r;

using namespace Eigen;
using namespace Sophus;
typedef Eigen::Matrix<double,9,9> Matrix9d;


class IMUPreintegrator{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  IMUPreintegrator();
  IMUPreintegrator(const IMUPreintegrator& pre);
  
  //reset to empty
  void reset();
  
  //incremental update 1、delta measurements, 2.jacobians 3. covariance matrix
  void update(const Eigen::Vector3d &omega,const Eigen::Vector3d& acc,const double &dt);
  
  //delta measurements position/velocity/rotation
  inline Eigen::Vector3d getDeltaP(void)const{
    return _delta_P;
  }
  inline Eigen::Vector3d getDeltaV(void)const{
    return _delta_V;
  }
  inline Eigen::Vector3d getDeltaR(void)const{
    return _delta_R;
  }
  
  //jacobian of delta measurements w.r.t bias of gyro/acc
  inline Eigen::Matrix3d getJPBiasg(void)const{
    return _J_P_Biasg;
  }
  inline Eigen::Matrix3d getJPBiasa(void)const{
    return _J_P_Biasa;
  }
  inline Eigen::Matrix3d getJVBiasg(void)const{
    return _J_V_Biasg;
  }
  inline Eigen::Matrix3d getJVBiasa(void)const{
    return _J_V_Biasa;
  }
  inline Eigen::Matrix3d getJRBiasg(void)const{
    return _J_R_Biasg;
  }
  //noise covariance propagation of delta measurements
  //note order here is rotation-velocity-position;
  inline Matrix9d getCOVPVPhi(void)const{
    return _cov_P_V_Phi;
  }
  inline double getDeltaTime()const{
    return _delta_time;
  }
  
  //skew symmetric matrix
  static Eigen::Matrix3d skew(const Eigen::Vector3d& v){
    return Sophus::SO3::hat(v);
  }
  //according to Rodrigues formula exponential map from vec3 to matrix3*3
  static Eigen::Matrix3d Expmap(const Eigen::Vector3d& v){
    return Sophus::SO3::exp(v).matrix();
  }
  static Eigen::Matrix3d JacobianR(const Eigen::Vector3d& w){
    Eigen::Matrix3d Jr = Eigen::Matrix3d::Identity();
    double theta = w.norm();
    double theta2 = theta*theta;
    double theta3 = theta2*theta;
    //very small angle
    if(theta<0.00001){
      return Jr;
    }else{
      Eigen::Vector3d k = w.normalized();
      Eigen::Matrix3d K = skew(k);
      //Jr = Eigen::Matrix3d::Indentity()-(1-cos(theta))/theta*K+(1-sin(theta)/theta)*K*K;
      Jr = Eigen::Matrix3d::Identity()- (1-cos(theta))/theta2*K+(theta-sin(theta))/theta3*K*K;
    }
    return Jr;
  }
  static Eigen::Matrix3d JacobianRInv(const Eigen::Vector3d& w){
    Eigen::Matrix3d Jrinv = Eigen::Matrix3d::Indentity();
    double theta = w.norm();
    double theta2= theta*theta;
    if(w<0.00001){
      return Jrinv;
    }else{
      Eigen::Vector3d k = w.normalized();
      Eigen::Matrix3d K = skew(k);
      Jrinv = Eigen::Matrix3d::Identity()+0.5*K
	      +(1.0-(1.0+cos(theta))*theta/(2.0*sin(theta)))*K*K;
    }
  }
  //left jacobian of SO3,J_l(x) = J_r(-x),because for SO3 R*R_t=I
  static Eigen::Matrix3d JacobianL(const Eigen::Vector3d& w){
    return JacobianR(-w);
  }
  static Eigen::Matrix3d JacobianLInv(const Eigen::Vector3d& w){
    return JacobianRInv(-w);
  }
  inline Eigen::Quaterniond normalizeRotationQ(const Eigen::Quaterniond& r){
    Eigen::Quaterniond _r(r);
    if(_r.w()<0)
      _r.coeffs()*=-1;
    return _r.normalized();
  }
  inline Eigen::Matrix3d normalizeRotationM(const Eigen::Matrix3d& R){
    Eigen::Quaterniond qr(R);
    return normalizeRotationQ(qr).toRotationMatrix();
  }
private:
  /*
   * Don't use pointer as a member
   * operator = is used in g2o,so don't override it
   */
  //delta measurements positions/velocity/rotation(matrix)
  Eigen::Vector3d _delta_P;	//P_k+1 = P_k + V_k*dt+0.5*R_k*a_k*dt*dt
  Eigen::Vector3d _delta_V;	//V_k+1 = V_k + R_k*a_k*dt
  Eigen::Vector3d _delta_R;	//R_k+1 = R_k + exp(w_k*dt)	note:Rwc,Rwc'=Rwc*[w_body]x
  
  //jacobian of delta measurements w.r.t bias of gyr/acc
  Eigen::Matrix3d _J_P_Biasg;		// position / gyro
  Eigen::Matrix3d _J_P_Biasa;		// position / acc
  Eigen::Matrix3d _J_V_Biasg;		// velocity / gyro
  Eigen::Matrix3d _J_V_Biasa;		// velocity / acc
  Eigen::Matrix3d _J_R_Biasg;		// Rotation / gyro
  
  //noise covariance propagation of delta measurements
  Matrix9d _cov_P_V_Phi;
  double _delta_time;
  
};
