#ifndef IMU_DATA_H
#define IMU_DATA_H

#include <sophus/so3.h>
#include <sophus/se3.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;


class IMUData{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  IMUData(const double& gx,const double& gy,const double& gz,
	  const double& ax,const double& ay,const double& az,
	  const double& t):_g(gx,gy,gz),_a(ax,ay,az),_t(t){}
  
  //covariance of measurement
  static Matrix3d _gyrMeasCOV;
  static Matrix3d _accMeasCOV;
  
  static Matrix3d getGyrMeasCov(void){return _gyrMeasCOV;}
  static Matrix3d getAccMeasCov(void){return _accMeasCOV;}
  
  //covariance of bias random walk
  static Matrix3d _gyrBiasRWCOV;
  static Matrix3d _accBiasRWCOV;
  
  static Matrix3d getGyrBiasRWCov(void){return _gyrBiasRWCOV;}
  static Matrix3d getAccBiasRWCov(void){return _accBiasRWCOV;}
  
  //bias
  static double _gyrBiasRW2;
  static double _accBiasRW2;
  
  static double getGyrBiasRW2(void){return _gyrBiasRW2;}
  static double getAccBiasRW2(void){return _accBiasRW2;}
  ////raw data of imu
  Eigen::Vector3d _g;//陀螺仪数据
  Eigen::Vector3d _a;//加速度数据
  double _t ;	     //timestamp
};


#endif