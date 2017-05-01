#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init() {
  //create a 4D state vector, we don't know yet the values of the x state
  x_= VectorXd(4);

  // the initial state covariance matrix P
  P_ = MatrixXd(4, 4);
  P_ <<   1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  //the initial transition matrix F_
  F_ = MatrixXd(4, 4);
  F_ <<   1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  //the process covariance matrix Q
  Q_ = MatrixXd(4, 4);

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht_laser_ = H_laser_.transpose();
  MatrixXd PHt = P_ * Ht_laser_;
  MatrixXd S = H_laser_ * PHt + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd z_pred = Tools::CalculateHprime(x_);
  VectorXd y = z - z_pred;

  while( y(1) < -M_PI)
  {
    y(1) += 2*M_PI;
  }

  while ( y(1) > M_PI)
  {
    y(1) -= 2*M_PI;
  }


  MatrixXd Hj_ = Tools::CalculateJacobian(x_);
  MatrixXd Hjt = Hj_.transpose();

  MatrixXd PHt = P_ * Hjt;
  MatrixXd S = Hj_ * PHt + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;
}
