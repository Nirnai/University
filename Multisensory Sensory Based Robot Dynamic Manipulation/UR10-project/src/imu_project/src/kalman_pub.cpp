#include "ros/ros.h"
#include "std_msgs/String.h"
#include "kalman.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include "geometry_msgs/Vector3Stamped.h"
#include "geometry_msgs/PoseStamped.h"
#include <tf/transform_broadcaster.h>
#include "visualization_msgs/Marker.h"
#include <math.h>
#include <imu_project/DesiredState.h>
#include <std_msgs/Bool.h>
#include <sensor_msgs/Imu.h>

Eigen::VectorXd sensor_measurments_x(2);
Eigen::VectorXd sensor_measurments_y(2);
Eigen::VectorXd sensor_measurments_z(2);
Eigen::VectorXd ang_vel_x(1);
Eigen::VectorXd ang_vel_y(1);
Eigen::VectorXd ang_vel_z(1);

double ang_vel_old_y = 0.0;
double num_deriv_ang_vel_y;

tf::Quaternion imu_orientation;

int shake_counter = 0;

bool flag_init_pos_measurement = false;

void shakeDetect(const std_msgs::Bool &msg)
{
  if (shake_counter == 0) {
    shake_counter++;
  }
}


void angVelExtract(const geometry_msgs::Vector3Stamped &msg)
{
  ang_vel_x(0) = msg.vector.x;
  ang_vel_y(0) = msg.vector.y;
  ang_vel_z(0) = msg.vector.z;

  num_deriv_ang_vel_y = (ang_vel_y(0) - ang_vel_old_y)/(0.008);
}


void acc_extract(const geometry_msgs::Vector3Stamped &msg)
{
  //Accel Input
  sensor_measurments_x(1) = msg.vector.x;
  sensor_measurments_y(1) = msg.vector.y;
  sensor_measurments_z(1) = msg.vector.z;

}

void kinect_extract(const geometry_msgs::PointStamped &msg)
{
  //Position Input
  sensor_measurments_x(0) = msg.point.x;
  sensor_measurments_y(0) = msg.point.y;
  sensor_measurments_z(0) = msg.point.z;
}

void orientationExtract(const sensor_msgs::ImuConstPtr& msg) {
  tf::quaternionMsgToTF(msg->orientation, imu_orientation);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "kalman_pub");
  ros::NodeHandle nh;

  double dt = 1.0/125; // Time step
  ros::Rate loop_rate(125);

  ros::Publisher accel_pub = nh.advertise<geometry_msgs::Vector3Stamped>("/kalman/filtered_accel", 1000);
  ros::Publisher vel_pub = nh.advertise<geometry_msgs::Vector3Stamped>("/kalman/filtered_vel", 1000);
  ros::Publisher pos_pub = nh.advertise<geometry_msgs::Vector3Stamped>("/kalman/filtered_pos", 1000);

  ros::Publisher state_pub = nh.advertise<imu_project::DesiredState>("/robot_des_state", 1000);
  imu_project::DesiredState desired_state;

  ros::Subscriber sub_imu = nh.subscribe("/acc_base", 1000, &acc_extract);
  ros::Subscriber sub_angular = nh.subscribe("/ang_vel_base", 1000, &angVelExtract);
  ros::Subscriber sub_kinect = nh.subscribe("/kinect_transformed", 1000, &kinect_extract);
  ros::Subscriber sub_shake = nh.subscribe("/shaking_detected", 1000, &shakeDetect);
  ros::Subscriber sub_imu_orientation = nh.subscribe("/imu/data", 1000, &orientationExtract);

  //******************* LINEAR STATE ESTIMATIONS *********************** //
  int n = 4; // Number of states
  int m = 2; // Number of measurements

  Eigen::MatrixXd A(n, n); // System dynamics matrix
  Eigen::MatrixXd C(m, n); // Output matrix from state to measurement space
  Eigen::MatrixXd Q(n, n); // Process noise covariance
  Eigen::MatrixXd R(m, m); // Measurement noise covariance
  Eigen::MatrixXd P(n, n); // Estimate error covariance

  // Discrete constant acceleration motion
  A = Eigen::MatrixXd::Identity(n, n);
  A(0,1) = A(1,2) = dt;
  A(0,2) = dt*dt/2;

  C = Eigen::MatrixXd::Zero(m, n);
  C(0,0) = 1;
  C(1,2) = C(1,3) = 1;
  //C(0,0) = C(1,1) = C(2,2) = C(3,6) = C(4,7) = C(5,8) = 0.0004;


  Q = Eigen::MatrixXd::Zero(n, n); // Process noise covariance
  Q(0,0) = dt*dt*dt*dt/4;
  Q(0,1) = Q(1,0) = dt*dt*dt/2;
  Q(1,1) = dt*dt;
  Q(2,2) = 1;
  Q(0,2) = Q(2,0) = dt*dt/2;
  Q(2,1) = Q(1,2) = dt;
  Q = Q*0.0001;

  R = Eigen::MatrixXd::Identity(m, m); // Measurement noise covariance
  R = R*0.003;
  R(0,0) = 0.003;
  P = Eigen::MatrixXd::Identity(n, n); // Estimate error covariance
  P = P*0.1;

  std::cout << "A: \n" << A << std::endl;
  std::cout << "C: \n" << C << std::endl;
  std::cout << "Q: \n" << Q << std::endl;
  std::cout << "R: \n" << R << std::endl;
  std::cout << "P: \n" << P << std::endl;

  // Construct the filter
  KalmanFilter kf_x(dt, A, C, Q, R, P);
  KalmanFilter kf_y(dt, A, C, Q, R, P);
  KalmanFilter kf_z(dt, A, C, Q, R, P);

  Eigen::VectorXd x0(n);
  x0 << 0, 0 , 0, 0;
  kf_x.init(0.0,x0);
  kf_y.init(0.0,x0);
  kf_z.init(0.0,x0);

  Eigen::VectorXd x_hat_x(n);
  Eigen::VectorXd x_hat_y(n);
  Eigen::VectorXd x_hat_z(n);

  geometry_msgs::Vector3Stamped filtered_accel;
  filtered_accel.header.frame_id = "/robot_base";

  geometry_msgs::Vector3Stamped filtered_vel;
  filtered_vel.header.frame_id = "/robot_base";

  geometry_msgs::Vector3Stamped filtered_pos;
  filtered_pos.header.frame_id = "/robot_base";

  //******************* END INIT LINEAR STATE ESTIMATIONS *********************** //




  //******************* ANGULAR STATE ESTIMATIONS *********************** //
  int n_ang = 2; // Number of states
  int m_ang = 1; // Number of measurements

  Eigen::MatrixXd A_ang(n_ang, n_ang); // System dynamics matrix
  Eigen::MatrixXd C_ang(m_ang, n_ang); // Output matrix from state to measurement space
  Eigen::MatrixXd Q_ang(n_ang, n_ang); // Process noise covariance
  Eigen::MatrixXd R_ang(m_ang, m_ang); // Measurement noise covariance
  Eigen::MatrixXd P_ang(n_ang, n_ang); // Estimate error covariance

  // Discrete constant acceleration motion
  A_ang = Eigen::MatrixXd::Identity(n_ang, n_ang);
  A_ang(0,0) = A_ang(1,1) = 1;
  A_ang(0,1) = dt;

  C_ang = Eigen::MatrixXd::Zero(m_ang, n_ang);
  C_ang(0,0) = 1;

  Q_ang = Eigen::MatrixXd::Identity(n_ang, n_ang); // Process noise covariance
  Q_ang = 0.000001*Q_ang;
  Q_ang(0,0) = 0.001;
  Q_ang(1,1) = 0.2;

  R_ang = Eigen::MatrixXd::Identity(m_ang, m_ang); // Measurement noise covariance
  R_ang(0,0) = 0.03;

  P_ang = Eigen::MatrixXd::Identity(n_ang, n_ang); // Estimate error covariance
  P_ang = P_ang*0.1;

  std::cout << "Angular Estimation" << std::endl;
  std::cout << "A_ang: \n" << A_ang << std::endl;
  std::cout << "C_ang: \n" << C_ang << std::endl;
  std::cout << "Q_ang: \n" << Q_ang << std::endl;
  std::cout << "R_ang: \n" << R_ang << std::endl;
  std::cout << "P_ang: \n" << P_ang << std::endl;

  // Construct the filter
  KalmanFilter kf_x_ang(dt, A_ang, C_ang, Q_ang, R_ang, P_ang);
  KalmanFilter kf_y_ang(dt, A_ang, C_ang, Q_ang, R_ang, P_ang);
  KalmanFilter kf_z_ang(dt, A_ang, C_ang, Q_ang, R_ang, P_ang);

  Eigen::VectorXd x0_ang(n_ang);
  x0_ang << 0, 0;
  kf_x_ang.init(0.0,x0_ang);
  kf_y_ang.init(0.0,x0_ang);
  kf_z_ang.init(0.0,x0_ang);

  Eigen::VectorXd x_hat_x_ang(n_ang);
  Eigen::VectorXd x_hat_y_ang(n_ang);
  Eigen::VectorXd x_hat_z_ang(n_ang);

  geometry_msgs::Vector3Stamped filtered_accel_ang;
  filtered_accel.header.frame_id = "/robot_base";

  geometry_msgs::Vector3Stamped filtered_vel_ang;
  filtered_vel.header.frame_id = "/robot_base";


  //******************* END INIT ANGULAR STATE ESTIMATIONS *********************** //



  tf::TransformBroadcaster br;
  tf::Transform transform;
  tf::Quaternion q;

  q.setRPY(0,0,0);
  q.normalize();
  transform.setRotation(q);
  geometry_msgs::Point p;

  tf::Transform transform_original_kinect;
  transform_original_kinect.setRotation(q);

  geometry_msgs::Point initial_position;
  tf::Quaternion initial_orientation;


  tf::Transform transform_relative_eef;
  tf::Quaternion transform_relative_eef_quat;

  while (ros::ok())
  {

    // ************* Linear Estimations ************* //
    kf_x.update(sensor_measurments_x);
    x_hat_x = kf_x.state();
    kf_y.update(sensor_measurments_y);
    x_hat_y = kf_y.state();
    kf_z.update(sensor_measurments_z);
    x_hat_z = kf_z.state();

    filtered_accel.vector.x = x_hat_x(2);
    filtered_accel.vector.y = x_hat_y(2);
    filtered_accel.vector.z = x_hat_z(2);

    filtered_vel.vector.x = x_hat_x(1);
    filtered_vel.vector.y = x_hat_y(1);
    filtered_vel.vector.z = x_hat_z(1);

    filtered_pos.vector.x = x_hat_x(0);
    filtered_pos.vector.y = x_hat_y(0);
    filtered_pos.vector.z = x_hat_z(0);

    // ************* END ************* //

    // ************* Angular Estimations ************* //
    kf_x_ang.update(ang_vel_x);
    x_hat_x_ang = kf_x_ang.state();
    kf_y_ang.update(ang_vel_y);
    x_hat_y_ang = kf_y_ang.state();
    kf_z_ang.update(ang_vel_z);
    x_hat_z_ang = kf_z_ang.state();

    filtered_accel_ang.vector.x = x_hat_x_ang(1);
    filtered_accel_ang.vector.y = x_hat_y_ang(1);
    filtered_accel_ang.vector.z = x_hat_z_ang(1);

    filtered_vel_ang.vector.x = x_hat_x_ang(0);
    filtered_vel_ang.vector.y = x_hat_y_ang(0);
    filtered_vel_ang.vector.z = x_hat_z_ang(0);

    // ************* END ************* //

    filtered_accel.header.stamp = ros::Time::now();
    filtered_vel.header.stamp = ros::Time::now();
    filtered_pos.header.stamp = ros::Time::now();

    accel_pub.publish(filtered_accel);
    vel_pub.publish(filtered_vel);
    pos_pub.publish(filtered_pos);

    if (shake_counter == 1) {
      ROS_INFO_STREAM("Start Publishing Desired Position");
      ros::Duration(3).sleep();
      initial_position.x = filtered_pos.vector.x;
      initial_position.y = filtered_pos.vector.y;
      initial_position.z = filtered_pos.vector.z;
      initial_orientation = imu_orientation;
      shake_counter++;
      ROS_INFO_STREAM(initial_position);


    }
    else if (shake_counter > 1) {
      desired_state.pos.position.x = x_hat_x(0) - initial_position.x;
      desired_state.pos.position.y = x_hat_y(0) - initial_position.y;
      desired_state.pos.position.z = x_hat_z(0) - initial_position.z;

      transform_relative_eef_quat = imu_orientation*initial_orientation.inverse();
      desired_state.pos.orientation.x = transform_relative_eef_quat.getX();
      desired_state.pos.orientation.y = transform_relative_eef_quat.getY();
      desired_state.pos.orientation.z = transform_relative_eef_quat.getZ();
      desired_state.pos.orientation.w = transform_relative_eef_quat.getW();

      desired_state.vel.linear = filtered_vel.vector;
      desired_state.acc.linear = filtered_accel.vector;
      // desired_state.vel.angular = filtered_vel_ang.vector;
      // desired_state.acc.angular = filtered_accel_ang.vector;

      // desired_state.vel.linear.x = 0;
      // desired_state.vel.linear.y = 0;
      // desired_state.vel.linear.z = 0;
      desired_state.vel.angular.x = 0;
      desired_state.vel.angular.y = 0;
      desired_state.vel.angular.z = 0;
      // desired_state.acc.linear.x = 0;
      // desired_state.acc.linear.y = 0;
      // desired_state.acc.linear.z = 0;
      desired_state.acc.angular.x = 0;
      desired_state.acc.angular.y = 0;
      desired_state.acc.angular.z = 0;


      state_pub.publish(desired_state);

      transform_relative_eef.setOrigin(tf::Vector3(desired_state.pos.position.x,desired_state.pos.position.y,desired_state.pos.position.z));
      transform_relative_eef.setRotation(transform_relative_eef_quat);

      br.sendTransform(tf::StampedTransform(transform_relative_eef,ros::Time::now(),"/test_eef","/rel_imu"));
    }


    transform.setOrigin(tf::Vector3(filtered_pos.vector.x,filtered_pos.vector.y,filtered_pos.vector.z));
    br.sendTransform(tf::StampedTransform(transform,ros::Time::now(),"/robot_base","/imu_estimate"));

    transform_original_kinect.setOrigin(tf::Vector3(sensor_measurments_x(0),sensor_measurments_y(0),sensor_measurments_z(0)));
    br.sendTransform(tf::StampedTransform(transform_original_kinect,ros::Time::now(),"/robot_base","/kinect_data"));

    ros::spinOnce();
    loop_rate.sleep();
  }

  return 0;
}
