/*
 * Node responsible for publishing the tfs in the project
 */

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "sensor_msgs/Imu.h"
#include "geometry_msgs/QuaternionStamped.h"
#include "tf/transform_listener.h"
#include "tf/transform_broadcaster.h"
#include "tf_conversions/tf_eigen.h"
#include "visualization_msgs/Marker.h"
#include "visualization_msgs/MarkerArray.h"
#include <geometry_msgs/Vector3Stamped.h>


geometry_msgs::Quaternion imu_orientation;
tf::Quaternion calibration_orientation;

visualization_msgs::Marker arrow_x;
visualization_msgs::Marker arrow_y;
visualization_msgs::Marker arrow_z;


void poseToRviz(const sensor_msgs::ImuConstPtr& msg)
{
  imu_orientation = msg->orientation;
}

void updateArrows(const geometry_msgs::Vector3StampedConstPtr& msg)
{
  arrow_x.scale.x = msg->vector.x;
  arrow_y.scale.x = msg->vector.y;
  arrow_z.scale.x = msg->vector.z;
}

void setCalibrationOrientation(const geometry_msgs::QuaternionConstPtr &msg)
{
  calibration_orientation.setW(msg->w);
  calibration_orientation.setX(msg->x);
  calibration_orientation.setY(msg->y);
  calibration_orientation.setZ(msg->z);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "tf_imu_sender");
  ros::NodeHandle nh;

  ros::Rate loop_rate(200);

  ros::Subscriber sub = nh.subscribe("imu/data", 1000, &poseToRviz);
  ros::Subscriber sub_acc = nh.subscribe("/kalman/filtered_accel", 1000, &updateArrows);
  ros::Subscriber sub_calibrate = nh.subscribe("calibration_orientation", 1000, &setCalibrationOrientation);

  ros::Publisher vis_pub = nh.advertise<visualization_msgs::Marker>("visualization_marker",0);

  // Arrows representing the accleration
  arrow_x.header.frame_id = "/robot_base";
  arrow_x.header.stamp = ros::Time();
  arrow_x.id = 1;
  arrow_x.type = visualization_msgs::Marker::ARROW;
  arrow_x.action = visualization_msgs::Marker::ADD;
  arrow_x.scale.x = 1;
  arrow_x.scale.y = 0.05;
  arrow_x.scale.z = 0.05;
  arrow_x.pose.orientation.x = 0.0;
  arrow_x.pose.orientation.y = 0.0;
  arrow_x.pose.orientation.z = 0.0;
  arrow_x.pose.orientation.w = 1.0;
  arrow_x.pose.position.x = 0;
  arrow_x.pose.position.y = 0;
  arrow_x.pose.position.z = 0;
  arrow_x.color.a = 1.0;

  arrow_y.header.frame_id = "/robot_base";
  arrow_y.header.stamp = ros::Time();
  arrow_y.id = 2;
  arrow_y.type = visualization_msgs::Marker::ARROW;
  arrow_y.action = visualization_msgs::Marker::ADD;
  arrow_y.scale.x = 1;
  arrow_y.scale.y = 0.05;
  arrow_y.scale.z = 0.05;
  arrow_y.pose.orientation.x = 0.0;
  arrow_y.pose.orientation.y = 0;
  arrow_y.pose.orientation.z = 0.7071068;
  arrow_y.pose.orientation.w = 0.7071068;
  arrow_y.pose.position.x = 0;
  arrow_y.pose.position.y = 0;
  arrow_y.pose.position.z = 0;
  arrow_y.color.a = 1.0;

  arrow_z.header.frame_id = "/robot_base";
  arrow_z.header.stamp = ros::Time();
  arrow_z.id = 3;
  arrow_z.type = visualization_msgs::Marker::ARROW;
  arrow_z.action = visualization_msgs::Marker::ADD;
  arrow_z.scale.x = 1;
  arrow_z.scale.y = 0.05;
  arrow_z.scale.z = 0.05;
  arrow_z.pose.orientation.x = 0.0;
  arrow_z.pose.orientation.y = -0.7071068;
  arrow_z.pose.orientation.z = 0.0;
  arrow_z.pose.orientation.w = 0.7071068;
  arrow_z.pose.position.x = 0;
  arrow_z.pose.position.y = 0;
  arrow_z.pose.position.z = 0;
  arrow_z.color.a = 1.0;


  tf::TransformBroadcaster br;

  tf::Transform transform;
  transform.setOrigin(tf::Vector3(0,0,0));

  tf::Quaternion imu_orientation_tf;

  calibration_orientation.setX(0.620);
  calibration_orientation.setY(-0.036);
  calibration_orientation.setZ(0.679);
  calibration_orientation.setW(-0.391);
  calibration_orientation.normalize();

  tf::Transform transform2;

  tf::Transform transform_kinect_top;
  transform_kinect_top.setOrigin(tf::Vector3(0,0,0));
  tf::Quaternion kinect_top_quat;
  kinect_top_quat.setX(0.5);
  kinect_top_quat.setY(-0.5);
  kinect_top_quat.setZ(0.5);
  kinect_top_quat.setW(-0.5);
  transform_kinect_top.setRotation(kinect_top_quat);


  tf::Transform transform_robot_base;
  transform_robot_base.setOrigin(tf::Vector3(1,0,0));
  tf::Quaternion robot_base_quat;
  robot_base_quat.setX(0);
  robot_base_quat.setY(0.1736482);
  robot_base_quat.setZ(0);
  robot_base_quat.setW(0.9848078);

  transform_robot_base.setRotation(robot_base_quat.inverse());

  tf::Transform transform_robot_base_sim;
  transform_robot_base_sim.setOrigin(tf::Vector3(1,0,0));
  tf::Quaternion robot_base_quat_sim;
  robot_base_quat_sim.setX(0);
  robot_base_quat_sim.setY(0);
  robot_base_quat_sim.setZ(0);
  robot_base_quat_sim.setW(1);
  transform_robot_base_sim.setRotation(robot_base_quat_sim);

  tf::Transform transform_test_eef;
  transform_test_eef.setOrigin(tf::Vector3(0.5,0,1));
  tf::Quaternion transform_test_eef_quat;
  transform_test_eef_quat.setX(0);
  transform_test_eef_quat.setY(0);
  transform_test_eef_quat.setZ(0);
  transform_test_eef_quat.setW(1);
  transform_test_eef.setRotation(transform_test_eef_quat);

  while (ros::ok())
  {

    br.sendTransform(tf::StampedTransform(transform_robot_base_sim,ros::Time::now(),"/dh_arm_joint_0","/robot_base"));
    br.sendTransform(tf::StampedTransform(transform_robot_base,ros::Time::now(),"/robot_base","/kinect_top"));
    br.sendTransform(tf::StampedTransform(transform_kinect_top,ros::Time::now(),"/kinect_top","/kinect"));



    transform.setRotation(calibration_orientation.inverse());
    br.sendTransform(tf::StampedTransform(transform,ros::Time::now(),"/robot_base","/enu"));


    tf::quaternionMsgToTF(imu_orientation,imu_orientation_tf);
    imu_orientation_tf.normalize();

    transform2.setRotation(imu_orientation_tf);
    transform2.setOrigin(tf::Vector3(0,0,0));
    br.sendTransform(tf::StampedTransform(transform2,ros::Time::now(),"/enu","/imu"));

    br.sendTransform(tf::StampedTransform(transform_test_eef,ros::Time::now(),"/robot_base","/test_eef"));

    vis_pub.publish(arrow_x);
    vis_pub.publish(arrow_y);
    vis_pub.publish(arrow_z);

    ros::spinOnce();

    loop_rate.sleep();
  }

  return 0;
}
