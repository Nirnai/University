/*
 * Run this node when the IMU is placed on the kinect aligned with top
 * kinect coordinate systems. After calibration the IMU coordinate system should be aligned.
 *
 */

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "geometry_msgs/Quaternion.h"
#include "sensor_msgs/Imu.h"

geometry_msgs::Quaternion imu_orientation;

void poseToRviz(const sensor_msgs::ImuConstPtr& msg)
{
  imu_orientation = msg->orientation;
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "calibration");
  ros::NodeHandle nh;

  ros::Publisher calibration_pub = nh.advertise<geometry_msgs::Quaternion>("/calibration_orientation", 1000);

  ros::Subscriber sub = nh.subscribe("imu/data", 1000, &poseToRviz);

  ros::Rate loop_rate(10);

  while (ros::ok())
  {

    calibration_pub.publish(imu_orientation);

    ros::spinOnce();

    loop_rate.sleep();
  }

  return 0;
}
