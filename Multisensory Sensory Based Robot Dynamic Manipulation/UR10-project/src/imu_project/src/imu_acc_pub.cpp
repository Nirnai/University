/**
 *
 * This node translates the imu acceleration from the enu/world frame to the kinect_top
 * where it is processed.
 */

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "sensor_msgs/Imu.h"
#include "geometry_msgs/QuaternionStamped.h"
#include "tf/transform_listener.h"
#include "tf/transform_broadcaster.h"
#include "tf_conversions/tf_eigen.h"




geometry_msgs::Vector3Stamped acc_enu;
geometry_msgs::Vector3Stamped vel_angular_imu;
geometry_msgs::PointStamped kinect_point_raw;

void kinectGrabbing(const geometry_msgs::Point::ConstPtr &msg)
{
  if (!(isnan(msg->x) || isnan(msg->y) || isnan(msg->z))) {
    kinect_point_raw.point.x = msg->x;
    kinect_point_raw.point.y = msg->y;
    kinect_point_raw.point.z = msg->z;
  }
}

void acc_grabbing(const sensor_msgs::ImuConstPtr& msg)
{

  acc_enu.vector = msg->linear_acceleration;
  vel_angular_imu.vector = msg->angular_velocity;
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "imu_acc_pub");
  ros::NodeHandle nh;

  ros::Rate loop_rate(125);

  ros::Subscriber sub = nh.subscribe("imu/data",1000,&acc_grabbing);
  ros::Subscriber sub_kinect = nh.subscribe("/objectTracker/position",1000,&kinectGrabbing);
  ros::Publisher pub = nh.advertise<geometry_msgs::Vector3Stamped>("/acc_base",100);
  ros::Publisher ang_vel_pub = nh.advertise<geometry_msgs::Vector3Stamped>("/ang_vel_base",100);
  ros::Publisher pub_kinect = nh.advertise<geometry_msgs::PointStamped>("/kinect_transformed",100);

  acc_enu.header.frame_id = "/world";
  vel_angular_imu.header.frame_id = "/imu";

  geometry_msgs::Vector3Stamped acc_base;
  geometry_msgs::Vector3Stamped vel_angular_base;


  geometry_msgs::PointStamped kinect_robot_base;
  kinect_robot_base.header.frame_id = "/robot_base";
  kinect_point_raw.header.frame_id = "/kinect";

  acc_base.header.frame_id = "/robot_base";
  vel_angular_base.header.frame_id = "/robot_base";

  tf::TransformListener listener;

  listener.waitForTransform("/robot_base","/imu",ros::Time::now(), ros::Duration(0.1));
  listener.waitForTransform("/robot_base","/kinect",ros::Time::now(), ros::Duration(0.1));

  while (ros::ok())
  {

    try{
      ros::Time now = ros::Time::now();
      listener.transformVector("/robot_base",acc_enu,acc_base);
      listener.transformPoint("/robot_base",kinect_point_raw,kinect_robot_base);
      listener.transformVector("/robot_base",vel_angular_imu,vel_angular_base);
    }
    catch (tf::TransformException ex){
      ROS_INFO_STREAM("Error");
    }

    pub.publish(acc_base);
    pub_kinect.publish(kinect_robot_base);
    ang_vel_pub.publish(vel_angular_base);

    ros::spinOnce();
    loop_rate.sleep();
  }

  return 0;
}
