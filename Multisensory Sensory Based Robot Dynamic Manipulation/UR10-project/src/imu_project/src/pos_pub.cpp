#include "ros/ros.h"
#include "std_msgs/String.h"
#include <math.h>

#include "geometry_msgs/PointStamped.h"
#include "tf/transform_broadcaster.h"

#define PI 3.14159265

int main(int argc, char **argv)
{
  ros::init(argc, argv, "pos_pub");
  ros::NodeHandle nh;

  ros::Publisher pos_pub = nh.advertise<geometry_msgs::Point>("/objectTracker/position", 1000);

  ros::Rate loop_rate(30);

  geometry_msgs::Point point_kinect;

  point_kinect.x = 0;
  point_kinect.y = 0;
  point_kinect.z = 1;

  tf::TransformBroadcaster br;
  tf::Transform transform;
  tf::Quaternion q;

  q.setRPY(0,0,0);
  q.normalize();
  transform.setRotation(q);


  while (ros::ok())
  {

    transform.setOrigin(tf::Vector3(point_kinect.x,point_kinect.y,point_kinect.z));
    br.sendTransform(tf::StampedTransform(transform,ros::Time::now(),"/kinect","/kinect_simulation"));

    pos_pub.publish(point_kinect);

    ros::spinOnce();

    loop_rate.sleep();
  }

  return 0;
}
