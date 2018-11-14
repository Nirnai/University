#include "ros/ros.h"
#include "std_msgs/String.h"
#include "geometry_msgs/Vector3Stamped.h"
#include <queue>
#include "std_msgs/Bool.h"
#include <math.h>

std::queue<double> measurement_queue;
int counter_high = 0;
int counter_low = 0;

void shakyCallback(const geometry_msgs::Vector3Stamped &msg)
{
  double accel_treshold = 3;
  double sum_accel = msg.vector.x;

  measurement_queue.push(sum_accel);

  if (sum_accel > accel_treshold)
  {
    counter_high++;
  }
  else if (sum_accel < -accel_treshold)
  {
    counter_low++;
  }

  if (measurement_queue.size() > 150) {
    measurement_queue.pop();
    double old_accel = measurement_queue.front();
    // reduce counters if element is popped
    if (old_accel > accel_treshold)
    {
      counter_high--;
    }
    else if (old_accel < -accel_treshold)
    {
      counter_low--;
    }
  }
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "shake_sensor");
  ros::NodeHandle nh;

  ros::Rate loop_rate(100);
  ros::Subscriber sub_accel = nh.subscribe("/kalman/filtered_accel", 1000, shakyCallback);
  ros::Publisher shake_pub = nh.advertise<std_msgs::Bool>("/shaking_detected", 1000);


  std_msgs::Bool shake_bool;
  shake_bool.data = false;

  while (ros::ok())
  {
    ROS_INFO_STREAM(counter_high);
    ROS_INFO_STREAM(counter_low);

    if ((counter_low > 30) && (counter_high > 30) && abs(counter_high - counter_low) < 5)
    {
      ros::Duration(2).sleep();
      shake_bool.data = true;
      shake_pub.publish(shake_bool);
      shake_bool.data = false;
      ros::Duration(1).sleep();
      measurement_queue = std::queue<double>();
      counter_high = 0;
      counter_low = 0;
    }

    ros::spinOnce();
    loop_rate.sleep();

  }

  return 0;
}
