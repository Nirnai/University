#include <ros/ros.h>
#include "object_detection.h"

int main(int argc, char	**argv)	{
  // Initiate new ROS node named "ObjectDetection"
  ros::init(argc, argv,	"object_detection");

  // Create new stopper object
  ObjectDetection ObjectDetection;

  // Start the movement
  ObjectDetection.startObjectDetection();

  return 0;
}
