#include <ros/ros.h>
#include "object_tracker.h"

int main(int argc, char	**argv)	{

    ros::init(argc, argv,"object_tracker");

    ObjectTracker ObjectTracker;

    ObjectTracker.run();


    return 0;
}