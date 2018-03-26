#include <ros/ros.h>
#include "image_processing.h"

int main(int argc, char	**argv)	{
    ros::init(argc, argv,	"image_processing");
    ImageProcessing ImageProcessing;

    ImageProcessing.startBV();

    return 0;
}
