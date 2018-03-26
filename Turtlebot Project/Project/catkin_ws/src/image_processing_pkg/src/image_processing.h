#ifndef IMAGE_PROCESSING_H
#define IMAGE_PROCESSING_H

#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/CompressedImage.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <std_msgs/String.h>
#include <std_msgs/Byte.h>
#include <std_msgs/Bool.h>

#include <cv_bridge/cv_bridge.h>
#include <cv_bridge/rgb_colors.h>
#include <image_transport/image_transport.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d.hpp>

#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/Pose2D.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseArray.h>

#include <image_processing_pkg/obstacle_poses.h>

#include <vector>
#include <algorithm>


typedef struct {
    int iLowH;
    int iHighH;
    int iLowS;
    int iHighS;
    int iLowV;
    int iHighV;
} HSVColorRange;



using namespace cv;

class ImageProcessing {
  public:
    // Tunable parameters
    const static int yellowFieldDetectionThres = 1000;

    ImageProcessing();
    void startBV();

    bool calculateTeamColor();
    void saveIncomingImage();
    int detectGreenLimitPost();

  private:
    ros::NodeHandle node;
    ros::Subscriber img;    // Subscriber to the robots image topic
    ros::Subscriber sub_points;
    ros::Publisher  teamColorPub;
    ros::Publisher  drivingGoalPub;
    ros::Publisher  alivePub;
    ros::Publisher initTeamColorSuccessPub;
    ros::Publisher puckAtMidPub;
    ros::Subscriber modiSub;
    ros::Subscriber correctTeamColorSub;
    ros::Subscriber obstaclesSub;

    std_msgs::String       teamColor;
    std_msgs::String       correctTeamColor;
    std_msgs::ByteConstPtr currentCalcMode;
    std_msgs::Bool         isAlive;
    std_msgs::Bool         initTCOkay;
    std_msgs::Bool         puckAtMidFound;

    geometry_msgs::Pose2D bestPuck;

    image_processing_pkg::obstacle_poses obstaclesArround;

    sensor_msgs::Image::ConstPtr     image;
    cv_bridge::CvImageConstPtr       cv_ptr;
    sensor_msgs::PointCloud2ConstPtr global_cloud_msg;

    std::vector<float> puckFeatures;

    bool yellowTeam;
    bool blueTeam;
    bool receivecFirstMode;
    bool testWithPictureMode;   // If 1, no connection to turtlebot is needed, instead previous taken images are loaded
    bool correctTeamColorDefined;
    const char* testPictureFilename;
    bool flag_image_ok_point;
    bool receivedFirstCloud;
    bool receivecFirstObstacles;


    bool shapeDetection(std::vector<std::vector<Point> > &contours, std::vector<std::vector<Point> > &approxPolygons,
                        std::vector<std::vector<double> > &huDescriptor, std::vector<Point2f> &mc, HSVColorRange colour,
                        int cutImage, bool showResults);

    void shapeMatching(std::vector<std::vector<Point> > &contours, std::vector<std::vector<Point> > &approxPolygons,
                       std::vector<std::vector<double> > &huDescriptor, const char *currentColour);

    void drainClassifier();

    geometry_msgs::Point pixelTo3DPoint(int row_num,  int col_num);
    void update();
    void scanImageCallback( const sensor_msgs::Image::ConstPtr  & img);
    void cloudCallback(const sensor_msgs::PointCloud2ConstPtr& cloud_msg);
    void modeCallback(const std_msgs::ByteConstPtr &mode);
    void correctTeamColorCallback(const std_msgs::String &str);
    void obsCallback(const image_processing_pkg::obstacle_poses &obs);
    bool findNearestObstacle(geometry_msgs::Pose2D &bestMatch, float x, float y, float theta);

};

#endif // IMAGE_PROCESSING_H
