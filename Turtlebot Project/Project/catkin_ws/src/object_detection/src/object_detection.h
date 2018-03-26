#ifndef OBJECT_DETECTION_H
#define OBJECT_DETECTION_H

#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/CompressedImage.h>
#include <sensor_msgs/Image.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Byte.h>
#include <std_msgs/Float32.h>
#include <nav_msgs/OccupancyGrid.h>
#include <nav_msgs/MapMetaData.h>
#include <vector>
#include <algorithm>

#include <object_detection/goal_edge_points.h>
#include <object_detection/obstacle_poses.h>
#include <object_detection/game_field.h>

typedef struct {
    int index_min;
    int index_max;
    float angle_min;
    float angle_max;
    float angle_mean;
    float distance;
    float x;
    float y;
//    float distancesToOtherObj[];
    std::vector<float> distancesToOtherObj;
} ObstacleInLidar;

struct Quaternion {
  double x;
  double y;
  double z;
  double w;
};

class ObjectDetection {
  public:
    // Tunable parameters

    const static int yellowFieldDetectionThres = 5000;
//    std::vector<float> rightStackA;
//    std::vector<float> leftStackA;
//    std::vector<float> stackB;
//    int measurements;

    std::vector<ObstacleInLidar> obstacleStack;
    std::vector<ObstacleInLidar> obstacleStackCompl;

    ObjectDetection();
    void startObjectDetection();

  private:
    ros::NodeHandle node;
    ros::Publisher alivePub; // Publisher to the robot's velocity command topic
    ros::Publisher abRatioPub;
    ros::Publisher ownPositionPub;
    ros::Publisher ownPosition2DPub;
    ros::Publisher enemyGoalPub;
    ros::Publisher gameFieldPub;
    ros::Publisher obstaclePosesPub;
    ros::Publisher initABSuccessPub;

    ros::Subscriber laserSub; // Subscriber to the robot's laser scan topic
    ros::Subscriber gmapSub;
    ros::Subscriber modiSub;

    sensor_msgs::LaserScan::ConstPtr lidarScan;

    std_msgs::Bool isAlive;
    std_msgs::Bool initABOkay;
    std_msgs::ByteConstPtr currentCalcMode; // different modi: 1 = try Measurement from start Position, 2: measure in front,
    std_msgs::Float32 abRatio;

    object_detection::goal_edge_points enemyGoal;
    object_detection::obstacle_poses obsPoses;
    object_detection::game_field gameFieldCorners;

    bool receivecFirstLidar;
    bool receivecFirstMode;

    bool foundSomethingOnLeftSide;
    bool foundSomethingOnRightSide;

    // Variables for a,b calculation
    std::vector<ObstacleInLidar> rightObstacleStack;
    std::vector<ObstacleInLidar> leftObstacleStack;
    std::vector<float> rightStackA;
    std::vector<float> leftStackA;
    std::vector<float> stackB;

    std::vector<float> posX;
    std::vector<float> posY;
    std::vector<float> quatZ;
    std::vector<float> quatW;
    std::vector<float> posTheta;

    geometry_msgs::Pose ownPosition;
    geometry_msgs::Pose2D ownPosition2D;

    int outputMode; // 0 = none, 1 = some, 2 = all
    float rangeMultMin; // Scan radius multiplicator from min
    float objectScanOffset; // Distance difference threshold where lidar points are counted to one object
    int indexOffsetToObjRow1; // index offset to adress distances between obj correct (workaround)
    int indexOffsetToObjRow2;
    int measurements;
    bool measurementsDone;

    float meanRA;
    float varRA;
    float medianRA;
    float meanLA;
    float varLA;
    float medianLA;
    float meanB;
    float varB;
    float medianB;

    int leftIndexForB;
    int rightIndexForB;
    double begin;
    double runtime;
    double diff;

    void scanCallback(const sensor_msgs::LaserScan::ConstPtr& scan);
    void mapCallback(const nav_msgs::OccupancyGridConstPtr &map);
    void modeCallback(const std_msgs::ByteConstPtr &mode);

    geometry_msgs::Quaternion toQuaternion(float pitch, float roll, float yaw);
    void scanForObjects(float scanInDist);
    bool initLocalisation(std::vector<ObstacleInLidar> &postObstacleStack, int &indexOffsetToObjRowSave, float minAngle, float maxAngle,
                          int minAmountOfPoints, float rangeMultMin, float objectScanOffset, int outputMode);
    bool plausibilityCheckForThreeObj(int fstIndex, int sndIndex, int thrIndex, int outputMode, int offset);
    void searchMinDistanceInLidarRange(float &minRange, float minAngle, float maxAngle);
    void calculateMeanVarMedian(float &mean, float &var, float &median, std::vector<float> values, bool sorted);

    void calculateAandBFromStartPos();
    void update();
};

#endif // OBJECT_DETECTION_H
