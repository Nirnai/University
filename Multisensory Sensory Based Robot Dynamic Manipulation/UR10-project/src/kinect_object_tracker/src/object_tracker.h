#ifndef OBJECT_TRACKING_H
#define OBJECT_TRACKING_H


#include <ros/ros.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/Point.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PCLPointCloud2.h>

#include <cv_bridge/cv_bridge.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <tf/transform_broadcaster.h>

#include <vector>
#include <Eigen/Core>
#include <numeric>
#include <algorithm>
#include <stdint.h>



class ObjectTracker {
    public:
        ObjectTracker();

        void run();
        void findObjectPosition(int meas_count);  
        void findObjectContour();

    private:
        // Nodehandle
        ros::NodeHandle node;

        // Subscriber
        ros::Subscriber imgSub;
        ros::Subscriber pointSub;

        // Publisher
        ros::Publisher hsvPub;
        ros::Publisher colorPub;
        ros::Publisher circlePub;
        ros::Publisher positionPub;

        // Member Variables
        sensor_msgs::Image::ConstPtr image;

        // CV Brdge Variables
        cv_bridge::CvImageConstPtr cv_ptr;
        cv_bridge::CvImage hsv_msg;
        cv_bridge::CvImage color_msg;

        // Opencv Variables
        cv::Mat hsvImg;
        cv::Mat redImgUpper;
        cv::Mat redImgLower;
        cv::Mat redHueImg;
        cv::Mat redColorImg;
        cv::Mat circImg;

        // Variables needed to find Object
        // std::vector<std::vector<cv::Point>> contours_poly;
        // std::vector<cv::Rect> boundRect;
        // std::vector<cv::Point2f> center;
        // std::vector<float> radius;

        cv::Rect largestRect;
        cv::Point centerRect;
        //cv_bridge::CvImageConstPtr hsvImg;

        // PCL Variables
        // pcl::PCLPointCloud2* pcl_cloud = new pcl::PCLPointCloud2;;
        // pcl::PCLPointCloud2ConstPtr pcl_cloudPtr;
        // pcl::PointCloud<pcl::PointXYZ> XYZcloud;

        // sensor_msgs/Pointcloud2 data
        sensor_msgs::PointCloud2 pcloud;

        // Published Msgs
        
        geometry_msgs::Point position;

        // tf variables
        tf::TransformBroadcaster *br = new tf::TransformBroadcaster();
        tf::Transform transform;

        int measurmentFrame;
        std::vector<float> xvec;
        std::vector<float> yvec;
        std::vector<float> zvec;
        Eigen::Vector3d meas;
        Eigen::Vector3d filtered_pos;

        // Subscriber Callbacks
        void imageCallback(const sensor_msgs::Image::ConstPtr  &img);
        void pointsCallback(const sensor_msgs::PointCloud2::ConstPtr &cloud);
        void EMA(Eigen::Vector3d &filtered_pos, Eigen::Vector3d &meas, float alpha, int count);

};

#endif // OBJECT_TRACKING_H