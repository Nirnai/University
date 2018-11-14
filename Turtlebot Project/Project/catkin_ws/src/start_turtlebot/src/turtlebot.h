#ifndef TURTLEBOT_H
#define TURTLEBOT_H

#include <ros/ros.h>
#include <geometry_msgs/TransformStamped.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float32.h>

class Turtlebot{
    public:
        // Tunable parameters
        double FORWARD_SPEED_MPS;
        double ANGULAR_SPEED_MPS;
        // K value for p-control
        const static double K_FORWARD = 0.8;
        const static double K_ANGULAR = 8;
        const static double Kp = 0.1;
        const static double Kd = 0.1;


        const static double MIN_SCAN_ANGLE_RAD = -30.0/180*M_PI;
        const static double MAX_SCAN_ANGLE_RAD = +30.0/180*M_PI;
        const static float MIN_PROXIMITY_RANGE_M = 0.2; // Should be smaller than sensor_msgs::LaserScan::range_max

        Turtlebot();
        void TurtlebotStart();
        // getter
        double getTheta(const nav_msgs::Odometry::ConstPtr& odom);

        // setter
        void setOrientation(double theta);

    private:
        ros::NodeHandle node;
        ros::Publisher commandPub;
        ros::Publisher reachedGoalPub;
        ros::Subscriber laserSub;
        ros::Subscriber odomSub;
        ros::Subscriber goalSub;
        ros::Subscriber stopSub;
        ros::Subscriber aSub;
        ros::Subscriber bSub;
        ros::Subscriber shotGoalSub;
        tf::TransformListener* listner;


        double x;
        double y;
        double theta;

        double x_vel;
        double y_vel;
        double theta_vel;

        double a;
        double b;
        double xmax;
        double xmin;
        double ymax;
        double ymin;

        double x_goal;
        double y_goal;
        double theta_goal;
        double x_globalgoal;
        double y_globalgoal;
        double theta_globalgoal;
        double x_localgoal;
        double y_localgoal;
        double theta_localgoal;

        bool avoid;
        bool reachedlocal;
        std_msgs::Bool reachedglobal;
        bool move;
        bool reverse;
        bool stopImmediately;
        bool turn_right;
        bool turn_left;
        bool obsticle;

        void moveForward();
        void turnRight();
        void turnLeft();
        void moveBackward();
        void stop();
        void checkGoal(double x_goal, double y_goal, double theta_goal);
        void scanCallback(const sensor_msgs::LaserScan::ConstPtr& scan);
        void odomCallback(const nav_msgs::Odometry::ConstPtr& odom);
        void stopCallback(const std_msgs::Bool::ConstPtr& stopImmediately);
        void aCallback(const std_msgs::Float32::ConstPtr& A);
        void bCallback(const std_msgs::Float32::ConstPtr& B);
        void shotGoalCallback(const std_msgs::Bool::ConstPtr& msg);

        void goalCallback(const geometry_msgs::Pose::ConstPtr& goal);
        void decideTurnDirection();
};

#endif
