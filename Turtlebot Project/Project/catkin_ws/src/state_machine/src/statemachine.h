#ifndef STATEMACHINE_H
#define STATEMACHINE_H

#include <ros/ros.h>
#include <geometry_msgs/Pose.h>
#include <std_msgs/Byte.h>
#include <std_msgs/Bool.h>
#include <geometry_msgs/TransformStamped.h>
#include <tf/transform_listener.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Pose2D.h>

class Statemachine{
    public:
      Statemachine();

      int state;
      bool start;
      bool startGame;
      bool abDone;
      bool teamfarbeDone;
      bool reachedGoal;
      bool puckFound;



      std_msgs::Byte OD_mode;
      std_msgs::Byte IP_mode;
      geometry_msgs::Pose goal;
      std_msgs::Bool isAlive;
      std_msgs::Bool shotGoal;

      ros::Subscriber startSub;
      ros::Subscriber startGameSub;
      ros::Subscriber abDoneSub;
      ros::Subscriber teamfarbeDoneSub;
      ros::Subscriber reachedSub;
      ros::Subscriber odomSub;

      ros::Subscriber puckAtMidSub;
      ros::Subscriber puckGoal;

      ros::Publisher goalPub;
      ros::Publisher ODPub;
      ros::Publisher IPPub;
      ros::Publisher isAlivePub;
      ros::Publisher shotGoalPub;

      double x_vel;
      double theta_vel;



    private:
      double x;
      double y;
      double theta;



    ros::NodeHandle n;


    void startgameCallback(const std_msgs::Bool::ConstPtr& startGame);
    void goalCallback(const std_msgs::Bool::ConstPtr& reached);
    void startCallback(const std_msgs::Bool::ConstPtr& go);
    void abDoneCallback(const std_msgs::Bool::ConstPtr& done);
    void teamfarbeCallback(const std_msgs::Bool::ConstPtr& done);
    void odomCallback(const nav_msgs::Odometry::ConstPtr& odom);
    void puckAtMidCallback(const std_msgs::Bool::ConstPtr& msg);
    void puckGoalCallback(const geometry_msgs::Pose2D::ConstPtr& msg);

};

#endif
