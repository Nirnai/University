//#include <math.h>
#include <algorithm>
#include <geometry_msgs/Twist.h>
#include <nav_msgs/Odometry.h>
#include <std_msgs/Bool.h>
#include "turtlebot.h"


#define PI 3.141592


Turtlebot::Turtlebot()
{
    listner = new tf::TransformListener;
    move = false;
    reverse = false;
    turn_left = false;
    turn_right = false;
    obsticle = false;
    reachedlocal = true;
    reachedglobal.data = false;
    xmin = -100;
    xmax = 100;
    ymin = -100;
    ymax = 100;

    FORWARD_SPEED_MPS = 0;
    ANGULAR_SPEED_MPS = 0;

    // Advertise a new publisher for the simulated robot's velocity command topic
    commandPub = node.advertise<geometry_msgs::Twist>("/cmd_vel_mux/input/navi", 10);
    reachedGoalPub = node.advertise<std_msgs::Bool>("ReachedGoal", 1000);
    odomSub = node.subscribe("/odom", 1000, &Turtlebot::odomCallback, this);
    laserSub = node.subscribe("/lidarscan", 1, &Turtlebot::scanCallback, this);
    goalSub = node.subscribe("/GoalPosition", 1, &Turtlebot::goalCallback, this);
    stopSub = node.subscribe("/communication/stopImmediately",1, &Turtlebot::stopCallback,this);
    aSub = node.subscribe("/communication/correctA", 1, &Turtlebot::aCallback,this);
    bSub = node.subscribe("/communication/correctB", 1, &Turtlebot::bCallback,this);
    shotGoalSub = node.subscribe("/statemachine/shotGoal", 1, &Turtlebot::shotGoalCallback,this);

}


void Turtlebot::shotGoalCallback(const std_msgs::Bool::ConstPtr& msg){
  if(msg->data){
    reverse = true;
  }else{
    reverse = false;
  }
}

void Turtlebot::aCallback(const std_msgs::Float32::ConstPtr& A){
  this->a = A->data;
  this->xmax = 100;// 5/2 * a - 0.05;
  this->xmin = -100; //- (1/2 * a) + 0.05;
  this->ymax = 100; //b/2 - 0.05;
  this->ymin = -100; //-b/2 + 0.05;
}
void Turtlebot::bCallback(const std_msgs::Float32::ConstPtr& B){
  this->b = B->data;
}

void Turtlebot::stopCallback(const std_msgs::Bool::ConstPtr& stopImmediately){
  if(stopImmediately->data){
    this->stopImmediately = true;
  }else{
    this->stopImmediately = false;
  }
}


void Turtlebot::odomCallback(const nav_msgs::Odometry::ConstPtr& odom){

  // tf::StampedTransform transform1;
  // tf::StampedTransform transform2;
  //
  // try{
  //   this->listner->lookupTransform("/map","/base_footprint",
  //                            ros::Time(0), transform1);/GoalPosition
  // }
  // catch (tf::TransformException &ex) {
  //   ROS_ERROR("%s",ex.what());
  //   ros::Duration(1.0).sleep();
  // }
  // try{
  //   this->listner->lookupTransform("/odom","/map",
  //                            ros::Time(0), transform2);
  // }
  // catch (tf::TransformException &ex) {
  //   ROS_ERROR("%s",ex.what());
  //   ros::Duration(1.0).sleep();
  // }
  this->x = odom->pose.pose.position.x;
  this->y = odom->pose.pose.position.y;
  double quatx= odom->pose.pose.orientation.x;
  double quaty= odom->pose.pose.orientation.y;
  double quatz= odom->pose.pose.orientation.z;
  double quatw= odom->pose.pose.orientation.w;

  tf::Quaternion q(quatx, quaty, quatz, quatw);
  tf::Matrix3x3 m(q);
  double roll, pitch, yaw;
  m.getRPY(roll, pitch, yaw);

  this->theta = -yaw;

  this->x_vel = odom->twist.twist.linear.x;
  this->y_vel = odom->twist.twist.linear.y;
  this->theta_vel = odom->twist.twist.angular.z;
  //ROS_INFO("OdomPos = [%f],[%f],[%f]", this->y,this->x,this->theta);

  // this->x = transform1.getOrigin().x() - transform2.getOrigin().x();
  // this->y = transform1.getOrigin().y() - transform2.getOrigin().y();
  // this->theta = -tf::getYaw(transform1.getRotation()) - tf::getYaw(transform2.getRotation());
  //ROS_INFO("TfPos = [%f],[%f],[%f]", this->y,this->x,this->theta);
}

void Turtlebot::goalCallback(const geometry_msgs::Pose::ConstPtr& goal){

  // this->x_globalgoal = goal->position.x;
  // this->y_globalgoal = goal->position.y;

  if(goal->position.x > xmin && goal->position.x < xmax){
    this->x_globalgoal = goal->position.x;
  }else{
    this->x_globalgoal = goal->position.x + 0.01;
    ROS_ERROR("X Goal is out of bounds!!!!");
  }
  if(goal->position.y > ymin && goal->position.y < ymax){
    this->y_globalgoal = goal->position.y;
  }else{
    this->y_globalgoal = goal->position.y;
    ROS_ERROR("Y Goal is out of bounds!!!!");
  }

  ROS_INFO("Goal = [%f],[%f]", this->x_globalgoal, this->y_globalgoal);

}

// Process the incoming laser scan message
void Turtlebot::scanCallback(const sensor_msgs::LaserScan::ConstPtr& scan){
  //ROS_INFO("Scan Callback!");

  int minIndex1 = 0;
  int maxIndex1 = 30;
  int minIndex2 = 210;
  int maxIndex2 = 359;

  //int minIndex = ceil((MIN_SCAN_ANGLE_RAD - scan->angle_min) / scan->angle_increment);
  //int maxIndex = floor((MAX_SCAN_ANGLE_RAD - scan->angle_min) / scan->angle_increment);
  float closestRange = scan->ranges[minIndex1];
  for (int currIndex = minIndex1; currIndex <= maxIndex1; currIndex++) {
    if (scan->ranges[currIndex] < closestRange) {
      closestRange = scan->ranges[currIndex];
    }
  }
  for (int currIndex = minIndex2; currIndex <= maxIndex2; currIndex++) {
    if (scan->ranges[currIndex] < closestRange) {
      closestRange = scan->ranges[currIndex];
    }
  }

  //ROS_INFO("closestRange = [%f]", closestRange);


  if (closestRange < MIN_PROXIMITY_RANGE_M) {

    ROS_INFO("Obsticle!");
    move = false;
    obsticle = true;
    reachedlocal = false;
    x_localgoal = this->x+0.4;
    y_localgoal = this->y+0.4;
  }else{
    //ROS_INFO("NO Obsticle!");
    obsticle = false;
  }
}


double Turtlebot::getTheta(const nav_msgs::Odometry::ConstPtr& odom){
  // yaw (z-axis rotation)
  double siny = +2.0 * (odom->pose.pose.orientation.w * odom->pose.pose.orientation.z + odom->pose.pose.orientation.x * odom->pose.pose.orientation.y);
  double cosy = +1.0 - 2.0 * (odom->pose.pose.orientation.y * odom->pose.pose.orientation.y + odom->pose.pose.orientation.z * odom->pose.pose.orientation.z);
  return atan2(siny, cosy) + PI;
}

void Turtlebot::decideTurnDirection(){
  if(this->theta > 0 && this->theta_goal > 0 || this->theta < 0 && this->theta_goal < 0){
    if(this->theta>this->theta_goal){
      turn_left = true;
      turn_right = false;
    }else{
      turn_left = false;
      turn_right = true;
    }
  }
  else if(this->theta > 0 && this->theta_goal < 0){
    if((PI - this->theta + PI + this->theta_goal) < PI){
      turn_left = false;
      turn_right = true;
    }else{
      turn_left = true;
      turn_right = false;
    }
  }
  else if(this->theta < 0 && this->theta_goal > 0){
    if(PI + this->theta + PI-this->theta_goal < PI){
      turn_left = true;
      turn_right = false;
    }else{
      turn_left = false;
      turn_right = true;
    }
  }
  else{
    turn_left = true;
    turn_right = false;
  }
}

void Turtlebot::checkGoal(double x_goal, double y_goal, double theta_goal){
  // Check orientation goal
  double error_theta = fabsf(this->theta_goal-this->theta);
  double error_theta_vel = fabsf(0-this->theta_vel);
  double error_xy = sqrt(pow(this->x_goal - this->x,2) + pow(this->y_goal- (-this->y),2));


  ROS_INFO("Error_theta = [%f]", error_theta);
  ROS_INFO("Error_xy = [%f]", error_xy);
  //ROS_INFO("theta = [%f]", theta);
  //ROS_INFO("theta_goal = [%f]", theta_goal);

  // Check Orientational error
  if(error_theta < 0.03)
  {
    turn_left = false;
    turn_right = false;
    // Check Positional Error
    if(error_xy < 0.03)
    {
      move = false;
      this->reachedglobal.data=true;
      ROS_INFO("REACHED GOAL!!!!!!");
    }else{
      move = true;
      this->reachedglobal.data=false;
      if(error_xy < 0.2)
      {
        this-> FORWARD_SPEED_MPS = (error_xy * K_FORWARD);
              //ROS_INFO("Jumped in P-Control");
              if(error_xy<0.03 && this->avoid){
                this->reachedlocal = true;
              }else if(error_xy<0.03 && !this->avoid){

              }
      }else{
        this-> FORWARD_SPEED_MPS = 0.1;
      }
    }
  // Turn until orientaion goal is reached.
  }else{
    this->reachedglobal.data=false;
    // Both are same sign
    decideTurnDirection();
    // Set turn speed
    move = false;
    //this-> ANGULAR_SPEED_MPS = 0.3;
    //this-> ANGULAR_SPEED_MPS = Kp * error_theta - Kd * error_theta_vel;

    if(error_theta < 0.1)
    {
      this-> ANGULAR_SPEED_MPS = (error_theta * K_ANGULAR);
      ROS_INFO("Jumped in P-Control");

    }else{
      this-> ANGULAR_SPEED_MPS = 0.4;
    }

    //ROS_INFO("ANGULAR_SPEED_MPS = [%f]", this-> ANGULAR_SPEED_MPS);


  }


}

void Turtlebot::moveForward() {
    geometry_msgs::Twist msg;           // The default constructor will set all commands to 0
    msg.linear.x = FORWARD_SPEED_MPS;
    commandPub.publish(msg);
}

void Turtlebot::turnRight() {
  geometry_msgs::Twist msg;
  msg.angular.z = -ANGULAR_SPEED_MPS;
  commandPub.publish(msg);

}

void Turtlebot::turnLeft() {
    geometry_msgs::Twist msg;
    msg.angular.z = ANGULAR_SPEED_MPS;
    commandPub.publish(msg);
}

void Turtlebot::moveBackward() {
  geometry_msgs::Twist msg;
  msg.linear.x = - FORWARD_SPEED_MPS;
  commandPub.publish(msg);
  ros::Duration(3).sleep();
  msg.linear.x = 0;
  commandPub.publish(msg);
}

void Turtlebot::stop(){
    geometry_msgs::Twist msg;
    msg.angular.z = 0;
    msg.linear.x = 0;
    commandPub.publish(msg);
}


void Turtlebot::TurtlebotStart() {
  ros::Rate rate(10);
  ROS_INFO("Start");
  // Keep spinning loop until user presses Ctrl+C
  //this->y_globalgoal = this->y;
  //this->x_globalgoal = this->x+0.0001;
  //this->theta_goal = PI/2;

  static tf::TransformBroadcaster br;


  while (ros::ok()) {
      ROS_INFO("x_goal, y_goal, x, y,theta_goal, theta = [%f],[%f],[%f],[%f],[%f],[%f]", this->x_goal,this->y_goal, this->x,this->y,this->theta_goal,this->theta);
      tf::Transform transform;
      transform.setOrigin( tf::Vector3(this->x, this->y, 0.0) );
      tf::Quaternion q;
      q.setRPY(0, 0, this->theta);
      transform.setRotation(q);
      br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "/map", "/turtle"));
      //if(obsticle){
        //this->avoid = true;
        //this->x_goal = this->x_localgoal;
        //this->y_goal = this->y_localgoal;
      //}else if(!obsticle && reachedlocal){
        this->x_goal = this->x_globalgoal;
        this->y_goal = this->y_globalgoal;
        this->avoid = false;
      //}

        // double ydiff = roundf(this->y_goal-(-this->y) * 1000) / 1000;
        // double xdiff = roundf(this->x_goal - this->x * 1000) / 1000;
        double ydiff = this->y_goal-(-this->y);
        double xdiff = this->x_goal - this->x;
        // if(ydiff == -0){
        //   ydiff =0;
        // }
        // if(xdiff == -0){
        //   xdiff =0;
        // }

        ROS_INFO("ydiff = [%f], xdiff=[%f]", ydiff, xdiff);
        this->theta_goal = atan2(ydiff,xdiff);
      // Continue moving if scan changes (we or the obstacle gets moved)

      checkGoal(x_goal, y_goal, theta_goal);
      if (move && !stopImmediately/*&& !obsticle*/){
        moveForward();
        ROS_INFO("Moving!");
      }else if(turn_right){
        turnRight();
        ROS_INFO("Turning right!");
      }else if(turn_left){
        turnLeft();
        ROS_INFO("Turning left!");
      }else if(reverse){
        moveBackward();
        ROS_INFO("reversing");
      }else{
        stop();
        ROS_INFO("Stopped!");

      }
      //ROS_INFO("Reached [%d]", this->reachedglobal.data);
      this->reachedGoalPub.publish(this->reachedglobal);

      ros::spinOnce(); // Need to call this function often to allow ROS to process incoming messages
      rate.sleep();
    }
}
