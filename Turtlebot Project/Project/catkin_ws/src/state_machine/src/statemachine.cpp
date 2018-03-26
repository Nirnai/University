#include <ros/ros.h>
#include "statemachine.h"



Statemachine::Statemachine(){

  state = 4;
  start = true;
  startGame = false;
  reachedGoal = false;
  abDone = false;
  teamfarbeDone = false;
  goal.position.x = this->x+0.01;
  goal.position.y = this->y;
  OD_mode.data = 0;
  IP_mode.data = 0;
  isAlive.data = false;


  startSub = n.subscribe("/communication/startInit", 1, &Statemachine::startCallback, this);
  startGameSub = n.subscribe("/communication/startGame", 1, &Statemachine::startgameCallback,this);
  abDoneSub = n.subscribe("/object_detection/initABSuccess", 1, &Statemachine::abDoneCallback, this);
  teamfarbeDoneSub = n.subscribe("/image_processing/initTCSuccess", 1, &Statemachine::teamfarbeCallback, this);

  puckAtMidSub = n.subscribe("/image_processing/puckAtMid", 1, &Statemachine::puckAtMidCallback, this);
  puckGoal = n.subscribe("/image_processing/puckGoal", 1, &Statemachine::puckGoalCallback, this);

  reachedSub = n.subscribe("ReachedGoal", 1, &Statemachine::goalCallback, this);
  odomSub = n.subscribe("/odom", 1000, &Statemachine::odomCallback, this);
  goalPub = n.advertise<geometry_msgs::Pose>("/GoalPosition", 10);
  ODPub = n.advertise<std_msgs::Byte>("/object_detection/mode", 10);
  IPPub = n.advertise<std_msgs::Byte>("/image_processing/mode", 10);
  isAlivePub = n.advertise<std_msgs::Bool>("/state_machine/isAlive", 10);
  shotGoalPub = n.advertise<std_msgs::Bool>("/statemachine/shotGoal", 10);
}

void Statemachine::puckAtMidCallback(const std_msgs::Bool::ConstPtr& msg){
  if(msg->data){
    this->puckFound = true;
  }else{
    this->puckFound = false;
  }
}

void Statemachine::puckGoalCallback(const geometry_msgs::Pose2D::ConstPtr& msg){
  if(this->puckFound && this->state == 6){
    this->goal.position.x = msg->x;
    this->goal.position.y = msg->y;
  }
}

void Statemachine::startgameCallback(const std_msgs::Bool::ConstPtr& startGame){
  if(startGame->data){
    this->startGame = startGame->data;
  }else{
    this->startGame = false;
  }
}


void Statemachine::odomCallback(const nav_msgs::Odometry::ConstPtr& odom){
  this->x = odom->pose.pose.position.x;
  this->y = odom->pose.pose.position.y;
  this->x_vel = odom->twist.twist.linear.x;
  this->theta_vel = odom->twist.twist.angular.z;
  double quatx= odom->pose.pose.orientation.x;
  double quaty= odom->pose.pose.orientation.y;
  double quatz= odom->pose.pose.orientation.z;
  double quatw= odom->pose.pose.orientation.w;
  tf::Quaternion q(quatx, quaty, quatz, quatw);
  tf::Matrix3x3 m(q);
  double roll, pitch, yaw;
  m.getRPY(roll, pitch, yaw);
  this->theta = -yaw;
}

void Statemachine::teamfarbeCallback(const std_msgs::Bool::ConstPtr& done){
  if(done->data){
    teamfarbeDone = true;
  }else{
    teamfarbeDone = false;
  }
}

void Statemachine::abDoneCallback(const std_msgs::Bool::ConstPtr& done){
  if(done->data){
    abDone = true;
  }else{
    abDone = false;
  }
}
void Statemachine::goalCallback(const std_msgs::Bool::ConstPtr& reached){
  if(reached->data && state >= 2){
    reachedGoal = true;
    ROS_INFO("reached: [%d]", reachedGoal);
  }else{
    reachedGoal = false;
  }

}

void Statemachine::startCallback(const std_msgs::Bool::ConstPtr& go){
  if(go->data){
    start = true;
  }
}

int main(int argc, char** argv){
  int count = 0;
  ros::init(argc, argv, "state_machine");
  Statemachine statemachine = Statemachine();
  ros::Rate rate(10);
  while(ros::ok()){
    ROS_INFO("State:[%i]", statemachine.state);
      // Start game
    if(statemachine.start && statemachine.state == 0){
        statemachine.state = 1;
    }else if(!statemachine.start){
        statemachine.state = 0;
    }
    // Check for current state
    if(statemachine.state == 1){
      statemachine.OD_mode.data = 1;
      if(statemachine.abDone){
        statemachine.OD_mode.data = 0;
        statemachine.state = 2;
        // Setze Spielfeld grenzen
      }
    }
    if(statemachine.state == 2){
      if(count == 0){
        statemachine.goal.position.x = 0;
        statemachine.goal.position.y = -1;
        statemachine.goalPub.publish(statemachine.goal);
        ros::Duration(2).sleep();
      }
      ROS_INFO("count: [%i]", count);
      if(statemachine.reachedGoal){
        if(count < 3){
          count++;
        }
        ROS_INFO("count: [%i]", count);
        if(count == 1){
          statemachine.goal.position.x = 0;
          statemachine.goal.position.y = -0.9;
          statemachine.goalPub.publish(statemachine.goal);
          ros::Duration(2).sleep();
          //statemachine.reachedGoal = false;
        }

        if(count == 2){
          ROS_INFO("bin rein gesprungen");
          statemachine.IP_mode.data = 1;
          count = 0;
        }
      }
      if(statemachine.teamfarbeDone){
        statemachine.state = 3;
      }
    }
    if(statemachine.state == 3){
      statemachine.goal.position.x = 0;
      statemachine.goal.position.y = 0;
      statemachine.goalPub.publish(statemachine.goal);
      ros::Duration(2).sleep();
      if(statemachine.reachedGoal){
        count++;
        statemachine.goal.position.x = 0.01;
        statemachine.goal.position.y = -0.015;
        statemachine.goalPub.publish(statemachine.goal);
        ros::Duration(2).sleep();
        if(count == 2){
          statemachine.state = 4;
          count = 0;
        }
      }
    }
    if(statemachine.state == 4){
      if(statemachine.startGame){
        statemachine.state = 5;
      }
    }
    if(statemachine.state == 5){
      statemachine.IP_mode.data = 2;
      statemachine.OD_mode.data = 2;
      statemachine.IPPub.publish(statemachine.IP_mode);
      statemachine.ODPub.publish(statemachine.OD_mode);
      ros::Duration(3).sleep();
      statemachine.goal.position.x = 0.01;
      statemachine.goal.position.y = 0.01;
      statemachine.goalPub.publish(statemachine.goal);
      ros::Duration(2).sleep();
      if(fabsf(statemachine.theta_vel) < 0.01){
        count++;
        statemachine.goal.position.x = 0.01;
        statemachine.goal.position.y = -0.01;
        statemachine.goalPub.publish(statemachine.goal);
        ros::Duration(2).sleep();
        if(count == 2){
          statemachine.state = 6;
        }
      }
    if(statemachine.state == 6){
      statemachine.goalPub.publish(statemachine.goal);
      ros::Duration(2).sleep();
      if(statemachine.reachedGoal && statemachine.puckFound){
        statemachine.state = 7;
      }
    }
    if(statemachine.state == 7){
      statemachine.goal.position.x = 2.4;
      statemachine.goal.position.y = 0;
      statemachine.goalPub.publish(statemachine.goal);
      ros::Duration(2).sleep();
      if(statemachine.reachedGoal){
        ROS_INFO("GOALLLLLL!!!!!!!!");
        statemachine.shotGoal.data = true;
        statemachine.shotGoalPub.publish(statemachine.shotGoal);
        ros::Duration(1).sleep();
        statemachine.shotGoal.data = false;
        statemachine.shotGoalPub.publish(statemachine.shotGoal);
        statemachine.state = 8;
      }
    }
    if(statemachine.state == 8){
    }



    }
    statemachine.goalPub.publish(statemachine.goal);
    statemachine.ODPub.publish(statemachine.OD_mode);
    statemachine.IPPub.publish(statemachine.IP_mode);
    statemachine.isAlive.data = true;
    statemachine.isAlivePub.publish(statemachine.isAlive);

    ros::spinOnce();               // check for incoming messages
    rate.sleep();
  }
}
