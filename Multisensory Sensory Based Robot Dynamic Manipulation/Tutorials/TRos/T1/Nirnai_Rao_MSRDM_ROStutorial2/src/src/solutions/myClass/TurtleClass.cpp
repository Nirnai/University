#include<iostream>
#include<turtle_vis/myClass/TurtleClass.h>
#include<Eigen/Dense>


namespace turtleSpace {

TurtleClass::TurtleClass()
//: turtlePose_g(turtlePose_g(0,0,0))
{
    count_mutex = PTHREAD_MUTEX_INITIALIZER;
}
TurtleClass::~TurtleClass()
{
    turtlePose_g << 1,0,0;
}


void TurtleClass::getPose(const turtle_vis::DesiredPose::ConstPtr &msg)
{
    pthread_mutex_lock( &this->count_mutex );
    //#>>>>TODO: COPY THE MSG TO A LOCAL VARIABLE
    Vector3d pose(msg->x, msg->y, msg->theta);
    turtlePose_g = pose;
    pthread_mutex_unlock( &this->count_mutex );

    //#>>>>TODO:PLOT THE OBTAINED DATA
    ROS_INFO_STREAM("x = "<< pose(0));
    ROS_INFO_STREAM("x = "<< pose(1));
    ROS_INFO_STREAM("x = "<< pose(2));
}


bool TurtleClass::getDPose(turtle_vis::send_desired_pose::Request &req, turtle_vis::send_desired_pose::Response &res)
{
    pthread_mutex_lock( &this->count_mutex );
    //#>>>>TODO:COPY THE REQUEST MSG TO A LOCAL VARIABLE
    Vector3d dpose(req.x, req.y, req.theta);

    turtlePose_desired_g = dpose;

    pthread_mutex_unlock( &this->count_mutex );

    //#>>>>TODO:PLOT THE OBTAINED DATA
    ROS_INFO_STREAM("x = "<< dpose(0));
    ROS_INFO_STREAM("y = "<< dpose(1));
    ROS_INFO_STREAM("theta = "<< dpose(2));
    res.reply=1;

    return true;
}


Vector3d TurtleClass::getLocalPose()
{
    Vector3d local;
    pthread_mutex_lock( &this->count_mutex );
    local=this->turtlePose_g;
    pthread_mutex_unlock( &this->count_mutex );

    return local;
}

Vector3d TurtleClass::getLocalDesiredPose()
{
    Vector3d local;
    pthread_mutex_lock( &this->count_mutex );
    local=this->turtlePose_desired_g;
    pthread_mutex_unlock( &this->count_mutex );

    return local;
}



}
