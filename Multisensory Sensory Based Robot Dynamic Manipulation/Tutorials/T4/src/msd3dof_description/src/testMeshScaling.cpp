/*********************************************************************
* STD INCLUDES
********************************************************************/
#include <iostream>
#include <fstream>
#include <pthread.h>


/*********************************************************************
* ROS INCLUDES
********************************************************************/
#include <ros/ros.h>
#include <std_msgs/String.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <tf/transform_listener.h>
    #include <tf/transform_broadcaster.h>
#include <tf_conversions/tf_eigen.h>
#include <sensor_msgs/JointState.h>

/*********************************************************************
* EIGEN INCLUDES
********************************************************************/
#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <Eigen/Geometry>
#include <Eigen/Core>

using namespace Eigen;

double gPose;

void springScale(const sensor_msgs::JointState::ConstPtr& msg)
{
  //Get the position of joint 1
   
  gPose = msg->position[0];
  //ROS_INFO_STREAM("I heard: ["<<gPose<<"]");
}

const visualization_msgs::Marker createMarkerMesh(std::string frame_id, int id, int shape,
                                                  double x, double y, double z, /*position*/
                                                  double q_w, double q_x, double q_y, double q_z, /*orientation in quatern*/
                                                  double s_x, double s_y, double s_z, std::string meshFile/*scale*/)
{
    visualization_msgs::Marker marker;

    marker.header.frame_id = frame_id;
    marker.header.stamp = ros::Time();
    marker.ns = "tracker_markers";
    marker.id = id;
    marker.type = shape;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.position.x = x;
    marker.pose.position.y = y;
    marker.pose.position.z = z;
    marker.pose.orientation.x = q_x;
    marker.pose.orientation.y = q_y;
    marker.pose.orientation.z = q_z;
    marker.pose.orientation.w = q_w;
    marker.scale.x = s_x;
    marker.scale.y = s_y;
    marker.scale.z = s_z;
    marker.mesh_resource = meshFile;

    marker.color.r = 1.0;
    marker.color.g = 0.0;
    marker.color.b = 0.3;
    marker.color.a = 1;
    marker.lifetime = ros::Duration();

    return marker;
}



int main( int argc, char** argv )
{

    ros::init(argc, argv, "spring_visualization",ros::init_options::AnonymousName);

    ROS_INFO_STREAM("**Publishing spring for rviz..");

    ros::NodeHandle n;
    ros::Rate r(60);

    static tf::TransformBroadcaster br;


    tf::Transform transform;

    ros::Publisher marker_pub = n.advertise<visualization_msgs::Marker>("visualization_marker", 1);



    visualization_msgs::Marker msd3DOF;

	//Create a marker using the function createMarkerMesh (see above), you need to load the spring.stl CAD file. 
    //Connect this marker to the tf "/spring"
    msd3DOF = createMarkerMesh("/spring",0,10,0,0,0,0,0,-1,0,0.3,0.3,0.3,"package://msd3dof_description/meshes/spring.stl");


    //Subscribe to the joint_state topic using the function springScale (see above) as a callback function
    ros::Subscriber sub = n.subscribe("/joint_states",10,&springScale);

    tf::Quaternion qtf;


    while(ros::ok())
    {

        //Spring
        transform.setOrigin( tf::Vector3(0.0, 0.0, 0.0) );
        tf::Quaternion q;
        q.setRPY(0, 0, 0);
        transform.setRotation(q);
        //publish  a tf named "/spring" with respect to the Link_0
        br.sendTransform(tf::StampedTransform(transform, ros::Time(),"link_0", "/spring"));
    
        //Define the scale of the spring with the position of joint 1 
        msd3DOF.scale.z = 0.3 + gPose;

        //Publish the marker
        marker_pub.publish(msd3DOF);

        ros::spinOnce();

        r.sleep();
    }
    //////////////////////////////////Skin Cells Visualization//////////////////////////////////
    return 0;
}

