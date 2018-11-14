#include <controller_interface/controller.h>
#include <hardware_interface/force_torque_sensor_interface.h>
#include <pluginlib/class_list_macros.h>

#include <geometry_msgs/Wrench.h>
#include <geometry_msgs/Vector3Stamped.h>
#include <visualization_msgs/Marker.h>
#include <tf/transform_listener.h>


namespace reemc_zmp_estimator{


  class ZmpController : public controller_interface::Controller<hardware_interface::ForceTorqueSensorInterface>
  {
  public:
    bool init(hardware_interface::ForceTorqueSensorInterface* hw, ros::NodeHandle &n)
    {
      // zmp_L_pub_ = n.advertise<geometry_msgs::Point>("zmp_L", 1000);
      // marker_L_pub_ = n.advertise<visualization_msgs::Marker>("zmp_L_marker", 10);
      // marker_R_pub_ = n.advertise<visualization_msgs::Marker>("zmp_R_marker", 10);
      // get force torque sensors names
      const std::vector<std::string>& sensor_names = hw->getNames();
      for (unsigned i=0; i<sensor_names.size(); i++)
        ROS_INFO("Got sensor %s", sensor_names[i].c_str());

      // Store left ankle ft sensor
      sensors_.push_back(hw->getHandle(sensor_names[0]));
      // Store right ankle ft sensor
      sensors_.push_back(hw->getHandle(sensor_names[2]));

      for(unsigned int s = 0; s<2; ++s){
        ROS_INFO_STREAM( "Sensor Names: " << sensors_[s].getName() << "Index: " << s);
      }



      return true;
    }

    void update(const ros::Time& time, const ros::Duration& period)
    {
     
      zmp.header.stamp = zmp_L.header.stamp = zmp_R.header.stamp = ros::Time::now();
      zmp.header.frame_id = "/base_link";
      zmp_L.header.frame_id = "/left_sole_link";
      zmp_R.header.frame_id = "/right_sole_link";

      // Get current sensor readings
      const double* torque_L = sensors_[0].getTorque();
      const double* torque_R = sensors_[1].getTorque();
      const double* force_L = sensors_[0].getForce();
      const double* force_R = sensors_[1].getForce();

      // Compute zmp for each foot
      zmp_L.vector.x = - torque_L[1]/force_L[2];
      zmp_L.vector.y =   torque_L[0]/force_L[2];
      zmp_L.vector.z =                        0;

      zmp_R.vector.x = - torque_R[1]/force_R[2];
      zmp_R.vector.y =   torque_R[0]/force_R[2];
      zmp_R.vector.z =                        0;

      // Transform to reference frame
      geometry_msgs::Vector3Stamped tmp_L, tmp_R;

      tf::StampedTransform transform;
      try{
        listener.lookupTransform("/base_link", "/left_sole_link",  
                                 ros::Time(0), transform);
      }
      catch (tf::TransformException ex){
        ROS_ERROR("%s",ex.what());
      }

      // try
      // {
      // listener.waitForTransform("/base_link", "/left_sole_link", zmp.header.stamp, ros::Duration(3.0));
      // listener.waitForTransform("/base_link", "/right_sole_link", zmp.header.stamp, ros::Duration(3.0));
      listner.transformVector("/base_link", zmp_L, tmp_L);
      listener.transformVector("/base_link", zmp_R, tmp_R);
      // } 
      // catch (tf::TransformException ex)
      // {
      //   ROS_ERROR("%s",ex.what());
      //   ros::Duration(1.0).sleep();
      // }



      ROS_INFO_STREAM_THROTTLE(1.0, "sole: [" << zmp_L.vector.x << "," << zmp_L.vector.y << "," << zmp_L.vector.z << "]\n" <<
                                    "base: [" << tmp_L.vector.x << "," << tmp_L.vector.y << "," << tmp_L.vector.z << "]");

    
      // visualization_msgs::Marker zmp_L_marker, zmp_R_marker;

      // zmp_L_marker.header.stamp = zmp_R_marker.header.stamp = ros::Time::now();

      // zmp_L_marker.type = zmp_R_marker.type = visualization_msgs::Marker::SPHERE;
      // zmp_L_marker.action = zmp_R_marker.action = visualization_msgs::Marker::ADD;
      // zmp_L_marker.pose.orientation.w = zmp_R_marker.pose.orientation.w = 1.0;

      // zmp_L_marker.scale.x = zmp_R_marker.scale.x = 0.05;
      // zmp_L_marker.scale.y = zmp_R_marker.scale.y = 0.05;
      // zmp_L_marker.scale.z = zmp_R_marker.scale.z = 0.001;
      // zmp_L_marker.color.a = zmp_R_marker.color.a = 1.0;
      // zmp_L_marker.color.r = zmp_R_marker.color.r = 0.0;
      // zmp_L_marker.color.g = zmp_R_marker.color.g = 1.0;
      // zmp_L_marker.color.b = zmp_R_marker.color.b = 0.0;

      // zmp_L_marker.header.frame_id = "/left_sole_link";
      // zmp_R_marker.header.frame_id = "/right_sole_link";

      // zmp_L_marker.pose.position.x = zmp_L.x;
      // zmp_L_marker.pose.position.y = zmp_L.y;
      // zmp_L_marker.pose.position.z = 0;

      // zmp_R_marker.pose.position.x = zmp_R.x;
      // zmp_R_marker.pose.position.y = zmp_R.y;
      // zmp_R_marker.pose.position.z = 0;




      // ROS_INFO_STREAM_THROTTLE(1.0, "ZMP_L x = "<<zmp_L.x<<" ZMP_L y = "<< zmp_L.y);
      // // zmp_pub_.publish(zmp);
      // marker_L_pub_.publish(zmp_L_marker);
      // marker_R_pub_.publish(zmp_R_marker);



    }

    void starting(const ros::Time& time) { }
    void stopping(const ros::Time& time) { }

  private:
    std::vector<hardware_interface::ForceTorqueSensorHandle> sensors_;
    geometry_msgs::Vector3Stamped zmp;
    geometry_msgs::Vector3Stamped zmp_L; 
    geometry_msgs::Vector3Stamped zmp_R;
    tf::TransformListener listener;

    // ros::Publisher zmp_L_pub_;
    // ros::Publisher marker_L_pub_;
    // ros::Publisher marker_R_pub_;
    
    

  };
  PLUGINLIB_DECLARE_CLASS(reemc_zmp_estimator, ZmpController, reemc_zmp_estimator::ZmpController, controller_interface::ControllerBase);

}
