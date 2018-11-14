#ifndef UR_ROBOT_LLI_SIMPLEEFFORTCONTROL_H
#define UR_ROBOT_LLI_SIMPLEEFFORTCONTROL_H

#include<tum_ics_ur_robot_lli/RobotControllers/ControlEffort.h>
#include<imu_project/DesiredState.h>
#include <tf/transform_broadcaster.h>
#include "tf_conversions/tf_eigen.h"
#include <tf/transform_datatypes.h>
#include <math.h>
#include <Eigen/QR>


namespace tum_ics_ur_robot_lli{
namespace RobotControllers{

class SimpleEffortControl: public ControlEffort
{
    // Member variables
private:

    bool m_startFlag;

    //////////////////
    bool startPosition;
    bool rememberR;
    //////////////////

    Vector6d m_qStart;
    JointState m_qInit;
    JointState m_qHome;
    JointState m_qPark;

    ros::NodeHandle n;
    ros::Publisher pubCtrlData;
    ros::Subscriber desiredStateSub;

    Matrix6d m_Kp;
    Matrix6d m_Kd;
    Matrix6d m_Ki;
    Vector6d m_goal;
    double m_totalTime;

    Vector6d m_DeltaQ;
    Vector6d m_DeltaQp;
    Vector6d m_iDeltaQ;

    ///////////////////////////
    Matrix6d JA;
    Matrix6d JAp;
    Matrix6d JDamped;

    Vector6d X;
    // Rotion Base 
    Matrix3d R_old;
    Matrix3d R;

    Matrix3d R_home;
    Matrix3d R_offset2;
    Vector3d euler_offset;

    Vector6d Xp;
    Vector6d Xpp;

    Vector6d Xd;
    Vector6d Xdp;
    Vector6d Xdpp;

    Vector6d DeltaX;
    Vector6d DeltaXp;

    double old_time;

    Matrix<double,95, 1> Theta_init;
    Matrix<double,95, 1> Theta_estimate;
    Matrix<double,6,95> Yr;

    // MatrixXd Gamma(95,95);

    Vector6d filtered_tau;
    Vector6d tau_original;
    int count;

    ///////////////////////////

    double m_controlPeriod;     //[s]
    double m_controlPeriod_2;   //[s]

    // Member Methods
public:
    SimpleEffortControl(double weight=1.0,
                        const QString& name="SimpleEffortCtrl");
    ~SimpleEffortControl();

    void setQInit(const JointState& qinit);
    void setQHome(const JointState& qhome);
    void setQPark(const JointState& qpark);

    ////////////////////////////////////////////////////////
    Vector3d rot2eul(Matrix3d &R);
    void getJacobiAnalytical(const JointState &current);
    void getJacobiAnalyticalp(const JointState &current);
    void getJacobiDamped();
    void getForwardKinematics(const JointState &current);
    void getDifferentialKinematics(const JointState &current);
    void getDesiredPose(const imu_project::DesiredState::ConstPtr &state);
    void dynamicModel(const JointState &current, const JointState &ref);
    void updateParameters(Vector6d &Sq, Matrix<double,95, 1> &Theta_estimate);
    void EMA(Vector6d &filtered_tau, Vector6d &tau_original, float alpha, int count);
    ///////////////////////////////////////////////////////


    // Memeber Methods
private:
    bool init();
    bool start();
    Vector6d update(const RobotTime& time, const JointState &current);
    bool stop();


};

}
}



#endif // UR_ROBOT_LLI_SIMPLEEFFORTCONTROL_H
