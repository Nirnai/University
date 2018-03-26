/***************************************************************************
 *   Server GUI "angelina" for "Projektkurs C++" WS 07/08		   *
 *   									   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <QtGui>
#include "testgui.h"
#include "ros/ros.h"
#include "std_msgs/String.h"
#include <sstream>
#include <geometry_msgs/Pose2D.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Bool.h>
#include <nav_msgs/Odometry.h>


bool newPositionReceived;
bool receivedABCallback;
bool receivedTeamColorCallback;
bool aBsendingDone;
bool teamColorSendingDone;
geometry_msgs::Pose2D currentPosition;
std_msgs::Bool abCalcRdyFlag;
std_msgs::Bool ipAliveFlag;
std_msgs::Bool odAliveFlag;
std_msgs::Bool teamColorRdyFlag;
std_msgs::Float32 abRatio;
std_msgs::Float32 a;
std_msgs::Float32 b;
std_msgs::String teamColor;

void ownPositionCallback(const geometry_msgs::Pose2D::ConstPtr &msg)
{
    // TODO: AUS ODOM EINBAUEN Nirnai statemachine zeile 31
    newPositionReceived = true;
    currentPosition.x = msg->x;
    currentPosition.y = msg->y;
    currentPosition.theta = msg->theta;
}
void abCalcRdyCallback(const std_msgs::Bool::ConstPtr &msg)
{
    abCalcRdyFlag.data = msg->data;
}
void ipAliveCallback(const std_msgs::Bool::ConstPtr &msg)
{
    ipAliveFlag.data = msg->data;
}
void odAliveCallback(const std_msgs::Bool::ConstPtr &msg)
{
    odAliveFlag.data = msg->data;
}
void abRatioCallback(const std_msgs::Float32::ConstPtr &msg)
{
    receivedABCallback = true;
    abRatio.data = msg->data;
}
void odomCallback(const nav_msgs::Odometry::ConstPtr& odom){

  currentPosition.x = odom->pose.pose.position.x;
  currentPosition.y = odom->pose.pose.position.y;
}

//void teamColorRdyCallback(const std_msgs::Bool::ConstPtr &msg)
//{
//    teamColorRdyFlag.data = msg->data;
//}ros::SubscaBsendingDoneriber
void teamColorCallback(const std_msgs::String::ConstPtr &msg)
{
    teamColor.data = msg->data;
    receivedTeamColorCallback = true;
}

int main(int argc, char** argv)
{
	QApplication app(argc, argv);
	TestGui gui;
	gui.show();
	TestGui gui2;
	gui2.show();
	// return app.exec();

    abCalcRdyFlag.data = false;
    aBsendingDone = false;
    teamColorSendingDone = false;
    receivedTeamColorCallback = false;
    abRatio.data = 0.0;

    bool counterSet = false;

    ros::init(argc, argv, "communication");
    ros::NodeHandle node;
    ros::Publisher stopImmediatelyPub = node.advertise<std_msgs::Bool>("/communication/stopImmediately", 10);
    ros::Publisher startInitPub = node.advertise<std_msgs::Bool>("/communication/startInit", 10);
    ros::Publisher startGamePub = node.advertise<std_msgs::Bool>("/communication/startGame", 10);
    ros::Publisher correctTeamColorPub = node.advertise<std_msgs::String>("/communication/correctTeamColor", 10);
    ros::Publisher correctABRatioPub = node.advertise<std_msgs::Float32>("/communication/correctABRatio", 10);
    ros::Publisher correctAPub = node.advertise<std_msgs::Float32>("/communication/correctA", 10);
    ros::Publisher correctBPub = node.advertise<std_msgs::Float32>("/communication/correctB", 10);

    ros::Subscriber ipAliveSub    = node.subscribe("/image_processing/isAlive", 1, ipAliveCallback);
    ros::Subscriber odAliveSub    = node.subscribe("/object_detection/isAlive", 1, odAliveCallback);
    ros::Subscriber odomSub       = node.subscribe("/odom", 1000, odomCallback);
    ros::Subscriber tellPositionSub = node.subscribe("/localization/ownPosition", 1, ownPositionCallback);
    ros::Subscriber rdyAbCalcSub = node.subscribe("/object_detection/initABSuccess", 1, abCalcRdyCallback);
    ros::Subscriber getABRatioSub = node.subscribe("/object_detection/abRatio", 1, abRatioCallback);
//    ros::Subscriber rdyTeamColorCalcSub = node.subscribe("/image_processing/xxxx", 1, teamColorRdyCallback);
    ros::Subscriber teamColorSub = node.subscribe("/image_processing/teamColor", 1, teamColorCallback);

     // send from angelina
	// ros::Publisher connectedPub_ = nh_.advertise<std_msgs::Bool>("/communication/connected", 1);                                      // publish connection status
	// ros::Publisher stopPub_ = nh_.advertise<std_msgs::Bool>("/communication/stop_signal", 1);                                        // publish stop command
	// ros::Publisher gameOverPub_ = nh_.advertise<std_msgs::Bool>("communication/game_over_signal", 1);                               // publish game over command


    ros::Rate rate(30);
    double beginTime = ros::Time::now().toSec();
    double currentTimeAB = ros::Time::now().toSec();
    double currentTimeAB2 = ros::Time::now().toSec();
    double currentTime = ros::Time::now().toSec(); double oldTime = currentTime;
    double oldTime2 = currentTime; double currentTime2 = currentTime;
    while (ros::ok()){
        app.processEvents();
        if (gui.startInitPhase.data)
        {
            if(!counterSet)
            {
                currentTimeAB = ros::Time::now().toSec();
                counterSet = true;
            }
            currentTimeAB2 = ros::Time::now().toSec();
        }
        // ------------------------------ Start --------------------------------------------------

        if( ipAliveFlag.data && odAliveFlag.data )
    		if(!gui.sendRdyCorrect)
    		{
    			gui.reportRdy();
    		}

        // ------------------------------ Initialisierung ----------------------------------------

        // Sending --------
        // Check if AB Calculation is rdy

        // ROS_INFO("Bool: %d, Time: %f, %f", aBsendingDone, currentTimeAB, currentTimeAB2);
        if( (abCalcRdyFlag.data && !aBsendingDone) || ((abs(currentTimeAB2-currentTimeAB) > 20) && !aBsendingDone))
        {
            // Calculation seems to be rdy, sent AB to Angelina
            if(receivedABCallback)
            {
                gui.sendABRatio(abRatio.data);
                aBsendingDone = true;
            }
        }
        if( receivedTeamColorCallback && !teamColorSendingDone )
        {
            // send team Color
            gui.sendTeamColor( teamColor);
            teamColorSendingDone = true;
        }
        // Receiving --------
        if( gui.initialisierungTeamColorReceived )
        {
            // Publishen!!!
            correctTeamColorPub.publish(gui.correctTeamColor);
        }
        if( gui.initialisierungABReceived )
        {
            correctABRatioPub.publish(gui.correctABRatio);
            correctAPub.publish(gui.correctA);
            correctBPub.publish(gui.correctB);
        }
        // -------------------------- kontinuierlicher Teil --------------------------------------
		currentTime = ros::Time::now().toSec();
        if( (currentTime - oldTime) > 30.0)
		{
            gui.sendAlive();
			oldTime = currentTime;
		}

        currentTime2 = ros::Time::now().toSec();
        if((currentTime2 - oldTime2) > 10.0)
        {
            if( gui.initialisierungABReceived )
            {
                if( gui.correctTeamColor.data == "yellow" )
                {
                    currentPosition.x = currentPosition.x + 2.75*gui.correctA.data;
                    currentPosition.y = currentPosition.y + gui.correctB.data/2.0;
                    gui.tellPosition(currentPosition.x,currentPosition.y);
                    oldTime2 = currentTime2;
                    // newPositionReceived = false;
                }
                else if( gui.correctTeamColor.data == "blue" )
                {
                    currentPosition.x = currentPosition.x + 3/8*gui.correctA.data;
                    currentPosition.y = currentPosition.y + gui.correctB.data/2.0;
                }
            }
        }
		// ROS_INFO("current: %f, old: %f", currentTime, oldTime );

        // Publish data

        stopImmediatelyPub.publish(gui.stopImmediately);
        startGamePub.publish(gui.startGamePhase);
        startInitPub.publish(gui.startInitPhase);
        ros::spinOnce();
        rate.sleep();
    }

    return app.exec(); // the program will keep running
}
