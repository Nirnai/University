/***************************************************************************
 *   Server GUI "angelina" for "Projektkurs C++" WS 07/08				   *
 *   																	   *
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

#ifndef TESTGUI_H
#define TESTGUI_H

#include <QLineEdit>
#include <QListWidget>
#include <QRadioButton>
#include <QTimer>
#include <QWidget>
#include <QMap>

#include <std_msgs/Bool.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Bool.h>
#include <std_msgs/String.h>

#include <QString>
#include "referee.h"


class TestGui: public QWidget
{
	Q_OBJECT

	public:
		TestGui(QWidget *parent = 0);
		std_msgs::Bool stopImmediately;
        std_msgs::Bool startInitPhase;
        std_msgs::Bool startGamePhase;
        std_msgs::String correctTeamColor;
        std_msgs::Float32 correctABRatio;
        std_msgs::Float32 correctA;
        std_msgs::Float32 correctB;

		bool sendRdyCorrect;
        bool initialisierungABReceived;
        bool initialisierungTeamColorReceived;


	private:
		int myPosID;
		int myOriID;
		QLineEdit *edID;
		QLineEdit *edPort;
		QLineEdit *edIP;
		QLineEdit *edRatio;
		QLineEdit *edX;
		QLineEdit *edY;
		QListWidget *list;
		QRadioButton *yellowRbtn;
		QRadioButton *blueRbtn;
		QTimer *aliveTimer;
		Referee *referee;

		TeamColor finalTeamColor;

		QString ip;
		int port;
		int id;

		//int colorGetAngelina = 0; //1 = Farbe gerade bekommen


	private Q_SLOTS:
		void slotConnect();
		void slotReportReady();
		void slotReportDone();
		void slotReportGoal();
		void slotSendAlive();
		void slotToggleAliveTimer(bool);
		void slotTellAbRatio();
		void slotTellEgoPos();
		void slotGameStart();
		void slotDetectionStart();
		void slotGameOver();
		void slotAbValues(double a, double b);
		void slotStopMovement();
		void slotTellTeamColor();
		void slotTeamColor(TeamColor);


	public:
		void connectReferee();
		void sendTeamColor(int color);
		void sendRatio(double a, double b);
		void sendInit(int color, double a, double b);
		void tellPosition(float x, float y);
		void reportRdy();
		void sendAlive();
        void sendABRatio(double ratio);
        void sendTeamColor(std_msgs::String color);
};


#endif /* TESTGUI_H */
