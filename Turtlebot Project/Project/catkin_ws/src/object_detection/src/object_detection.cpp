#include <geometry_msgs/Twist.h>
#include "object_detection.h"

ObjectDetection::ObjectDetection() {

  // Subscribe to the simulated robot's laser scan topic
    laserSub = node.subscribe("/lidarscan", 1, &ObjectDetection::scanCallback, this);
//    gmapSub = node.subscribe("/map", 1, &ObjectDetection::mapCallback, this);
    alivePub = node.advertise<std_msgs::Bool>("/object_detection/isAlive",10);
    modiSub = node.subscribe("/object_detection/mode", 1, &ObjectDetection::modeCallback, this);
    abRatioPub = node.advertise<std_msgs::Float32>("/object_detection/abRatio",10);
    ownPositionPub = node.advertise<geometry_msgs::Pose>("/object_detection/ownPosition",10);
    ownPosition2DPub = node.advertise<geometry_msgs::Pose2D>("/object_detection/ownPosition_2D",10);
    enemyGoalPub = node.advertise<object_detection::goal_edge_points>("/object_detection/enemyGoalPosition",10);
    gameFieldPub = node.advertise<object_detection::game_field>("/object_detection/game_field",10);
    obstaclePosesPub = node.advertise<object_detection::obstacle_poses>("/object_detection/obstacle_poses",10);
    initABSuccessPub = node.advertise<std_msgs::Bool>("/object_detection/initABSuccess",10);

    isAlive.data = false;
    abRatio.data = 0.0;
    initABOkay.data = false;

    receivecFirstLidar = false;
    receivecFirstMode = false;

    foundSomethingOnLeftSide = false;
    foundSomethingOnRightSide = false;

//    measurements = 0;
    outputMode = 0; // 0 = none, 1 = some, 2 = all
    rangeMultMin = 1.4; // Scan radius multiplicator from min
    objectScanOffset = 0.15; // Distance difference threshold where lidar points are counted to one object
    indexOffsetToObjRow1 = 0; // index offset to adress distances between obj correct (workaround)
    indexOffsetToObjRow2 = 0;
    measurements = 0;
    measurementsDone = false;

    meanRA = 0.0;
    varRA = 0.0;
    medianRA = 0.0;
    meanLA = 0.0;
    varLA = 0.0;
    medianLA = 0.0;
    meanB = 0.0;
    varB = 0.0;
    medianB = 0.0;

    leftIndexForB = 0;
    rightIndexForB = 0;

    begin = ros::Time::now().toSec();
    runtime;
    diff;
}

// Process the incoming laser scan message
void ObjectDetection::scanCallback(const sensor_msgs::LaserScan::ConstPtr& scan)
{
    lidarScan = scan;
    if (!receivecFirstLidar)
        receivecFirstLidar = true;
}
void ObjectDetection::mapCallback(const nav_msgs::OccupancyGridConstPtr &map)
{
//    geometry_msgs/Pose origin
//    geometry_msgs::Pose origin = map->origin;
//    geometry_msgs::Point position = origin.position;
    nav_msgs::MapMetaData info = map->info;

    geometry_msgs::Pose origin = info.origin;
    geometry_msgs::Point position = origin.position;

    ROS_INFO("Map received: %f, %f", position.x, position.y);

    for (int i=0; i<4000*4000; ++i)
    {
//        if( map->data[i] != -1)
//            ROS_INFO("KEIN MINUS EINS");
    }

    ROS_INFO("done");
//        ROS_INFO("Wert: %i", map->data[i]);

    //    geometry_msgs/Pose origin = map->origin;
}
void ObjectDetection::modeCallback(const std_msgs::ByteConstPtr &mode)
{
    currentCalcMode = mode;
    if(!receivecFirstMode)
        receivecFirstMode = true;

}
geometry_msgs::Quaternion ObjectDetection::toQuaternion(float pitch, float roll, float yaw)
{
    geometry_msgs::Quaternion q;
    // Abbreviations for the various angular functions
    double cy = cos(yaw * 0.5);
    double sy = sin(yaw * 0.5);
    double cr = cos(roll * 0.5);
    double sr = sin(roll * 0.5);
    double cp = cos(pitch * 0.5);
    double sp = sin(pitch * 0.5);

    q.w = cy * cr * cp + sy * sr * sp;
    q.x = cy * sr * cp - sy * cr * sp;
    q.y = cy * cr * sp + sy * sr * cp;
    q.z = sy * cr * cp - cy * sr * sp;
    return q;
}

void ObjectDetection::scanForObjects(float scanInDist)
{
    if(!lidarScan)
    {
      return;
    }

    // float angleRange = lidarScan->angle_max - lidarScan->angle_min;
    int numberScansCompl = ceil((lidarScan->angle_max - lidarScan->angle_min) / lidarScan->angle_increment);
    // int minIndex = ceil((minAngle/180*M_PI - lidarScan->angle_min) / lidarScan->angle_increment);
    // int maxIndex = floor((maxAngle/180*M_PI - lidarScan->angle_min) / lidarScan->angle_increment);
    // int numberScansEval = maxIndex - minIndex;
    // bool distances[numberScansEval] = { };
    int tmp_ctr=1;
    float scanDistances[numberScansCompl];
    int counterForSaving = 0;
    float meanDist = 0.0;
    float angle_mid;

    for(int l = 0; l < numberScansCompl; ++l)
    {
        if(lidarScan->ranges[l] == INFINITY)
            scanDistances[counterForSaving] = 20.0;
        else
            scanDistances[counterForSaving] = lidarScan->ranges[l];

        counterForSaving++;
    }
    // for( int haha=1; haha < numberScansCompl; ++haha)
    // {
    //     std::cout << scanDistances[haha] << " ";
    // }
    // std::cout << "\n";
    // ------------------------- Search for objects ---------------------------------------------

    int loop_counter = 0;
    int obs_counter = 0;
    while ( loop_counter < numberScansCompl)
    {
        bool nextIsOk = true;
        if (scanDistances[loop_counter] < scanInDist)
        {
            tmp_ctr=1;
            nextIsOk = true;
            while(nextIsOk)
            {
                if( abs(scanDistances[loop_counter] - scanDistances[loop_counter+tmp_ctr]) < 0.2)
                    tmp_ctr++; // Punkt gehört dazu
                else
                    nextIsOk = false;
            }
            if(tmp_ctr >= 1)
            {
                // save obstacle
                ObstacleInLidar newObstacle;
                geometry_msgs::Pose2D pose_of_obstacle;
                newObstacle.index_min = loop_counter;
                newObstacle.index_max = loop_counter + tmp_ctr;
                newObstacle.angle_min = 360.0*newObstacle.index_min/numberScansCompl;
                newObstacle.angle_max = 360.0*newObstacle.index_max/numberScansCompl;
                newObstacle.angle_mean = (newObstacle.angle_min + newObstacle.angle_max)/2;
                angle_mid = (newObstacle.angle_max-newObstacle.angle_min)/2;
                for(int a=0;a<tmp_ctr;++a)
                {
                    meanDist = meanDist + scanDistances[loop_counter+a];
                }
                newObstacle.distance = meanDist/tmp_ctr;
                newObstacle.x = newObstacle.distance*cos(newObstacle.angle_mean*M_PI/180.0);
                newObstacle.y = newObstacle.distance*sin(newObstacle.angle_mean*M_PI/180.0);
                obstacleStackCompl.push_back(newObstacle);

                pose_of_obstacle.x = -newObstacle.x;
                pose_of_obstacle.y = newObstacle.y;
                pose_of_obstacle.theta = newObstacle.angle_mean;
                obsPoses.obstacles.push_back(pose_of_obstacle);
                obs_counter++;

                meanDist = 0.0;
            }
        }
        loop_counter = tmp_ctr + loop_counter;
        tmp_ctr = 1;
    }
    // if (outputMode == 1 || outputMode == 2)
    //     ROS_INFO("Gefundene Objekte: %i",(int) obstacleStack.size());
    // if (outputMode == 1 || outputMode == 2)
    // {
    //     for(int ind = 0; ind < obstacleStackCompl.size(); ++ind)
    //     {
    //         std::cout << " Objekt: " << ind << " Angle: " << obstacleStackCompl[ind].angle_mean << " Distanz: " << obstacleStackCompl[ind].distance
    //                     << " Index min: " << obstacleStackCompl[ind].index_min << " X: " << obstacleStackCompl[ind].x << " Y: " << obstacleStackCompl[ind].y << "\n";
    //     }
    // }
    // std::cout << "\n";

    // Publish data
    ROS_INFO("Found obstacles: %i", (int)obsPoses.obstacles.size());
    obstaclePosesPub.publish(obsPoses);

    obsPoses.obstacles.erase(obsPoses.obstacles.begin(),obsPoses.obstacles.end());
    obstacleStackCompl.erase(obstacleStackCompl.begin(),obstacleStackCompl.end());
}

bool ObjectDetection::initLocalisation(std::vector<ObstacleInLidar> &postObstacleStack, int &indexOffsetToObjRowSave, float minAngle, float maxAngle,
                                       int minAmountOfPoints, float rangeMultMin, float objectScanOffset, int outputMode)
{
    if(!lidarScan)
    {
        return false;
    }
    if(outputMode >= 1)
        std::cout << "\n";

    float angleRange = maxAngle - minAngle;
    int numberScansCompl = ceil((lidarScan->angle_max - lidarScan->angle_min) / lidarScan->angle_increment);
    int minIndex = ceil((minAngle/180*M_PI - lidarScan->angle_min) / lidarScan->angle_increment);
    int maxIndex = floor((maxAngle/180*M_PI - lidarScan->angle_min) / lidarScan->angle_increment);
    int numberScansEval = maxIndex - minIndex;
    bool distances[numberScansEval] = { };
    int counter = 0;
    int tmp_ctr=1;
    int tmp_ind;
    float scanDistances[numberScansEval];
    int counterForSaving = 0;
    float meanDist = 0.0;
    float distRatio = 0.0;
    float angle_mid;
    bool foundSomethingRetval = false;

    // ------------------------- Save Lidar Points in the selected range -------------------------
    for(int l=minIndex; l<maxIndex; ++l)
    {
        if (l > numberScansCompl)
            tmp_ind = l - numberScansCompl - 1;
        else
            tmp_ind = l;

        if(lidarScan->ranges[tmp_ind] == INFINITY)
            scanDistances[counterForSaving] = 20.0;
        else
            scanDistances[counterForSaving] = lidarScan->ranges[tmp_ind];

        counterForSaving++;
    }
    // -------------------------- Search minimum distance in range ------------------------------
    float min_range = scanDistances[0];
    for( int k = 1; k < numberScansEval; ++k)
    {
        if ( scanDistances[k] < min_range)
            min_range = scanDistances[k];
    }
    if (outputMode == 1 || outputMode == 2)
        ROS_INFO("min:%f",min_range);

    // ------------------------- Search for objects ---------------------------------------------

    int loop_counter = 0;
    while ( loop_counter < numberScansEval)
    {
        bool nextIsOk = true;
        if (scanDistances[loop_counter] < rangeMultMin*min_range)
        {
            tmp_ctr=1;
            nextIsOk = true;
            while(nextIsOk)
            {
                if( abs(scanDistances[loop_counter] - scanDistances[loop_counter+tmp_ctr]) < objectScanOffset)
                    tmp_ctr++; // Punkt gehört dazu
                else
                    nextIsOk = false;
            }
            if(tmp_ctr >= minAmountOfPoints)
            {
                // save obstacle
                ObstacleInLidar newObstacle;
                newObstacle.index_min = loop_counter;
                newObstacle.index_max = loop_counter + tmp_ctr;
                newObstacle.angle_min = minAngle + angleRange*newObstacle.index_min/numberScansEval;
                newObstacle.angle_max = minAngle + angleRange*newObstacle.index_max/numberScansEval;
                newObstacle.angle_mean = (newObstacle.angle_min + newObstacle.angle_max)/2;
                angle_mid = (newObstacle.angle_max-newObstacle.angle_min)/2;
                for(int a=0;a<tmp_ctr;++a)
                {
                    meanDist = meanDist + scanDistances[loop_counter+a];
                }
                newObstacle.distance = meanDist/tmp_ctr;
                newObstacle.x = newObstacle.distance*cos(newObstacle.angle_mean*M_PI/180.0);
                newObstacle.y = newObstacle.distance*sin(newObstacle.angle_mean*M_PI/180.0);
                obstacleStack.push_back(newObstacle);

                meanDist = 0.0;
            }
        }
        loop_counter = tmp_ctr + loop_counter;
        tmp_ctr = 1;
    }
    if (outputMode == 1 || outputMode == 2)
        ROS_INFO("Gefundene Objekte: %i",(int) obstacleStack.size());
    if (outputMode == 1 || outputMode == 2)
    {
        for(int ind = 0; ind < obstacleStack.size(); ++ind)
        {
            std::cout << " Objekt: " << ind << " Angle: " << obstacleStack[ind].angle_mean << " Distanz: " << obstacleStack[ind].distance
                      << " Index min: " << obstacleStack[ind].index_min << " X: " << obstacleStack[ind].x << " Y: " << obstacleStack[ind].y << "\n";
        }
    }

    // -------------------------------- calculate Distances between Object -------------------------

    for (std::vector<ObstacleInLidar>::iterator it = obstacleStack.begin() ; it != obstacleStack.end(); ++it)
    {
        for (std::vector<ObstacleInLidar>::iterator it2 = obstacleStack.begin() ; it2 != obstacleStack.end(); ++it2)
        {
            if( (*it).index_min == (*it2).index_min)
            {
                (*it).distancesToOtherObj.push_back(0.0);
            }
            else
            {
                float a = (*it).distance;
                float b = (*it2).distance;
                float mid1 = ( (*it).angle_min + (*it).angle_max)/2;
                float mid2 = ( (*it2).angle_min + (*it2).angle_max)/2;
                float gamma = abs(mid1-mid2);

                float distToObj = sqrt( pow(a,2) + pow(b,2) - 2*a*b*cos(gamma*M_PI/180.0) );
                (*it).distancesToOtherObj.push_back(distToObj);
            }
        }
    }

    if (outputMode == 2)
    {
        ROS_INFO("Auflistung Distanzen");
        for (std::vector<ObstacleInLidar>::iterator it = obstacleStack.begin() ; it != obstacleStack.end(); ++it)
        {
            for(int i = 0; i<(int) obstacleStack.size(); ++i)
            {
                std::cout << (*it).distancesToOtherObj[i] << " ";
            }
            std::cout << "\n";
        }
    }

    // --------------------------- Plausibilitätsprüfung -------------------------------------------------

    bool plausibilityOkay = false;
    int indexOffsetToObjRow;

    if((int) obstacleStack.size() == 3)
    {
        // prüfe ob es sich bei den 3 Objekten um Pfosten halten kann
        // Ein Abstand zwischen 2 Objekten muss ca. 33% des darauffolgenden Abstandes sein, allerdings nicht exakt
        // Die 3 Objekte müssen auf einer Linie liegen.
        indexOffsetToObjRow = 0;
        plausibilityOkay = plausibilityCheckForThreeObj(0, 1, 2, outputMode, indexOffsetToObjRow);
        if (plausibilityOkay)
        {
            indexOffsetToObjRowSave = indexOffsetToObjRow;
            postObstacleStack.push_back(obstacleStack[0]);
            postObstacleStack.push_back(obstacleStack[1]);
            postObstacleStack.push_back(obstacleStack[2]);

            foundSomethingRetval = true;
        }
    }
    else if((int) obstacleStack.size() > 3)
    {
        for(int i = 0; i < (int) obstacleStack.size()-2; ++i)
        {
            indexOffsetToObjRow = i;
            plausibilityOkay = plausibilityCheckForThreeObj(i, i+1, i+2, outputMode, indexOffsetToObjRow);
            if (plausibilityOkay)
            {
                indexOffsetToObjRowSave = indexOffsetToObjRow;
                postObstacleStack.push_back(obstacleStack[i]);
                postObstacleStack.push_back(obstacleStack[i+1]);
                postObstacleStack.push_back(obstacleStack[i+2]);

                foundSomethingRetval = true;
            }
        }
    }

    obstacleStack.erase(obstacleStack.begin(),obstacleStack.end());
    return foundSomethingRetval;
}
bool ObjectDetection::plausibilityCheckForThreeObj(int fstIndex, int sndIndex, int thrIndex, int outputMode, int offset)
{
    float distRatio = 0.0;
    bool retval = false;
    // Verhältnis
    if( obstacleStack[sndIndex].distancesToOtherObj[0+offset] < obstacleStack[sndIndex].distancesToOtherObj[2+offset])
        distRatio = obstacleStack[sndIndex].distancesToOtherObj[0+offset]/obstacleStack[sndIndex].distancesToOtherObj[2+offset] * 100.0;
    else
        distRatio = obstacleStack[sndIndex].distancesToOtherObj[2+offset]/obstacleStack[sndIndex].distancesToOtherObj[0+offset] * 100.0;
    if (outputMode == 1 || outputMode == 2)
        ROS_INFO("dist ratio: %f",distRatio);

    // Auf einer Geraden
    // Die Unterscheidung ist dazu da, damit die Gerade immer aus den 2 Punkten gebildet wird,
    // die am weitesten von einander entfernt sind.
    if( obstacleStack[sndIndex].distancesToOtherObj[0+offset] < obstacleStack[sndIndex].distancesToOtherObj[2+offset])
    {
        float m1 = (obstacleStack[thrIndex].y-obstacleStack[sndIndex].y) / (obstacleStack[thrIndex].x-obstacleStack[sndIndex].x);
        float t1 = obstacleStack[thrIndex].y - m1*obstacleStack[thrIndex].x;
        float m2 = -1.0/m1;
        float t2 = obstacleStack[fstIndex].y - m2*obstacleStack[fstIndex].x;
        float SP_x = (t2-t1)/(m1-m2);
        float SP_y = m2*SP_x + t2;

        // Test first object
        float cmpPunktGerade = sqrt(pow(SP_y-obstacleStack[fstIndex].y ,2)+pow(SP_x-obstacleStack[fstIndex].x ,2));
        if (outputMode == 1 || outputMode == 2)
            ROS_INFO("Geradenvergleich: %f", cmpPunktGerade);


        // Prüfe auf Plausiblität
        if ( distRatio < 42.0 && distRatio >27.0 && cmpPunktGerade < 0.1 && cmpPunktGerade > - 0.1 )
        {
            // scheint zu passen
            if (outputMode == 2)
                ROS_INFO("Found something.");
            retval = true;
        }
    }
    else
    {
        float m1 = (obstacleStack[sndIndex].y-obstacleStack[fstIndex].y) / (obstacleStack[sndIndex].x-obstacleStack[fstIndex].x);
        float t1 = obstacleStack[sndIndex].y - m1*obstacleStack[sndIndex].x;
        float m2 = -1.0/m1;
        float t2 = obstacleStack[thrIndex].y - m2*obstacleStack[thrIndex].x;
        float SP_x = (t2-t1)/(m1-m2);
        float SP_y = m2*SP_x + t2;

        // Test first object
        float cmpPunktGerade = sqrt(pow(SP_y-obstacleStack[thrIndex].y ,2)+pow(SP_x-obstacleStack[thrIndex].x ,2));
        if (outputMode == 1 || outputMode == 2)
            ROS_INFO("Geradenvergleich:%f", cmpPunktGerade);


        // Prüfe auf Plausiblität
        if ( distRatio < 45.0 && distRatio >25.0 && cmpPunktGerade < 0.3 && cmpPunktGerade > -0.3 )
        {
            // scheint zu passen
            if (outputMode == 2)
                ROS_INFO("Found something.");
            retval = true;
        }
    }
    return retval;
}

void ObjectDetection::searchMinDistanceInLidarRange(float &minRange, float minAngle, float maxAngle)
{
    if(!lidarScan)
    {
        return;
    }

    int numberScansCompl = ceil((lidarScan->angle_max - lidarScan->angle_min) / lidarScan->angle_increment);
    int minIndex = ceil((minAngle/180*M_PI - lidarScan->angle_min) / lidarScan->angle_increment);
    int maxIndex = floor((maxAngle/180*M_PI - lidarScan->angle_min) / lidarScan->angle_increment);
    int numberScansEval = maxIndex - minIndex;

    int tmp_ind;
    float scanDistances[numberScansEval];
    int counterForSaving = 0;

    // ------------------------- Save Lidar Points in the selected range -------------------------
    for(int l=minIndex; l<maxIndex; ++l)
    {
        if (l > numberScansCompl)
            tmp_ind = l - numberScansCompl - 1;
        else
            tmp_ind = l;

        if(lidarScan->ranges[tmp_ind] == INFINITY)
            scanDistances[counterForSaving] = 20.0;
        else
            scanDistances[counterForSaving] = lidarScan->ranges[tmp_ind];

        counterForSaving++;
    }

    // -------------------------- Search minimum distance in range ------------------------------
    float min_range = scanDistances[0];
    for( int k = 1; k < numberScansEval; ++k)
    {
        if ( scanDistances[k] < min_range)
            min_range = scanDistances[k];
    }
//    ROS_INFO("Min:%f",min_range);
    minRange = min_range;
}
void ObjectDetection::calculateMeanVarMedian(float &mean, float &var, float &median, std::vector<float> values, bool sorted)
{
    int size = (int) values.size();

    // calc mean
    for ( int i=0; i < size; ++i)
       mean = mean + values[i];
    mean = mean / size;

    // calc variance
    for(int i=0; i < size; ++i)
        var = var + pow(values[i] - mean,2);
    var = var / size;

    // calc median
    if (!sorted)
        sort(values.begin(), values.end());

    // calc median
    if ((values.size() % 2) == 0)
        median = values[size/2-1];
    else
        median = values[(size+1)/2-1];
}
void ObjectDetection::calculateAandBFromStartPos()
{
    if(receivecFirstLidar) // Wait till first Lidar data received!
    {
        if (measurements == 0)
         ROS_INFO("MODE 1: START DETECTING ON RIGHT AND LEFT SIDE!");
        // do 10 measurements on both sides
        if(measurements < 11)   // do 10 measurements
        {
            bool foundSomethingRight = initLocalisation(rightObstacleStack, indexOffsetToObjRow1, 30.0, 130.0, 1, rangeMultMin, objectScanOffset, outputMode); // Pfosten rechts vom Turtlebot
            if(foundSomethingRight)
            {
                if (!rightObstacleStack.empty())
                {
                    rightStackA.push_back(rightObstacleStack[0].distancesToOtherObj[2+indexOffsetToObjRow1]);
                }
            }

            bool foundSomethingLeft = initLocalisation(leftObstacleStack, indexOffsetToObjRow2, 230.0, 330.0, 1,rangeMultMin, objectScanOffset, outputMode); // Pfosten links vom Turtlebot
            if(foundSomethingLeft)
            {
                if (!leftObstacleStack.empty())
                    leftStackA.push_back(leftObstacleStack[0].distancesToOtherObj[2+indexOffsetToObjRow2]);
            }
            // calculate b and own location
            if(foundSomethingLeft && foundSomethingRight)
            {
                // indices for the 'lower' post
                if( leftObstacleStack[0].angle_mean > leftObstacleStack[2].angle_mean)
                    leftIndexForB = 0;
                else
                    leftIndexForB = 2;
                if( rightObstacleStack[0].angle_mean < rightObstacleStack[2].angle_mean)
                    rightIndexForB = 0;
                else
                    rightIndexForB = 2;

                 // Kosinussatz
                 // ROS_INFO("B: Distances: %f, %f: %f,%f",rightObstacleStack[rightIndexForB].distance, leftObstacleStack[leftIndexForB].distance, rightObstacleStack[rightIndexForB].angle_mean, leftObstacleStack[leftIndexForB].angle_mean);
                float gamma;
                float b_tmp1 = pow(rightObstacleStack[rightIndexForB].distance,2) + pow(leftObstacleStack[leftIndexForB].distance,2);
                float b_tmp2 = 2.0*rightObstacleStack[rightIndexForB].distance*leftObstacleStack[leftIndexForB].distance;
                if( cos(rightObstacleStack[rightIndexForB].angle_mean < leftObstacleStack[leftIndexForB].angle_mean) )
                    gamma =  (360.0 - leftObstacleStack[leftIndexForB].angle_mean) + rightObstacleStack[rightIndexForB].angle_mean;
                else
                    gamma = (rightObstacleStack[rightIndexForB].angle_mean-leftObstacleStack[leftIndexForB].angle_mean);
                float b_tmp3 = cos(abs(gamma*M_PI/180.0));
                float b = sqrt( b_tmp1 - b_tmp2*b_tmp3 );
                // ROS_INFO("CALC: %f, %f, %f, %f", b_tmp1, b_tmp2, gamma, b_tmp3);
                stackB.push_back(b);

                // own Position
                // geometry_msgs::Pose pose;
                // geometry_msgs::Point position;
                float locAngle = rightObstacleStack[rightIndexForB].angle_mean;
                float x = cos(locAngle*M_PI/180.0) * rightObstacleStack[rightIndexForB].distance;
                float y = b - sin(locAngle*M_PI/180.0) * rightObstacleStack[rightIndexForB].distance;
                float thetaZ;
                thetaZ = leftObstacleStack[0].y - leftObstacleStack[2].y;
                thetaZ = M_PI - atan(thetaZ*M_PI/180.0);

                // float thetaZ = abs(360.0 - leftObstacleStack[leftIndexForB].angle_mean);
                // thetaZ = thetaZ - rightObstacleStack[rightIndexForB].angle_mean;
                // thetaZ = 180.0 - thetaZ;
                geometry_msgs::Quaternion q = toQuaternion(0.0, 0.0, thetaZ);
                // pose.position = position;
                // pose.orientation = q;
                // ownPosition.push_back(pose);

                posX.push_back(x);
                posY.push_back(y);
                quatZ.push_back(q.z);
                quatW.push_back(q.w);
                posTheta.push_back(thetaZ);
                // ROS_INFO("THETA: %f",thetaZ);
                // ROS_INFO("HEHE: %f, %f, %f, %f", q.x, q.y, q.z, q.w);
            }
            measurements++;
        }
        else
        {
            measurementsDone = true;
            measurements = 0;
        }

        if (measurementsDone)
        {
            ROS_INFO("Measurements done");

            if (!rightStackA.empty())
            {
                sort(rightStackA.begin(), rightStackA.end());
                 // Calculate Mean, Variance and median
//                 for(int i=0; i < (int)rightStackA.size(); ++i)
//                     std::cout << rightStackA[i] << " ";
//                 std::cout << "\n";
                calculateMeanVarMedian(meanRA, varRA, medianRA, rightStackA, true );
                ROS_INFO("Right  Mean: %f, Var: %f, Median: %f",meanRA, varRA, medianRA);
                foundSomethingOnRightSide = true;
            }
            else
                ROS_INFO("Nothing found on right side!");

            if (!leftStackA.empty())
            {
                sort(leftStackA.begin(), leftStackA.end());
//                 for(int i=0; i < (int)leftStackA.size(); ++i)
//                     std::cout << leftStackA[i] << " ";
//                 std::cout << "\n";
                calculateMeanVarMedian(meanLA, varLA, medianLA, leftStackA, true );
                ROS_INFO("Left  Mean: %f, Var: %f, Median: %f",meanLA, varLA, medianLA);
                foundSomethingOnLeftSide = true;
            }
            else
                ROS_INFO("Nothing found on left side!");

            if (!stackB.empty())
            {
                sort(stackB.begin(), stackB.end());
//                 for(int i=0; i < (int)stackB.size(); ++i)
//                     std::cout << stackB[i] << " ";
//                 std::cout << "\n";
                calculateMeanVarMedian(meanB, varB, medianB, stackB, true );
                ROS_INFO("B  Mean: %f, Var: %f, Median: %f",meanB, varB, medianB);
            }

            if( foundSomethingOnRightSide && foundSomethingOnLeftSide)
            {
                initABOkay.data = true;
                // calculate a/b
                float a_value_forCalc;
                if(varLA > varRA)
                {
                    if(varRA < 0.001 && varB < 0.001)
                    {
                        // Save a/b Ratio
                        abRatio.data = medianRA / medianB;
                        a_value_forCalc = medianRA;
                    }
                    else
                        abRatio.data = 0.0;
                }
                else
                {
                    if(varLA < 0.001 && varB < 0.001)
                    {
                        // Save a/b Ratio
                        abRatio.data = medianLA / medianB;
                        a_value_forCalc = medianLA;
                    }
                    else
                        abRatio.data = 0.0;
                }

                // calculate own position
                if(varRA < 0.001 && varLA < 0.001 && varB < 0.001)
                {
                    if(!posX.empty() && !posY.empty() && !quatZ.empty() && !quatW.empty() && !posTheta.empty())
                    {
                        sort(posX.begin(), posX.end());
                        sort(posY.begin(), posY.end());
                        sort(quatZ.begin(), quatZ.end());
                        sort(quatW.begin(), quatW.end());
                        sort(posTheta.begin(), posTheta.end());

                        int size = posX.size();

                        if ((size % 2) == 0)
                        {
                            ownPosition.position.x = posX[size/2-1];
                            ownPosition.position.y = posY[size/2-1];
                            ownPosition.orientation.z = quatZ[size/2-1];
                            ownPosition.orientation.w = quatW[size/2-1];

                            ownPosition2D.x = posX[size/2-1];
                            ownPosition2D.y = posY[size/2-1];
                            ownPosition2D.theta = posTheta[size/2-1];
                        }
                        else
                        {
                            ownPosition.position.x = posX[(size+1)/2-1];
                            ownPosition.position.y = posY[(size+1)/2-1];
                            ownPosition.orientation.z = quatZ[(size+1)/2-1];
                            ownPosition.orientation.w = quatW[(size+1)/2-1];

                            ownPosition2D.x = posX[(size+1)/2-1];
                            ownPosition2D.y = posY[(size+1)/2-1];
                            ownPosition2D.theta = posTheta[(size+1)/2-1];
                        }
                    }
                    ownPosition.position.z = 0.0;
                    ownPosition.orientation.x = 0.0;
                    ownPosition.orientation.y = 0.0;
                }

                // calculate edge points of the goals
                if(varRA < 0.001 && varLA < 0.001 && varB < 0.001)
                {
                    enemyGoal.center.x = 2.625 * a_value_forCalc;
                    enemyGoal.center.y = medianB/2;
                    enemyGoal.links_unten.x = 2.5 * a_value_forCalc;
                    enemyGoal.links_unten.y = medianB/3;
                    enemyGoal.rechts_unten.x = 2.5 * a_value_forCalc;
                    enemyGoal.rechts_unten.y = 2*medianB/3;
                    enemyGoal.links_oben.x = 2.75 * a_value_forCalc;
                    enemyGoal.links_oben.y = medianB/3;
                    enemyGoal.rechts_oben.x = 2.75 * a_value_forCalc;
                    enemyGoal.rechts_oben.y = 2*medianB/3;
                }
                // calculate game field corners
                if(varRA < 0.001 && varLA < 0.001 && varB < 0.001)
                {
                    gameFieldCorners.rechts_unten_GF.x = 0.0;
                    gameFieldCorners.rechts_unten_GF.y = medianB;
                    gameFieldCorners.links_unten_GF.x = 0.0;
                    gameFieldCorners.links_unten_GF.y = 0.0;
                    gameFieldCorners.rechts_oben_GF.x = a_value_forCalc * 3;
                    gameFieldCorners.rechts_oben_GF.y = medianB;
                    gameFieldCorners.links_oben_GF.x = a_value_forCalc * 3;
                    gameFieldCorners.links_oben_GF.y = 0.0;

                    gameFieldCorners.rechts_unten_RF.x = gameFieldCorners.rechts_unten_GF.x - ownPosition2D.x;
                    gameFieldCorners.rechts_unten_RF.y = gameFieldCorners.rechts_unten_GF.y - ownPosition2D.y;
                    gameFieldCorners.links_unten_RF.x = gameFieldCorners.links_unten_GF.x - ownPosition2D.x;
                    gameFieldCorners.links_unten_RF.y = gameFieldCorners.links_unten_GF.y - ownPosition2D.y;
                    gameFieldCorners.rechts_oben_RF.x = gameFieldCorners.rechts_oben_GF.x - ownPosition2D.x;
                    gameFieldCorners.rechts_oben_RF.y = gameFieldCorners.rechts_oben_GF.y - ownPosition2D.y;
                    gameFieldCorners.links_oben_RF.x = gameFieldCorners.links_oben_GF.x - ownPosition2D.x;
                    gameFieldCorners.links_oben_RF.y = gameFieldCorners.links_oben_GF.y - ownPosition2D.y;
                }

                // Publish Data
    //            abRatio.data = 0.0;
                abRatioPub.publish(abRatio);
                ownPositionPub.publish(ownPosition);
                ownPosition2DPub.publish(ownPosition2D);
                enemyGoalPub.publish(enemyGoal);
                gameFieldPub.publish(gameFieldCorners);


                // Reset Variables
                measurementsDone = false;
                rightStackA.erase(rightStackA.begin(), rightStackA.end());
                leftStackA.erase(leftStackA.begin(), leftStackA.end());
                stackB.erase(stackB.begin(), stackB.end());
                posX.erase(posX.begin(), posX.end());
                posY.erase(posY.begin(), posY.end());
                quatZ.erase(quatZ.begin(), quatZ.end());
                quatW.erase(quatW.begin(), quatW.end());
            }
            else
                initABOkay.data = false;

            initABSuccessPub.publish(initABOkay);
            // ownPosition.erase(ownPosition.begin(), ownPosition.end());

            runtime = ros::Time::now().toSec();
            diff = runtime - begin;
            ROS_INFO("Durchlaufzeit: %f", diff);
            begin = runtime;

            std::cout << "\n";
         }
     }

     // --------------- Erase erverything --------------------
     rightObstacleStack.erase(rightObstacleStack.begin(), rightObstacleStack.end());
     leftObstacleStack.erase(leftObstacleStack.begin(), leftObstacleStack.end());

     meanRA = 0;
     varRA = 0;
     medianRA = 0;
     meanLA = 0;
     varLA = 0;
     medianLA = 0;
     meanB = 0;
     varB = 0;
     medianB = 0;

     foundSomethingOnRightSide = false;
     foundSomethingOnLeftSide = false;
}
void ObjectDetection::update()
{
    if(receivecFirstMode)
    {
        // Do processing depending on current data in Topic
        if(currentCalcMode->data == 1)
        {
            calculateAandBFromStartPos();
        }
        else if (currentCalcMode->data == 2)
        {
            scanForObjects(2.0);
        }
    }

    isAlive.data = true;
    alivePub.publish(isAlive);
}

void ObjectDetection::startObjectDetection()
{
    ros::Rate rate(10);
    ROS_INFO("Object Detection Node started!");

    // Keep spinning loop until user presses Ctrl+C
    while (ros::ok())
    {
        update();
        ros::spinOnce(); // Need to call this function often to allow ROS to process incoming messages
        rate.sleep();
    }
}



























//         if(receivecFirstLidar) // Wait till first Lidar data received!
//         {
//             if (measurements == 0)
//                 ROS_INFO("START DETECTING ON RIGHT AND LEFT SIDE!");
//             // do 10 measurements on both sides
//             if(measurements < 11)   // do 10 measurements
//             {
//                 bool foundSomethingRight = initLocalisation(rightObstacleStack, indexOffsetToObjRow1, 20.0, 150.0, 1,rangeMultMin, objectScanOffset, outputMode); // Pfosten rechts vom Turtlebot
//                 if(foundSomethingRight)
//                 {
//                     if (!rightObstacleStack.empty())
//                     {
//                         rightStackA.push_back(rightObstacleStack[0].distancesToOtherObj[2+indexOffsetToObjRow1]);
//                     }
//                 }
//
//                 bool foundSomethingLeft = initLocalisation(leftObstacleStack, indexOffsetToObjRow2, 230.0, 320.0, 1,rangeMultMin, objectScanOffset, outputMode); // Pfosten links vom Turtlebot
//                 if(foundSomethingLeft)
//                 {
//                     if (!leftObstacleStack.empty())
//                         leftStackA.push_back(leftObstacleStack[0].distancesToOtherObj[2+indexOffsetToObjRow2]);
//                 }
//                 // calculate b
//                 if(foundSomethingLeft && foundSomethingRight)
//                 {
//                     if( leftObstacleStack[0].angle_mean > leftObstacleStack[2].angle_mean)
//                         leftIndexForB = 0;
//                     else
//                         leftIndexForB = 2;
//                     if( rightObstacleStack[0].angle_mean < rightObstacleStack[2].angle_mean)
//                         rightIndexForB = 0;
//                     else
//                         rightIndexForB = 2;
//
//                     // Kosinussatz
//                     // ROS_INFO("B: Distances: %f, %f: %f,%f",rightObstacleStack[rightIndexForB].distance, leftObstacleStack[leftIndexForB].distance, rightObstacleStack[rightIndexForB].angle_mean, leftObstacleStack[leftIndexForB].angle_mean);
//                     float gamma;
//                     float b_tmp1 = pow(rightObstacleStack[rightIndexForB].distance,2) + pow(leftObstacleStack[leftIndexForB].distance,2);
//                     float b_tmp2 = 2.0*rightObstacleStack[rightIndexForB].distance*leftObstacleStack[leftIndexForB].distance;
//                     if( cos(rightObstacleStack[rightIndexForB].angle_mean < leftObstacleStack[leftIndexForB].angle_mean) )
//                         gamma =  (360.0 - leftObstacleStack[leftIndexForB].angle_mean) + rightObstacleStack[rightIndexForB].angle_mean;
//                     else
//                         gamma = (rightObstacleStack[rightIndexForB].angle_mean-leftObstacleStack[leftIndexForB].angle_mean);
//                     float b_tmp3 = cos(abs(gamma*M_PI/180.0));
//                     float b = sqrt( b_tmp1 - b_tmp2*b_tmp3 );
//                     // ROS_INFO("CALC: %f, %f, %f, %f", b_tmp1, b_tmp2, gamma, b_tmp3);
//                     stackB.push_back(b);
//                 }
//                 measurements++;
//             }
//             else
//             {
//                 measurementsDone = true;
//                 measurements = 0;
//             }
//
//             if (measurementsDone)
//             {
//                 ROS_INFO("Measurements done");
//
//                 if (!rightStackA.empty())
//                 {
//                     sort(rightStackA.begin(), rightStackA.end());
//                     // Calculate Mean, Variance and median
//                     for(int i=0; i < (int)rightStackA.size(); ++i)
//                         std::cout << rightStackA[i] << " ";
//                     std::cout << "\n";
//                     calculateMeanVarMedian(meanRA, varRA, medianRA, rightStackA, true );
//                     ROS_INFO("Right  Mean: %f, Var: %f, Median: %f",meanRA, varRA, medianRA);
//                 }
//                 else
//                     ROS_INFO("Nothing found on right side!");
//
//                 if (!leftStackA.empty())
//                 {
//                     sort(leftStackA.begin(), leftStackA.end());
//                     for(int i=0; i < (int)leftStackA.size(); ++i)
//                         std::cout << leftStackA[i] << " ";
//                     std::cout << "\n";
//                     calculateMeanVarMedian(meanLA, varLA, medianLA, leftStackA, true );
//                     ROS_INFO("Left  Mean: %f, Var: %f, Median: %f",meanLA, varLA, medianLA);
//                 }
//                 else
//                     ROS_INFO("Nothing found on left side!");
//
//                 if (!stackB.empty())
//                 {
//                     sort(stackB.begin(), stackB.end());
//                     for(int i=0; i < (int)stackB.size(); ++i)
//                         std::cout << stackB[i] << " ";
//                     std::cout << "\n";
//                     calculateMeanVarMedian(meanB, varB, medianB, stackB, true );
//                     ROS_INFO("B  Mean: %f, Var: %f, Median: %f",meanB, varB, medianB);
//
//                 }
//
//                 measurementsDone = false;
//                 rightStackA.erase(rightStackA.begin(), rightStackA.end());
//                 leftStackA.erase(leftStackA.begin(), leftStackA.end());
//
// //                for(int i=0; i < (int) stackB.size();++i)
// //                    std::cout << stackB[i] << " ";
// //                std::cout < "\n";
//
//                 runtime = ros::Time::now().toSec();
//                 diff = runtime - begin;
//                 ROS_INFO("Durchlaufzeit: %f", diff);
//                 begin = runtime;
//
//                 std::cout << "\n";
//             }
//         }
//
//         // --------------- Erase erverything --------------------
//         rightObstacleStack.erase(rightObstacleStack.begin(), rightObstacleStack.end());
//         leftObstacleStack.erase(leftObstacleStack.begin(), leftObstacleStack.end());
//
//         meanRA = 0;
//         varRA = 0;
//         medianRA = 0;
//         meanLA = 0;
//         varLA = 0;
//         medianLA = 0;
//         meanB = 0;
//         varB = 0;
//         medianB = 0;

/*
 * //        std::vector<ObstacleInLidar>::iterator itr;
//        std::vector<ObstacleInLidar>::iterator itr2;
//        // Iterator auf mittleres Objekt
//        itr = obstacleStack.begin();
//        itr++;
//        // Verhältnis
//        if( (*itr).distancesToOtherObj[0] < (*itr).distancesToOtherObj[2])
//            distRatio = (*itr).distancesToOtherObj[0]/(*itr).distancesToOtherObj[2] * 100.0;
//        else
//            distRatio = (*itr).distancesToOtherObj[2]/(*itr).distancesToOtherObj[0] * 100.0;
//        ROS_INFO("dist ratio: %f",distRatio);

//        // Auf einer Geraden
//        itr = obstacleStack.begin();
//        itr2 = obstacleStack.begin();
//        itr2++;
//        float m = ((*itr2).y-(*itr).y) / ((*itr2).x-(*itr).x);
//        float t = (*itr).y - m*(*itr).x;
//        // Test third object
//        itr2++;
//        float cmpPunktGerade = (*itr2).y - (m*(*itr2).x + t);
//        ROS_INFO("Geradenvergleich:%f", cmpPunktGerade);

//        // Prüfe auf Plausiblität
//        if ( distRatio < 45.0 && distRatio >25.0 && cmpPunktGerade < 0.5 && cmpPunktGerade > -0.5 )
//        {
//            // scheint zu passen
//            ROS_INFO("PASST! Abstand a: %f", (*itr).distancesToOtherObj[2]);
//        }
*/
