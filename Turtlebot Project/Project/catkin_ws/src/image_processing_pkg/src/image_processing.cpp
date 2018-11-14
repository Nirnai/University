#include <geometry_msgs/Twist.h>

#include "image_processing.h"

ImageProcessing::ImageProcessing() {

    yellowTeam = false;
    blueTeam   = false;
    receivedFirstCloud = false;
    initTCOkay.data = false;
    receivecFirstObstacles = false;
    puckAtMidFound.data = false;

    bestPuck.x = 0.0;
    bestPuck.y = 0.0;
    bestPuck.theta = 0.0;
    isAlive.data = false;

    testWithPictureMode     = false;
    receivecFirstMode       = false;
    correctTeamColorDefined = false;
    testPictureFilename = "/home/ga84wih/Pictures/images/TurtleImg11.jpg";

    teamColorPub   = node.advertise<std_msgs::String>  ("/image_processing/teamColor",10);
    drivingGoalPub = node.advertise<geometry_msgs::Pose2D>("/image_processing/puckGoal",10);
    alivePub = node.advertise<std_msgs::Bool>("/image_processing/isAlive",10);
    initTeamColorSuccessPub = node.advertise<std_msgs::Bool>("/image_processing/initTCSuccess",10);
    puckAtMidPub = node.advertise<std_msgs::Bool>("/image_processing/puckAtMid",10);

    modiSub    = node.subscribe("/image_processing/mode", 1, &ImageProcessing::modeCallback, this);
    img        = node.subscribe("/camera/rgb/image_rect_color",10, &ImageProcessing::scanImageCallback, this);
    sub_points = node.subscribe("/camera/depth/points", 1, &ImageProcessing::cloudCallback, this);
    correctTeamColorSub = node.subscribe("/communication/correctTeamColor", 1, &ImageProcessing::correctTeamColorCallback, this);
    obstaclesSub = node.subscribe("/object_detection/obstacle_poses", 1, &ImageProcessing::obsCallback, this);


    puckFeatures.push_back(0.35);
    puckFeatures.push_back(0.09);
    puckFeatures.push_back(0.02);
    puckFeatures.push_back(0.01);
    puckFeatures.push_back(0.0);
    puckFeatures.push_back(0.0);
    puckFeatures.push_back(0.0);

}
void ImageProcessing::scanImageCallback( const sensor_msgs::Image::ConstPtr  & img)
{
    image  = img;
    cv_ptr = cv_bridge::toCvCopy(img,sensor_msgs::image_encodings::BGR8);
    // ROS_INFO("Bild erhalten!!!!!11111");
}

void ImageProcessing::correctTeamColorCallback(const std_msgs::String &str)
{
    correctTeamColor.data   = str.data;
    correctTeamColorDefined = true;
}

void ImageProcessing::cloudCallback(const sensor_msgs::PointCloud2ConstPtr& cloud_msg)
{
    // ROS_INFO("hi");
    global_cloud_msg = cloud_msg;
    //double stamp_3 = cloud_msg->header.stamp.toSec();
    //ROS_INFO("PointCloud Callback gets called\n TIME_STAMP_PointCloud: %f\n", stamp_3);
    if( (cloud_msg->width * cloud_msg->height) == 0 )
    {
       // not a dense cloud; return
       return;
    }
    try
    {

       flag_image_ok_point = true;
       //        geometry_msgs::Point p1;
       //        p1 = pixelTo3DPoint(cloud_msg, 240, 320);
       //        std::cout<<"x="<<p1.x<<", y="<<p1.y<<", z="<<p1.z<<std::endl;
    }
    catch( std::runtime_error e) {
       ROS_ERROR_STREAM("Error in converting point cloud to image: "<<e.what());
    }
    //       safety_ton_points = ros::Time::now();
    receivedFirstCloud = true;
}
void ImageProcessing::modeCallback(const std_msgs::ByteConstPtr &mode)
{
    currentCalcMode = mode;
    if(!receivecFirstMode)
        receivecFirstMode = true;

}
void ImageProcessing::obsCallback(const image_processing_pkg::obstacle_poses &obs)
{
    obstaclesArround = obs;
    if(!receivecFirstObstacles)
        receivecFirstObstacles = true;
}

bool ImageProcessing::calculateTeamColor()
{
    if (!testWithPictureMode)
    {
        if (!cv_ptr)
        {
            return true;
        }
    }
    saveIncomingImage();
    HSVColorRange yellow;
    // Thresholds for yellow
    yellow.iLowH = 21;
    yellow.iHighH = 39;
    yellow.iLowS = 40;
    yellow.iHighS = 255;
    yellow.iLowV = 40;
    yellow.iHighV = 255;

    std::vector<std::vector<Point> >  contours;
    std::vector<std::vector<Point> >  approxPolygons;
    std::vector<std::vector<double> > huDescriptor;
    std::vector<Point2f>              mc;
    int cutImg = 3;
    bool showResults = false;
    shapeDetection(contours, approxPolygons, huDescriptor, mc, yellow, cutImg, showResults);
    teamColor.data = "blue";
    for (int i = 0; i< contours.size(); i++)
    {
        if ((contourArea(contours[i]) > 8000) /*&& ((int)approxPolygons[i].size() <= 6)*/)
        {
            ROS_INFO("Yellow Goal detected! %f, %i",contourArea(contours[i]), (int)approxPolygons[i].size());
            teamColor.data = "yellow";
        }
    }
    teamColorPub.publish(teamColor);
    initTCOkay.data = true;
    initTeamColorSuccessPub.publish(initTCOkay);

//    shapeMatching(contours, approxPolygons, huDescriptor, "yellow");
    return true;
}
void ImageProcessing::saveIncomingImage()
{
    if (!cv_ptr)
    {
        return;
    }
    imwrite("/home/ga84wih/Pictures/TurtleImg.jpg",cv_ptr->image);
    ROS_INFO("Image saved");
    // imshow("aufgenommenes Bild",cv_ptr->image);
    // waitKey(0);
}
int ImageProcessing::detectGreenLimitPost()
{
    HSVColorRange green;
    // Thresholds for green
    green.iLowH = 40;
    green.iHighH = 80;
    green.iLowS = 50;
    green.iHighS = 255;
    green.iLowV = 90;
    green.iHighV = 255;

    std::vector<std::vector<Point> > contours;
    std::vector<std::vector<Point> > approxPolygons;
    std::vector<std::vector<double> > huDescriptor;
    std::vector<Point2f>              mc;
    shapeDetection(contours, approxPolygons, huDescriptor, mc, green, 0, true);
    //shapeDetection(contours, approxPolygons, huDescriptor, green, 0, true);

    return -1;
}

bool ImageProcessing::shapeDetection(std::vector<std::vector<Point> > &contours, std::vector<std::vector<Point> > &approxPolygons,
                                     std::vector<std::vector<double> > &huDescriptor, std::vector<Point2f> &mc, HSVColorRange colour,
                                     int cutImage, bool showResults)
{
    // Input handle to witch between test mode with using images and using turtlebot
    Mat dg_img;
    if (testWithPictureMode)
    {
        dg_img = imread(testPictureFilename);
        if(! dg_img.data)
        {
            ROS_INFO("Could not open or find the image");
            return false;
        }
    }
    else
    {
        if (!cv_ptr)
            return false;
        else
            dg_img = cv_ptr->image;
    }



    // Blur image
    GaussianBlur(dg_img, dg_img, Size(13, 13), 0);

    // Define ROI and cut
    if (cutImage == 2)
        dg_img = cv::Mat(dg_img, cv::Rect(0,dg_img.rows/2, dg_img.cols, dg_img.rows/2)).clone();
    else if(cutImage == 3)
        dg_img = cv::Mat(dg_img, cv::Rect(0,1*dg_img.rows/4, dg_img.cols, 3*dg_img.rows/4)).clone();

    // --------------------------- Detect regions by colour ------------------------------------------

    //        mc[i] = Point2f(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
    Mat imgHSV;
    cvtColor(dg_img, imgHSV, COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV
    Mat imgThresholded;
    inRange(imgHSV, Scalar(colour.iLowH, colour.iLowS, colour.iLowV),
            Scalar(colour.iHighH, colour.iHighS, colour.iHighV), imgThresholded); //Threshold the image

    // erode and dilate to remove noise
    // morphological opening (remove small objects from the foreground)
    erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
    dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
    // morphological closing (fill small holes in the foreground)
    dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
    erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );

    copyMakeBorder(imgThresholded, imgThresholded, 10, 10, 10, 10, BORDER_CONSTANT, Scalar(0,0,0));

    if (showResults)
        imshow("threshold",imgThresholded);

    // --------------------------- get moments (contours) --------------------------------------------
    RNG rng(12345);
    Mat srcThresh;
    double otsu;
    otsu = threshold(imgThresholded, srcThresh, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
    Mat cannyOut;
    Canny(imgThresholded, cannyOut, otsu, otsu * 1 / 2, 3, 1);
    // Find contours
//    std::vector<std::vector<Point> > contours;
    std::vector<Vec4i> hierarchy;

    findContours(cannyOut, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

    if(showResults)
        ROS_INFO("anzahl contours: %i",(int) contours.size());

    // Get the moments
    std::vector<Moments> mu(contours.size());
    for (int i = 0; i < contours.size(); i++)
    {
        mu[i] = moments(contours[i], false);
    }

    //  Get the mass centers:
    //std::vector<Point2f> mc(contours.size());
    for (int i = 0; i < contours.size(); i++)
    {
        mc.push_back(Point2f(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00));
    }

    Mat drawing = Mat::zeros(cannyOut.size(), CV_8UC3);
    std::string sObjectNumber;          // string which will contain the result
    std::ostringstream sContourNumber;   // stream used for the conversion
    if(showResults)
    {
        // Draw contours
        for (int i = 0; i< contours.size(); i++)
        {
            sContourNumber << i;
            sObjectNumber = sContourNumber.str();   // Convert int to string
            Point pCoordinates(mc[i].x + 3, mc[i].y - 3);   // Text's coordinates (A little bit off from mass center)
            Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
            drawContours(drawing, contours, i, color, 2, 8, hierarchy, 0, Point());
            circle(drawing, mc[i], 4, color, -1, 8, 0);     // Draw mass center
            putText(drawing, sObjectNumber, pCoordinates, CV_FONT_HERSHEY_COMPLEX, 1, color, 2, 8); // Write object number
            sContourNumber.str("");     // Clear string
            sContourNumber.clear();     // Clear any error flags
        }
    }

    // --------------------------- Build Descripor with HuMoments --------------------------------------------------
    std::vector<double> huPoints;
    double hu[7];
    for (int i = 0; i < contours.size(); i++)
    {
        HuMoments(mu[i], hu);
        huPoints.insert(huPoints.begin(), hu, hu + 7);

        if(showResults)
        {
            std::cout << "Contour: " << i << " Area: " << contourArea(contours[i]) << " Length: "
                      << arcLength(contours[i], true) << " Mass Center: " << mc[i] << "\n";
            for (int j = 0; j < 7; j++)
            {
                std::cout << "Contour: " << i << " Hu: " << j << " Result: " << hu[j] << "\n";
            }
            std::cout << "\n";
        }

        huDescriptor.push_back(huPoints);
        huPoints.erase(huPoints.begin(), huPoints.end());
    }
    if(showResults)
        imshow("Contours", drawing);

    // ----------------------- Shape Detection ---------------------------------------------------
    Mat drawShape = Mat::zeros(cannyOut.size(), CV_8UC3);
    std::vector<Point> result;
    for( int k =0; k<contours.size(); ++k)
    {
        approxPolyDP(contours[k], result, arcLength(Mat(contours[k]),true)*0.02, true);
        approxPolygons.push_back(result);
        if(showResults)
            ROS_INFO("result size: %i",(int) result.size());
    }

    if(showResults)
    {
        Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
        for (int i = 0; i < approxPolygons.size(); ++i)
        {
            drawContours(drawShape, approxPolygons, i, color, 2, 8, hierarchy, 0, Point());
        }
        imshow("approxPolyDP", drawShape);
        waitKey(0);
    }
}

void ImageProcessing::shapeMatching(std::vector<std::vector<Point> > &contours, std::vector<std::vector<Point> > &approxPolygons,
                                    std::vector<std::vector<double> > &huDescriptor, const char* currentColour)
{
    if(! strcmp(currentColour, "green"))
    {

    }
    else if(! strcmp(currentColour, "yellow"))
    {

    }
}

void ImageProcessing::drainClassifier()
{
    // Für gelbes tor:
    HSVColorRange yellow;
    // Thresholds for green
    yellow.iLowH = 22;
    yellow.iHighH = 38;
    yellow.iLowS = 80;
    yellow.iHighS = 255;
    yellow.iLowV = 80;
    yellow.iHighV = 255;

    std::vector<std::vector<Point> > contours;
    std::vector<std::vector<Point> > approxPolygons;
    std::vector<std::vector<double> > huDescriptor;
    std::vector<Point2f>              mc;
    shapeDetection(contours, approxPolygons, huDescriptor, mc, yellow, 3, true);
    //shapeDetection(contours, approxPolygons, huDescriptor, yellow, 3, true);

    shapeMatching(contours, approxPolygons, huDescriptor, "yellow");
}

geometry_msgs::Point ImageProcessing::pixelTo3DPoint(int row_num,  int col_num){
    geometry_msgs::Point p;
    float X=0.0;
    float Y=0.0;
    float Z=0.0;
    int arrayPosition = row_num  * 5120 + col_num * 16 ;
    int arrayPosX = arrayPosition + global_cloud_msg->fields[0].offset;
    int arrayPosY = arrayPosition + global_cloud_msg->fields[1].offset;
    int arrayPosZ = arrayPosition + global_cloud_msg->fields[2].offset;
    memcpy(&X, &global_cloud_msg->data[arrayPosX], sizeof(float));
    memcpy(&Y, &global_cloud_msg->data[arrayPosY], sizeof(float));
    memcpy(&Z, &global_cloud_msg->data[arrayPosZ], sizeof(float));
    p.x = X;
    p.y = Y;
    p.z = Z;
    return p;
}
bool ImageProcessing::findNearestObstacle(geometry_msgs::Pose2D &bestMatch, float x, float y, float theta)
{
    float errorX;
    float errorY;
    float error;
    float bestError = 10000.0;
    int bestIndex;

    if(!receivecFirstObstacles)
        return false;

    for( int i = 0; i < obstaclesArround.obstacles.size(); ++i)
    {
        if( obstaclesArround.obstacles[i].theta > 100 && obstaclesArround.obstacles[i].theta < 260)
        {
            // calculate x,y error
            errorX = abs( (x-0.15) - obstaclesArround.obstacles[i].x);
            errorY = abs( (y+0.03) - obstaclesArround.obstacles[i].y);
            error = errorX + errorY;
            if ( error < bestError)
            {
                bestIndex = i;
                bestError = error;
            }
        }
    }

    //
    bestMatch = obstaclesArround.obstacles[bestIndex];
    if(bestError < 1.0)
        return true;
    else
        return false;
}

void ImageProcessing::update()
{
    HSVColorRange color;
    std::vector<double> tmp;
    geometry_msgs::Pose2D puck;
    geometry_msgs::Point  puck3D;
    float dist;
    float bestDist = 10000.0;

    if (!correctTeamColorDefined){
      ROS_INFO("No correctTeamColor recieved, update return");
      return;
    }

    if (correctTeamColor.data == "yellow")
    {
        // yellow
        color.iLowH = 18;
        color.iHighH = 40;
        color.iLowS = 40;
        color.iHighS = 255;
        color.iLowV = 40;
        color.iHighV = 255;
    }
    else if (correctTeamColor.data == "blue")
    {
        // blue
        color.iLowH = 75;
        color.iHighH = 130;
        color.iLowS = 0;
        color.iHighS = 255;
        color.iLowV = 0;
        color.iHighV = 255;
    }

    std::vector<std::vector<Point> >  contours;
    std::vector<std::vector<Point> >  approxPolygons;
    std::vector<std::vector<double> > huDescriptor;
    std::vector<Point2f>              mc;

    int   cutImg      = 0;
    bool  showResults = false;
    float puckError;
    float distMidToX;

    shapeDetection(contours, approxPolygons, huDescriptor, mc, color, cutImg, showResults);

    //shapeDetection(contours, approxPolygons, huDescriptor, color, cutImg, showResults);

    // for (int i = 0; i< contours.size(); i++)
    // {
    //     if ((contourArea(contours[i]) > 10000) && ((int)approxPolygons[i].size() == 4))
    //     {
    //         ROS_INFO("Yellow Goal detected! %f, %i",contourArea(contours[i]), (int)approxPolygons[i].size());
    //         teamColor.data = "yellow";
    //     }
    // }
    // Option 1
    for( int i = 0; i < contours.size(); ++i)
    {
        if ((contourArea(contours[i]) > 500) )
        {
            puckError = 0.0;
            for( int j = 0; j < (int) tmp.size(); ++j)
                puckError = puckError + abs(tmp[j] - puckFeatures[j]);

            if ( puckError < 0.3)
            {
//                tempMC = mc[i];
//                ROS_INFO("contour %i: MC: %f, %f",i, mc[i].x, mc[i].y);
                distMidToX = mc[i].x-10.0 - 320.0;
//                ROS_INFO("dist: %f", distMidToX);
                if( abs(distMidToX) < 10)
                {
                    // Publish stop
                    puckAtMidFound.data = true;
                    if(receivecFirstObstacles)
                    {
                        for( int i = 0; i < obstaclesArround.obstacles.size(); ++i)
                        {
                            if( obstaclesArround.obstacles[i].theta < 190.0 && obstaclesArround.obstacles[i].theta > 170.0)
                            {
                                // könnte der Puck sein

                                puck.x = obstaclesArround.obstacles[i].x;
                                puck.y = obstaclesArround.obstacles[i].y;
                                dist = sqrt(pow(puck.x,2) + pow(puck.y,2));
                                if( dist < bestDist)
                                {
                                    bestPuck = puck;
                                    bestDist = dist;
                                }
                            }
                        }
                    }
                }
            }
        }
    }





    // Option 2
    //
    // for( int i = 0; i < contours.size(); ++i)
    // {
    //     tmp = huDescriptor[i];
    //
    //     if ((contourArea(contours[i]) > 500) )
    //     {
    //         for( int j = 0; j < (int) tmp.size(); ++j)
    //             puckError = puckError + abs(tmp[j] - puckFeatures[j]);
    //
    //         ROS_INFO("Puck Error: %f", puckError);
    //
    //         if(receivedFirstCloud)
    //         {
    //             if (true )//puckError < 2.0)
    //             {    //TODO Test Threshold for puckError here
    //                 puck3D = pixelTo3DPoint(mc[i].x, mc[i].y);
    //                 puck.x     = puck3D.z;
    //                 puck.y     = -puck3D.x;
    //                 puck.theta = atan2(puck.x, puck.y);
    //                 dist = sqrt(pow(puck.x,2) * pow(puck.y,2));
    //                 if( dist < bestDist)
    //                 {
    //                     bestPuck = puck;
    //                     bestDist = dist;
    //                 }
    //             }
    //         }
    //     }
    // }

    // Search best match
    // geometry_msgs::Pose2D puckToPublish;
    // bool retvalFNO;
    // retvalFNO = findNearestObstacle(puckToPublish, bestPuck.x, bestPuck.y, bestPuck.theta);
    // if(retvalFNO)
    //     drivingGoalPub.publish(puckToPublish);
    drivingGoalPub.publish(bestPuck);
    puckAtMidPub.publish(puckAtMidFound);
}

void ImageProcessing::startBV()
{
    ros::Rate rate(1);
    ROS_INFO("Start BV");

    if (testWithPictureMode)
        ROS_INFO("Ted mode: on, test images will be used");
    // Keep spinning loop until user presses Ctrl+C
    while (ros::ok())
    {
//        if (global_cloud_msg)position
//        {
//            geometry_msgs::Point p = pixelTo3DPoint(global_cloud_msg, 520, 240 );
//            ROS_INFO("X: %f, Y: %f, Z: %f", p.x, p.y, p.z);
//        }
        if(receivecFirstMode)
        {
            if ( currentCalcMode->data == 1)
            {
                ROS_INFO("Checking Team Color");
                calculateTeamColor();
            }
            if ( currentCalcMode->data == 2)
            {
                update();
            }
        }

        isAlive.data = true;
        alivePub.publish(isAlive);
        ros::spinOnce(); // Need to call this function often to allow ROS to process incoming messages
        rate.sleep();
    }
}






//        Mat imgHSV;
//        cvtColor(dg_img, imgHSV, COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV
//        Mat imgThresholded;
//        inRange(imgHSV, Scalar(colour.iLowH, colour.iLowS, colour.iLowV),
//                Scalar(colour.iHighH, colour.iHighS, colour.iHighV), imgThresholded); //Threshold the image

//        // erode and dilate to remove noise
//        // morphological opening (remove small objects from the foreground)
//        erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
//        dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
//        // morphological closing (fill small holes in the foreground)
//        dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
//        erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );

//        if (showResults)
//            imshow("threshold",imgThresholded);
//        saveIncomingImage();
//        detectGreenLimitPost();

//namedWindow("Control", CV_WINDOW_AUTOSIZE); //create a window called "Control"

//int iLowH = 75;
//int iHighH = 160;

//int iLowS = 40;
//int iHighS = 255;

//int iLowV = 40;
//int iHighV = 255;

////Create trackbars in "Control" window
//cvCreateTrackbar("LowH", "Control", &iLowH, 179); //Hue (0 - 179)
//cvCreateTrackbar("HighH", "Control", &iHighH, 179);

//cvCreateTrackbar("LowS", "Control", &iLowS, 255); //Saturation (0 - 255)
//cvCreateTrackbar("HighS", "Control", &iHighS, 255);

//cvCreateTrackbar("LowV", "Control", &iLowV, 255); //Value (0 - 255)
//cvCreateTrackbar("HighV", "Control", &iHighV, 255);

//while(true)
//{
////    Mat imgOriginal = cv_ptr->image;
// Mat imgOriginal = imread(testPictureFilename);


// Mat imgHSV;

// cvtColor(imgOriginal, imgHSV, COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV

// Mat imgThresholded;

// inRange(imgHSV, Scalar(iLowH, iLowS, iLowV), Scalar(iHighH, iHighS, iHighV), imgThresholded); //Threshold the image

// //morphological opening (remove small objects from the foreground)
// erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
// dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );

// //morphological closing (fill small holes in the foreground)
// dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
// erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );

// imshow("Thresholded Image", imgThresholded); //show the thresholded image
// imshow("Original", imgOriginal); //show the original image

// if (waitKey(30) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
//{
//   std::cout << "esc key is pressed by user" << std::endl;
//}

//}
