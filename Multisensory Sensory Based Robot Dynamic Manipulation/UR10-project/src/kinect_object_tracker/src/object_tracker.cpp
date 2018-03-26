#include "object_tracker.h"
#include <iostream>
#include <cmath>


ObjectTracker::ObjectTracker()
{

    imgSub = node.subscribe("/kinect2/qhd/image_color_rect",10, &ObjectTracker::imageCallback, this);
    pointSub = node.subscribe("/kinect2/qhd/points",10, &ObjectTracker::pointsCallback, this);
    hsvPub = node.advertise<sensor_msgs::Image>("/objectTracker/hsv_img", 1);
    colorPub = node.advertise<sensor_msgs::Image>("/objectTracker/color_img", 1);
    circlePub = node.advertise<sensor_msgs::Image>("/objectTracker/circle_img", 1);
    positionPub = node.advertise<geometry_msgs::Point>("/objectTracker/position", 1);

    measurmentFrame = 0;

    


}

void ObjectTracker::imageCallback(const sensor_msgs::Image::ConstPtr &img){
    image  = img;
    cv_ptr = cv_bridge::toCvCopy(img,sensor_msgs::image_encodings::BGR8);
    // int width = image->width;
    // int height= image->height;
    // ROS_INFO("Image : width, height: [%i][%i]", width, height);

}


void ObjectTracker::pointsCallback(const sensor_msgs::PointCloud2::ConstPtr& cloud){

    // if(measurmentFrame < 5){
        pcloud.header = cloud->header; 
        pcloud.height = cloud->height;
        pcloud.width = cloud->width;
        pcloud.fields = cloud->fields;
        pcloud.is_bigendian = cloud->is_bigendian;
        pcloud.point_step = cloud->point_step;
        pcloud.row_step = cloud->row_step;
        pcloud.data = cloud->data;
        pcloud.is_dense = cloud->is_dense;
        
        int arrayPosition = 0;
        int arrayPosX = 0;
        int arrayPosY = 0;
        int arrayPosZ = 0;

        // std::vector<float> xData;
        // std::vector<float> yData;
        // std::vector<float> zData;

        float tempx = 0.0;
        float tempy = 0.0;
        float tempz = 0.0;

        // for(int x = largestRect.tl().x; x <= largestRect.br().x; x++){
        //     for(int y = largestRect.tl().y; y <= largestRect.br().y; y++){
        //         // Get Arry Position of that pixel data
        //         arrayPosition = y * pcloud.row_step + x * pcloud.point_step;

        //         //Get x,y,z arry position of that pixel
        //         // arrayPosX = arrayPosition + pcloud.fields[0].offset; // X has an offset of 0
        //         // arrayPosY = arrayPosition + pcloud.fields[1].offset; // Y has an offset of 4
        //         arrayPosZ = arrayPosition + pcloud.fields[2].offset; // Z has an offset of 8
                



        //         // std::memcpy(&tempx, &pcloud.data[arrayPosX], sizeof(float));
        //         // std::memcpy(&tempy, &pcloud.data[arrayPosY], sizeof(float));
        //         std::memcpy(&tempz, &pcloud.data[arrayPosZ], sizeof(float));


        //         // xData.push_back(tempx);
        //         // yData.push_back(tempy);
        //         zData.push_back(tempz);
        //     }
        // }

        int x = centerRect.x;
        int y = centerRect.y;

        arrayPosition = y * pcloud.row_step + x * pcloud.point_step;

        //Get x,y,z arry position of that pixel
        arrayPosX = arrayPosition + pcloud.fields[0].offset; // X has an offset of 0
        arrayPosY = arrayPosition + pcloud.fields[1].offset; // Y has an offset of 4
        arrayPosZ = arrayPosition + pcloud.fields[2].offset;
        // arrayPosZ = arrayPosition + pcloud.fields[2].offset; // Z has an offset of 8



        std::memcpy(&tempx, &pcloud.data[arrayPosX], sizeof(float));
        std::memcpy(&tempy, &pcloud.data[arrayPosY], sizeof(float));
        std::memcpy(&tempz, &pcloud.data[arrayPosZ], sizeof(float));
        // std::memcpy(&tempz, &pcloud.data[arrayPosZ], sizeof(float));


        // xData.push_back(tempx);
        // yData.push_back(tempy);
        // zData.push_back(tempz);

        // // Compute avarage x and y data
        // float x = std::accumulate( xData.begin(), xData.end(), 0.0)/xData.size(); 
        // float y = std::accumulate( yData.begin(), yData.end(), 0.0)/yData.size(); 

        // // filter z data and compute avarage
        // std::sort(zData.begin(), zData.end());
        // zData.resize(1);
        // tempz = zData[0];//std::accumulate( zData.begin(), zData.end(), 0.0)/zData.size();

        // xvec.push_back(tempx);
        // yvec.push_back(tempy);
        // zvec.push_back(tempz);

        // measurmentFrame++;

        // std::cout<< measurmentFrame<<std::endl;
    
    //}else{

        
        // auto it1 = std::remove_if(xvec.begin(), xvec.end(), [](float d) {return std::isnan(d);});
        // auto it2 = std::remove_if(yvec.begin(), yvec.end(), [](float d) {return std::isnan(d);});
        // auto it3 = std::remove_if(zvec.begin(), zvec.end(), [](float d) {return std::isnan(d);});


        // xvec.erase(it1, xvec.end());
        // yvec.erase(it2, yvec.end());
        // zvec.erase(it3, zvec.end());

        // for (auto i: xvec)
        //     std::cout << i << ' '<<std::endl;
        //     std::cout<<"-----------------"<<std::endl;
        // for (auto i: yvec)
        //     std::cout << i << ' ';
        // for (auto i: zvec)
        //     std::cout << i << ' ';

        // float x = std::accumulate( xvec.begin(), xvec.end(), 0.0)/xvec.size();
        // float y = std::accumulate( yvec.begin(), yvec.end(), 0.0)/yvec.size();
        // float z = std::accumulate( zvec.begin(), zvec.end(), 0.0)/zvec.size();

        // Median Filter
        // std::sort(xvec.begin(),xvec.end());
        // std::sort(yvec.begin(),yvec.end());
        // std::sort(zvec.begin(),zvec.end());
        // float x = xvec[2];
        // float y = yvec[2];
        // float z = zvec[2];

        // ROS_INFO("x:[%f],y:[%f],z:[%f]",x,y,z);

        // empty vectors
        // xvec.clear();
        // yvec.clear();
        // zvec.clear();


        // if z is feasable -> publish
        if(tempz < 2){
            // transform.setOrigin( tf::Vector3(x, y, z) );
            // tf::Quaternion q;
            // q.setRPY(0, 0, 0);
            // transform.setRotation(q);

            // position.x = x;
            // position.y = y;
            // position.z = z;
            meas[0] = tempx;
            meas[1] = tempy;
            meas[2] = tempz;
            
        }



        // measurmentFrame = 0;

        
    // }

}

void ObjectTracker::findObjectPosition(int meas_count){

    EMA(filtered_pos, meas, 0.5, meas_count);
    
    position.x = filtered_pos[0];
    position.y = filtered_pos[1];
    position.z = filtered_pos[2];

    positionPub.publish(position);

}

void ObjectTracker::EMA(Eigen::Vector3d &filtered_pos, Eigen::Vector3d &meas, float alpha, int count){

        // ema
        filtered_pos = alpha*filtered_pos + (1-alpha)*meas;

        if(count < 100){
            filtered_pos = filtered_pos/(1-pow(alpha,count));
        }

}



void ObjectTracker::findObjectContour(){
    cv::Mat canny_output;
    std::vector<std::vector<cv::Point> > contours;
    std::vector<cv::Vec4i> hierarchy;
    if(cv_ptr)
        {
            // Image Conversions
            cv::cvtColor(cv_ptr->image, hsvImg, CV_BGR2HSV);

            // Search for object in HSV space
            cv::inRange(hsvImg, cv::Scalar(0, 100, 50), cv::Scalar(10, 255, 255), redImgLower);
            cv::inRange(hsvImg, cv::Scalar(160, 100, 50), cv::Scalar(179, 255, 255), redImgUpper);

            // Add all found pixels
            cv::addWeighted(redImgLower, 1.0, redImgUpper, 1.0, 0.0, redHueImg);

            // Filter noise
            cv::medianBlur (redHueImg, redHueImg, 27 );     
            // cv::GaussianBlur(redHueImg, redHueImg, cv::Size(9, 9), 0, 0);
            //cv::erode(redHueImg, redHueImg, 0);         // Erode Filter Effect
            //cv::dilate(redHueImg, redHueImg, 0);        // Dilate Filter Effect

            
            /// Detect edges using canny
            Canny(redHueImg, canny_output, 100, 200, 3);
            /// Find contours
            findContours( canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0));


            /// Approximate contours to polygons + get bounding rects and circles
            std::vector<std::vector<cv::Point>> contours_poly( contours.size() );
            std::vector<cv::Rect> boundRect( contours.size() );
            std::vector<cv::Point2f> center( contours.size() );
            std::vector<float> radius( contours.size() );

            int largest_area=0;
            int largest_contour_index=0;
            double a;

            for( int i = 0; i < contours.size(); i++ )
                {     
                approxPolyDP( cv::Mat(contours[i]), contours_poly[i], 3, true );
                boundRect[i] = boundingRect( cv::Mat(contours_poly[i]) );
                a = cv::contourArea( contours[i],false); 
                if(a>largest_area){
                    largest_area=a;
                    largest_contour_index=i;                //Store the index of largest contour
                    largestRect=boundingRect(cv::Mat(contours_poly[i])); // Find the bounding rectangle for biggest contour
                 }
                //minEnclosingCircle( cv::Mat(contours_poly[i]), center[i], radius[i] );
                }
            /// Draw polygonal contour + bonding rects + circles
            if(largest_area){
                centerRect = (largestRect.br() + largestRect.tl())*0.5;
            }else{
                centerRect.x = 0;
                centerRect.y = 0;
            }
            
            // ROS_INFO("RectCenter:[%i,%i]", centerRect.x, centerRect.y);


            //drawContours( cv_ptr->image, contours_poly, i, cv::Scalar(0, 255, 0), 1, 8, std::vector<cv::Vec4i>(), 0, cv::Point() );
            rectangle( cv_ptr->image, largestRect.tl(), largestRect.br(), cv::Scalar(0, 255, 0), 2, 8, 0 );
            cv::circle(cv_ptr->image, centerRect, 5, cv::Scalar(0, 0, 0), -1);
            // for( int i = 0; i< contours.size(); i++ )
            // {
            //     cv::Moments mu = cv::moments(contours[i]);
            //     cv::Point centroid = cv::Point (mu.m10/mu.m00 , mu.m01/mu.m00);
            //     cv::circle(cv_ptr->image, centroid, 5, cv::Scalar(0, 0, 0), -1);
            //     drawContours(cv_ptr->image, contours, i, cv::Scalar(0, 255, 0), 5, 8, hierarchy, 0, cv::Point(0,0));
            // }



            hsv_msg.header = cv_ptr->header;
            hsv_msg.encoding = sensor_msgs::image_encodings::MONO8;
            hsv_msg.image = redHueImg;

            hsvPub.publish(hsv_msg.toImageMsg());
            circlePub.publish(cv_ptr);

        }
}

void ObjectTracker::run(){

    ros::Rate rate(10);

     int count = 0;

    while (ros::ok())
    {
        count++;
        findObjectContour();
        findObjectPosition(count);
        //br->sendTransform(tf::StampedTransform(transform, ros::Time::now(), "kinect2_link", "imu"));
        ros::spinOnce(); // Need to call this function often to allow ROS to process incoming messages
        rate.sleep();
    }
}
