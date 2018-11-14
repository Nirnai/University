# Important commands 

rosrun image_view image_view image:=/camera/rgb/image_color

# iai Options

roslaunch kinect2_bridge kinect2_bridge.launch [options:=value]
base_name:=<string>
    default: kinect2
    info:    set base name for all topics
sensor:=<string>
    default:
    info:    serial of the sensor to use
fps_limit:=<double>
    default: -1.0
    info:    limit the frames per second
calib_path:=<string>
    default: /home/wiedemeyer/work/src/iai_kinect2/kinect2_bridge/data/
    info:    path to the calibration files
use_png:=<bool>
    default: false
    info:    Use PNG compression instead of TIFF
jpeg_quality:=<int>
    default: 90
    info:    JPEG quality level from 0 to 100
png_level:=<int>
    default: 1
    info:    PNG compression level from 0 to 9
depth_method:=<string>
    default: cuda
    info:    Use specific depth processing: default, cpu, opengl, opencl, cuda, clkde, cudakde
depth_device:=<int>
    default: -1
    info:    openCL device to use for depth processing
reg_method:=<string>
    default: opencl
    info:    Use specific depth registration: default, cpu, opencl
reg_device:=<int>
    default: -1
    info:    openCL device to use for depth registration
max_depth:=<double>
    default: 12.0
    info:    max depth value
min_depth:=<double>
    default: 0.1
    info:    min depth value
queue_size:=<int>
    default: 2
    info:    queue size of publisher
bilateral_filter:=<bool>
    default: true
    info:    enable bilateral filtering of depth images
edge_aware_filter:=<bool>
    default: true
    info:    enable edge aware filtering of depth images
publish_tf:=<bool>
    default: false
    info:    publish static tf transforms for camera
base_name_tf:=<string>
    default: as base_name
    info:    base name for the tf frames
worker_threads:=<int>
    default: 4
    info:    number of threads used for processing the images
