cmake_minimum_required(VERSION 2.8.3)
project(msd3dof_description)

find_package(catkin REQUIRED
    roscpp
    std_msgs
    sensor_msgs
    tf
    tf_conversions
    visualization_msgs
)

##############################################################################
# Eigen
##############################################################################
find_package(Eigen3 REQUIRED)
add_definitions(${EIGEN_DEFINITIONS})
include_directories(${EIGEN3_INCLUDE_DIR})

catkin_package(
        CATKIN_DEPENDS roscpp std_msgs tf tf_conversions visualization_msgs
)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${roscpp_INCLUDE_DIRS}
)

# OPTIONAL: You can create a node to visualize the spring!
add_executable(spring_node src/testMeshScaling.cpp)
target_link_libraries(spring_node 
		  	${catkin_LIBRARIES} 
			${roscpp_LIBRARIES})


foreach(dir meshes urdf)
    install(DIRECTORY ${dir}/
        DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/${dir})
endforeach(dir)
