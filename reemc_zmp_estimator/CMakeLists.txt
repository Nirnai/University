cmake_minimum_required(VERSION 2.8.3)
project(reemc_zmp_estimator)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  realtime_tools
  controller_interface
)

catkin_package()

include_directories(
  ${catkin_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME} src/zmp_controller.cpp)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

# install(TARGETS ${PROJECT_NAME}
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

# foreach(dir config launch)
#   install(DIRECTORY ${dir}
#     DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
#   )
# endforeach()

# install(FILES ros_tutorial_controller_plugins.xml
#         DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
