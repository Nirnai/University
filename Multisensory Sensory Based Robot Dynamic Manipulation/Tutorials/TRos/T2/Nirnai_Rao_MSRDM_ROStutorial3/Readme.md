# Homework 2

To build this project, copy the **src** directory into a working catkin worspace. To run it, use the provided launch file.

## Questions

1. D-H convention does not need to be followed in the URDF file. The coordinate frames depend on the orientation of the mesh.
2. URDF is a format to describe links and joints of a robot and visualize them. Xacro is a template for a URDF. The actual URDF code is created during runtime.
3. 
4. The Endeffector position in RViz was not the same as in matlab, since the sizes of the joint were not considered. This was corrected in the urdf.
