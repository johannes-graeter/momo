# Momo

* Momo is a mono-motion estimation library.
* Instead of calculating the full 5d movement (without scale) from matches, different motion-model based motion hypotheses are transformed into the camera coordinate system and the expected correspondences are compared with the measured correspondences
* Interpreting matching errors as residual blocks we can formulate the motion estimation as a non-linear optimization problem
* We solve it efficiently with google's ceres-solver
* Advantages:
    * We can use robustification by loss-functions, no need for RANSAC! 
    * We can use different error metrics (reprojection error, sampson distance, angle between epipolar plane and viewing ray), accounting for non-linear camera models
    * We can use any number of cameras we want, including all matches into one single optimization problem
    * We can reduce the parameter space so that hypothesis that are not feasible will be ignored (f.e. the main part of the image consists of a moving object).
    * We can easily implement new motion models
   
* The corresponding paper and additional information can be found in ./docu
* Videos of the system in action can be found in ./video, use git-lfs  

## Install
### Requirements
* Eigen3: 
<code>sudo apt-get install libeigen3-dev</code>
* catkin: 
    - follow the instructions on [http://wiki.ros.org/catkin](http://wiki.ros.org/catkin) or install ros
* ceres: 
    - follow the instructions on [http://ceres-solver.org/installation.html](http://ceres-solver.org/installation.html)
* googletest for unittests:
    - <code>sudo apt-get install libgtest-dev</code>

### installation
* initiate a catkin workspace:
    - <code>cd *your_catkin_workspace*</code>
    - <code>cd *your_catkin_workspace*/src</code>
    - <code>catkin_init_workspace</code>
* clone mrt_cmake_modules into src of workspace:
    - <code>cd *your_catkin_workspace*/src</code>
    - <code>git clone https://github.com/KIT-MRT/mrt_cmake_modules.git</code>
* clone momo into src of workspace:
    - <code>cd *your_catkin_workspace*/src</code>
    - <code>git clone https://github.com/johannes-graeter/momo.git</code>
* build it with catkin:
    - <code>cd *your_catkin_workspace*</code>
    - <code>catkin_make</code>
* unittests:
    - uncomment **<test_depend>gtest</test_depend>** in package.xml to activate unittests
    - <code>cd *your_catkin_workspace*</code>
    - <code>catkin_make run_tests</code>

* tested with docker ros image

## Usage
* Loss functions:
    * have a look on the [ceres header](https://github.com/kashif/ceres-solver/blob/master/include/ceres/loss_function.h)
    * We use cauchy loss, which downweights the residuals according to the cauchy distribution, see [wikipedia](https://en.wikipedia.org/wiki/Cauchy_distribution) or [wolfram alpha](http://mathworld.wolfram.com/CauchyDistribution.html)
* Implemented error functions:
    * Errors in image space:
        * Epipolar errors: 
            * Forward: "Normal" epipolar error, distance of image point to epipolar linear form old image to current image
            * Backward: Same as forward but from current image to old image
            * Symmetric: Adding both as residual blocks results in the symmetric epipolar error as in [Multiple View Geometry](http://www.robots.ox.ac.uk/~az/tutorials/tutoriala.pdf)
        * Sampson distance: First order approximation of the geometric error as defined in [Multplie View Geometry](http://www.robots.ox.ac.uk/~az/tutorials/tutoriala.pdf) (called Zhang's distance)
    * Errors in 3d space:
        * Ray to epipolar plane:
            * linear: This error is equal to sin(alpha) where alpha is the angle between the epipolar plane and the viewing ray. We call it "linear" since the sine is approximately linear around zero. Simulation in momo_simulation_tool show that normally distributed noise in the image space results in a cauchy-distribution-like looking distribution would be a fitting explanation according to [wolfram alpha](http://mathworld.wolfram.com/CauchyDistribution.html) 
            * quadratic: This error is equal to cos(alpha). We call it "quadratic" since the cosine is approximately quadratic around zero.  
## Limitations
* The camera model must be a single-view-point model
* All cameras must be calibrated to one point on the robot
* Motion-models implemented:
    * arc-length as control input, yaw as parameter (2d single-wheel)
    * arc-length as control, input, roll, pitch, yaw as paramters (3d single-wheel)

## Example
* this is exemplary code for executing momo, comments in capital letters are instructions you have to implement
```cpp
#include <momo/motion_optimizer.hpp>
#include <momo/commons.hpp>

Eigen::Vector3d match_to_viewing_ray(cv::Point a, Eigen::Matrix3d intrinsics){
    Eigen::Vector3d a_hom(a.x,a.y,1.);
    Eigen::Vector3d out=intrinsics.inverse()*a_hom;
    out.normalize()
    return out;
}
int main(int argc, char* argv[])
{
// dummy defines
// REPLACE F,CU,CV BY CAMERA INTRINSICS
// intrinsics
Eigen::Matrix3d intrinsics;
intrinsics<<f,0.,cu,0.,f,cv,0.,0.,1.;

// extrinsics set to zero
// REPLACE BY POSE FROM CAMERA TO MOTION CENTER
Eigen::Affine3d extrinsics = Eigen::Affine3d::Identity();

// camera storage class
momo::Camera cam(extrinsics, intrinsics);

// GET INPUT IMAGES I0 AND I1
// MATCH THEM WITH F.E. OPENCV AND STORE IN opencv_matches

// calculate viewing rays 
momo::MotionOptimizer::Matches input_matches;
for(const auto& el:opencv_matches){
    auto current_viewing_ray = match_to_viewing_ray(el[1], intrinsics);
    auto previous_viewing_ray = match_to_viewing_ray(el[0], intrinsics);
    input_matches.push_back(momo::MotionOptimizer::Match(previous_viewing_ray,current_viewing_ray));
}
momo::MotionOptimizer optimizer;

optimizer.AddData(input_matches, std::make_shared<momo::Camera>(cam), momo::commons::CostFunctionType::RayToEpipolarPlaneLinear)

// GET SCALE PRIOR FROM F.E. ODOMETRY AND STORE IN scale_prior

// calculate motion
Eigen::Affine3d current_motion; 
if (scale_prior > 0.01) {
        bool verbose = true;
        current_motion = optimizer.Solve(scale_prior, verbose);
}

}
```
## Credits

Johannes Gräter

## License
This software is under the LGPLv3 license
