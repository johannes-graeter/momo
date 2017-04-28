// google test docs
// wiki page: https://code.google.com/p/googletest/w/list
// primer: https://code.google.com/p/googletest/wiki/V1_7_Primer
// FAQ: https://code.google.com/p/googletest/wiki/FAQ
// advanced guide: https://code.google.com/p/googletest/wiki/V1_7_AdvancedGuide
// samples: https://code.google.com/p/googletest/wiki/V1_7_Samples
//
// List of some basic tests fuctions:
// Fatal assertion                      Nonfatal assertion                   Verifies / Description
//-------------------------------------------------------------------------------------------------------------------------------------------------------
// ASSERT_EQ(expected, actual);         EXPECT_EQ(expected, actual);         expected == actual
// ASSERT_NE(val1, val2);               EXPECT_NE(val1, val2);               val1 != val2
// ASSERT_LT(val1, val2);               EXPECT_LT(val1, val2);               val1 < val2
// ASSERT_LE(val1, val2);               EXPECT_LE(val1, val2);               val1 <= val2
// ASSERT_GT(val1, val2);               EXPECT_GT(val1, val2);               val1 > val2
// ASSERT_GE(val1, val2);               EXPECT_GE(val1, val2);               val1 >= val2
//
// ASSERT_FLOAT_EQ(expected, actual);   EXPECT_FLOAT_EQ(expected, actual);   the two float values
// are almost equal (4 ULPs)
// ASSERT_DOUBLE_EQ(expected, actual);  EXPECT_DOUBLE_EQ(expected, actual);  the two double values
// are almost equal (4 ULPs)
// ASSERT_NEAR(val1, val2, abs_error);  EXPECT_NEAR(val1, val2, abs_error);  the difference between
// val1 and val2 doesn't exceed the given absolute error
//
// Note: more information about ULPs can be found here:
// http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
//
// Example of two unit test:
// TEST(Math, Add) {
//    ASSERT_EQ(10, 5+ 5);
//}
//
// TEST(Math, Float) {
//	  ASSERT_FLOAT_EQ((10.0f + 2.0f) * 3.0f, 10.0f * 3.0f + 2.0f * 3.0f)
//}
//=======================================================================================================================================================
#include <Eigen/Eigen>
#include "camera.hpp"
#include "cost_functions.hpp"
#include "cost_functors_ceres.hpp"
#include "gtest/gtest.h"

Eigen::Vector3d back_project(Eigen::Vector3d p, Eigen::Matrix3d intrin) {
    p /= p[2];
    return intrin * p;
}

TEST(Momo, motion_models) {
    // test identities
    {
        momo::MotionModelOneWheel2d<double> model_2d(0., 0.);
        //        std::cout << "model_2d=\n" << model_2d.ToPose().matrix() << std::endl;

        ASSERT_EQ(model_2d.ToPose().isApprox(Eigen::Affine3d::Identity()), true);

        momo::MotionModelOneWheel3d<double> model_3d(0., 0., 0., 0.);
        //        std::cout << "model_3d=\n" << model_3d.ToPose().matrix() << std::endl;

        ASSERT_EQ(model_3d.ToPose().isApprox(Eigen::Affine3d::Identity()), true);
    }
    // test real movement
    {
        // input movement
        std::vector<double> roll_pitch_yaw{0., 0., 3. * M_PI / 180.};
        Eigen::Vector3d translation(2.99897, 0.0785308, 0.);

        // motion in vehicle cos
        Eigen::Affine3d m(Eigen::AngleAxisd(roll_pitch_yaw[2], Eigen::Vector3d::UnitZ()));
        m.translation() = translation;
        //        std::cout << "input movement=\n" << m.matrix() << std::endl;

        // calc arclength
        double secant = m.translation().norm();
        // secant=2*radius*sin(delta_yaw/2), arc_length=radius*delta_yaw ->
        // arc_length=secant/2.*delta_yaw/sin(delta_yaw/2.)
        // https://en.wikipedia.org/wiki/Circular_segment
        double yaw = roll_pitch_yaw[2];
        double arc_length = secant / 2. * yaw / std::sin(yaw / 2.);
        //        std::cout << "arc_length=" << arc_length << std::endl;

        // execute models
        momo::MotionModelOneWheel2d<double> model_2d(arc_length, yaw);
        //        std::cout << "model_2d=\n" << model_2d.ToPose().matrix() << std::endl;

        ASSERT_EQ(model_2d.ToPose().isApprox(m, 0.001), true);

        momo::MotionModelOneWheel3d<double> model_3d(arc_length, 0., 0., yaw);
        //        std::cout << "model_3d=\n" << model_3d.ToPose().matrix() << std::endl;

        ASSERT_EQ(model_3d.ToPose().isApprox(m, 0.001), true);
    }
}

TEST(Momo, cost_functors_ceres_no_movement) {
    using CostFuncType = momo::commons::CostFunctionType;

    // motion in vehicle cos
    std::vector<double> roll_pitch_yaw{0., 0., 0.};
    Eigen::Affine3d m(Eigen::AngleAxisd(roll_pitch_yaw[2], Eigen::Vector3d::UnitZ()));
    m.translation() = Eigen::Vector3d(0., 0., 0.);

    // convert vehicle motion to camera motion (camera views to the left)
    Eigen::Affine3d vehicle_to_cam(Eigen::AngleAxisd(-M_PI / 2., Eigen::Vector3d::UnitX()));

    // point to project in camera cos
    Eigen::Vector3d p_camera(1.5, 1., 3.);

    //    std::cout << "original motion=\n" << m.matrix() << std::endl;
    //    std::cout << "motion camera=\n" << m_camera.matrix() << std::endl;

    Eigen::Matrix3d intrin;
    intrin << 300., 0., 200., 0., 300., 100., 0., 0., 1.;

    //    std::cout << "intrinsics=\n" << intrin << std::endl;

    // instantiate functor
    auto camera = std::make_shared<momo::Camera>(vehicle_to_cam, intrin);

    auto match = std::make_pair(back_project(p_camera, intrin), back_project(p_camera, intrin));

    //    std::cout << "m_prev=" << match.first.transpose() << std::endl;
    //    std::cout << "m_cur=" << match.second.transpose() << std::endl;

    // calc arclength
    double arc_length = 0.001;

    // per default we evaluate the sampson distance, can be any of CostFunctionType. For symmetric
    // epipolar distance forward and backward epipolar distance have to be combined
    momo::CostFunctorCeresRollPitchYaw functor_sampson(match, arc_length, camera,
                                                       CostFuncType::SampsonDistance);
    momo::CostFunctorCeresRollPitchYaw functor_epi_forw(match, arc_length, camera,
                                                        CostFuncType::EpipolarDistanceForward);
    momo::CostFunctorCeresRollPitchYaw functor_epi_backw(match, arc_length, camera,
                                                         CostFuncType::EpipolarDistanceBackward);
    // non-ceres evaluation test with ideal values
    std::vector<double> res_sampson{1000.};
    std::vector<double> res_epi_forw{1000.};
    std::vector<double> res_epi_backw{1000.};

    // evaluate functors
    bool ret_sampson = functor_sampson(roll_pitch_yaw.data(), res_sampson.data());
    bool ret_epi_forw = functor_epi_forw(roll_pitch_yaw.data(), res_epi_forw.data());
    bool ret_epi_backw = functor_epi_backw(roll_pitch_yaw.data(), res_epi_backw.data());

    ASSERT_EQ(ret_sampson, true);
    ASSERT_EQ(ret_epi_forw, true);
    ASSERT_EQ(ret_epi_backw, true);

    std::cout << "residuals: sampson=" << res_sampson[0] << " forward=" << res_epi_forw[0]
              << " backward=" << res_epi_backw[0] << std::endl;

    ASSERT_DOUBLE_EQ(res_sampson[0], 0.);
    ASSERT_DOUBLE_EQ(res_epi_forw[0], 0.);
    ASSERT_DOUBLE_EQ(res_epi_backw[0], 0.);
}

std::tuple<bool, bool, bool, bool, bool, double, double, double, double, double>
execute_cost_functor_test(std::vector<double> roll_pitch_yaw, double arc_length) {
    using CostFuncType = momo::commons::CostFunctionType;

    // motion in vehicle cos
    momo::MotionModelOneWheel3d<double> model_3d(arc_length, roll_pitch_yaw[0], roll_pitch_yaw[1],
                                                 roll_pitch_yaw[2]);
    Eigen::Affine3d m = model_3d.ToPose();

    // convert vehicle motion to camera motion (camera views to the left)
    Eigen::Affine3d vehicle_to_cam(Eigen::AngleAxisd(-M_PI / 2., Eigen::Vector3d::UnitX()));
    Eigen::Affine3d m_camera = vehicle_to_cam.inverse() * m * vehicle_to_cam;

    // point to project in camera cos
    Eigen::Vector3d p_camera(1., 1., 1.);

    Eigen::Matrix3d intrin;
    intrin << 300., 0., 200., 0., 300., 100., 0., 0., 1.;

    // instantiate functor
    auto camera = std::make_shared<momo::Camera>(vehicle_to_cam.inverse(), intrin);

    Eigen::Vector3d transformed_point = m_camera.inverse() * p_camera;
    //    std::cout << "non-trasnformed point=" << p_camera.transpose() << std::endl;
    //    std::cout << "trasnformed point=" << transformed_point.transpose() << std::endl;
    auto match =
        std::make_pair(back_project(p_camera, intrin), back_project(transformed_point, intrin));

    auto match_geometric =
        std::make_pair(p_camera / p_camera.norm(), transformed_point / transformed_point.norm());

    // per default we evaluate the sampson distance, can be any of CostFunctionType. For symmetric
    // epipolar distance forward and backward epipolar distance have to be combined
    momo::CostFunctorCeresRollPitchYaw functor_sampson(match, arc_length, camera,
                                                       CostFuncType::SampsonDistance);
    momo::CostFunctorCeresRollPitchYaw functor_epi_forw(match, arc_length, camera,
                                                        CostFuncType::EpipolarDistanceForward);
    momo::CostFunctorCeresRollPitchYaw functor_epi_backw(match, arc_length, camera,
                                                         CostFuncType::EpipolarDistanceBackward);
    momo::CostFunctorCeresRollPitchYaw functor_geometric_linear(
        match_geometric, arc_length, camera, CostFuncType::RayToEpipolarPlaneLinear);
    momo::CostFunctorCeresRollPitchYaw functor_geometric_quadratic(
        match_geometric, arc_length, camera, CostFuncType::RayToEpipolarPlaneQuadratic);

    // non-ceres evaluation test with ideal values
    std::vector<double> res_sampson{1000.};
    std::vector<double> res_epi_forw{1000.};
    std::vector<double> res_epi_backw{1000.};
    std::vector<double> res_geometric_linear{1000.};
    std::vector<double> res_geometric_quadratic{1000.};

    // evaluate functors
    bool ret_sampson = functor_sampson(roll_pitch_yaw.data(), res_sampson.data());
    bool ret_epi_forw = functor_epi_forw(roll_pitch_yaw.data(), res_epi_forw.data());
    bool ret_epi_backw = functor_epi_backw(roll_pitch_yaw.data(), res_epi_backw.data());
    bool ret_geometric_linear =
        functor_geometric_linear(roll_pitch_yaw.data(), res_geometric_linear.data());
    bool ret_geometric_quadratic =
        functor_geometric_quadratic(roll_pitch_yaw.data(), res_geometric_quadratic.data());

    return std::make_tuple(ret_sampson, ret_epi_forw, ret_epi_backw, ret_geometric_linear,
                           ret_geometric_quadratic, res_sampson[0], res_epi_forw[0],
                           res_epi_backw[0], res_geometric_linear[0], res_geometric_quadratic[0]);
}

TEST(Momo, cost_functors_ceres) {
    std::vector<std::pair<std::vector<double>, double>> movements;
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -0.0001 / 180. * M_PI}, 0.1));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -0.1 / 180. * M_PI}, 0.1));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -1. / 180. * M_PI}, 0.1));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -3. / 180. * M_PI}, 0.1));

    movements.push_back(std::make_pair(std::vector<double>{0., 0., -0.0001 / 180. * M_PI}, 1.5));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -0.1 / 180. * M_PI}, 1.5));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -1. / 180. * M_PI}, 1.5));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -3. / 180. * M_PI}, 1.5));

    movements.push_back(std::make_pair(std::vector<double>{0., 0., -0.0001 / 180. * M_PI}, 3));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -0.1 / 180. * M_PI}, 3));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -1. / 180. * M_PI}, 3));
    movements.push_back(std::make_pair(std::vector<double>{0., 0., -3. / 180. * M_PI}, 3));

    for (const auto& el : movements) {
        bool ret_sampson, ret_epi_forw, ret_epi_backw, ret_geometric_linear,
            ret_geometric_quadratic;
        double res_sampson, res_epi_forw, res_epi_backw, res_geometric_linear,
            res_geometric_quadratic;
        std::tie(ret_sampson, ret_epi_forw, ret_epi_backw, ret_geometric_linear,
                 ret_geometric_quadratic, res_sampson, res_epi_forw, res_epi_backw,
                 res_geometric_linear, res_geometric_quadratic) =
            execute_cost_functor_test(el.first, el.second);

        //        momo::MotionModelOneWheel3d<double> model_3d_2(1.5, 0., 0., -3. / 180. * M_PI);
        //        std::cout << "model_3d, -3 degree=\n" << model_3d_2.ToPose().matrix() <<
        //        std::endl;

        ASSERT_EQ(ret_sampson, true);
        ASSERT_EQ(ret_epi_forw, true);
        ASSERT_EQ(ret_epi_backw, true);
        ASSERT_EQ(ret_geometric_linear, true);
        ASSERT_EQ(ret_geometric_quadratic, true);

        std::cout <<std::scientific<<std::setprecision(3)<< "residuals: sampson=" << res_sampson << " forward=" << res_epi_forw
                  << " backward=" << res_epi_backw << " geometric linear=" << res_geometric_linear
                  << " geometric quadratic=" << res_geometric_quadratic << std::endl;

        ASSERT_NEAR(res_sampson, 0., 1e-10);
        ASSERT_NEAR(res_epi_forw, 0., 1e-10);
        ASSERT_NEAR(res_epi_backw, 0., 1e-10);
        ASSERT_NEAR(res_geometric_linear, 0., 1e-10);
        ASSERT_NEAR(res_geometric_quadratic, 0., 1e-10);
    }
}
// A google test function (uncomment the next function, add code and
// change the names TestGroupName and TestName)
// TEST(TestGroupName, TestName) {
// TODO: Add your test code here
//}
