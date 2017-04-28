/*
 * Copyright 2017. All rights reserved.
 * Institute of Measurement and Control Systems
 * Karlsruhe Institute of Technology, Germany
 *
 * authors:
 *  Johannes Graeter (johannes.graeter@kit.edu)
 *  and others
 */
#pragma once

#include <Eigen/Eigen>

namespace motion_estimation_mono{
    namespace normal_prior_cost_functions{

using Tracklet=std::vector<Eigen::Vector<double,2,1>>;
using Tracklets=std::vector<Tracklet>;

// define vector and matrix just as ceres does
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<double,
                      Eigen::Dynamic,
                      Eigen::Dynamic,
                      Eigen::RowMajor>;
/**
* @class EpipolarErrorFunction
* Interface for epipolar error functions that satisfy the normal_prior interface -> defined as http://ceres-solver.org/nnls_modeling.html#normalprior
*/
class EpipolarErrorFunction
{
    public:
        EpipolarErrorFunction (){;}

        /**
         * @brief get error from fundamental matrix and two subsequent matches
         * @param fundamental matrix 
         * @return double error for this match correspondance
         * */
        virtual Eigen::VectorXd get_error(const Eigen::Matrix<double,1,9>& fundamentalMatrix)=0;
        
        /**
        *  @brief precompute stack
        */
        virtual void precompute(const Tracklets& input)=0;
};


/**
*  @class MatchesCostFunctorEpipolar
*  @par
*
*  get cost as epipolar error
*/
class EpipolarErrorFunctionSampsonInliers: public EpipolarErrorFunction
{
public: // public classes/enums/types etc...

public: //attributes
    double InlierThreshold;

    ///@brief store matches for nominator
    Eigen::MatrixXd MatchesStack;
    ///@brief store matches for denominator
    Eigen::MatrixXd MatchesDenominatorStack;

    double err_val;

public: // public methods

    /**
  * default constructor
  */
    EpipolarErrorFunctionSampsonInliers(double InlierThreshold):
        InlierThreshold(InlierThreshold)
    {
        MatchesStack= -1.*Eigen::MatrixXd::Ones(1,1);
        MatchesDenominatorStack= -1.*Eigen::MatrixXd::Ones(1,1);
        err_val=1000.;
    }

    inline double squ(const double& a) const { return a*a;}

    virtual void precompute(const CameraInput& In)
    {
        auto CurMatches=In.get_matches();

        if(CurMatches.size() == 0){ std::cerr<<"attention matches are empty. Functor will return "<<err_val<<std::endl;};

        Eigen::MatrixXd CurMatchesStack(9,CurMatches.size());
        Eigen::MatrixXd CurMatchesDenominatorStack(11,CurMatches.size());

        for (size_t i = 0; i < CurMatches.size(); ++i)
        {
            const auto & Tracklet=CurMatches[i];
            assert(Tracklet.size()>1);
            auto CurPoint=*(Tracklet.rbegin());
            auto LastPoint=*(Tracklet.rbegin()+1);

            Eigen::Vector3d XOld(LastPoint.first, LastPoint.second,1.);

            //assemble the matches as in Hartley and Zissermann for eight point algorithm
            CurMatchesStack.col(i) << XOld*CurPoint.first, XOld*CurPoint.second, XOld;

            
            //assemble matches for the denominator of Sampson Dist(Hartley)
            const double& xc=CurPoint.first;    const double& yc=CurPoint.second;
            const double& xp=LastPoint.first;   const double& yp=LastPoint.second;
            CurMatchesDenominatorStack.col(i) <<    squ(xc),   squ(yc),    xc*yc,  xc, yc,
                                                    squ(xp),   squ(yp),    xp*yp,  xp, yp, 
                                                    1.;
        }

        MatchesDenominatorStack=CurMatchesDenominatorStack;
        MatchesStack=CurMatchesStack;
//        std::cout<<"precomputed matches and denom"<<std::endl;
    }

    Eigen::Matrix<double,1,11> transformed_fundamental_matrix_denominator(const Eigen::Matrix<double,1,9>& f) const
    {
        //f is in form f11,f12,f13, f21,f22,f23, f31,f32,f33
        const double& f11=f(0,0); const double& f12=f(0,1); const double& f13=f(0,2);
        const double& f21=f(0,3); const double& f22=f(0,4); const double& f23=f(0,5);
        const double& f31=f(0,6); const double& f32=f(0,7); const double& f33=f(0,8);

        //this is JJt=squ(Fx[0])+squ(Fx[1])+squ(Ftxp[0])+squ(Ftxp[1]); ->denominator for samspon distance
        Eigen::Matrix<double,1,11> Out;
        Out<<   squ(f11)+squ(f21), squ(f12)+squ(f22), 2.*(f11*f12+f21*f22), 2.*(f11*f13+f21*f23), 2.*(f12*f13+f22*f23), 
                squ(f11)+squ(f12), squ(f21)+squ(f22), 2.*(f11*f21+f12*f22), 2.*(f11*f31+f12*f32), 2.*(f21*f31+f22*f32), 
                squ(f13)+squ(f23)+squ(f31)+squ(f32) ;
        return Out;
    }

    Eigen::VectorXd get_error(const Eigen::Matrix<double,1,9>& FundamentalVec)
    {
        if(MatchesStack.cols()==0) {
            Eigen::VectorXd a(1);
            a<<err_val;
            return a;
        }
        assert(MatchesStack.cols()>0);
        assert(MatchesStack.rows()==9);
        assert(MatchesStack(0,0)!=-1.&& "Matches need to be precomputed");
        assert(MatchesDenominatorStack.cols() > 0);
        assert(MatchesDenominatorStack.rows() == 11);
        assert(MatchesDenominatorStack(0,0)!=1.&& "Matches need to be precomputed");

        assert(FundamentalVec.allFinite());

        //calc epipolar error
        Eigen::VectorXd NominSqrt=(FundamentalVec*MatchesStack).transpose();
        Eigen::VectorXd Denomin=(transformed_fundamental_matrix_denominator(FundamentalVec)*MatchesDenominatorStack).transpose();
        assert((Denomin.array()>1e-12).all());
        assert(Denomin.rows()==NominSqrt.rows());
        
        //elementwise square and then division
        Eigen::VectorXd ErrorVec=NominSqrt.array().square() / Denomin.array();

        //threshold error and return it
        std::for_each(ErrorVec.data(),ErrorVec.data()+ErrorVec.size(),[&](double& a){
                if(a<InlierThreshold) a=0.;
                else a=1.;
                });
        return ErrorVec;
    }
};
}

/**
*  @class MatchesCostFunctorEpipolar
*  @par
*
*  get cost as epipolar error attention: this is wrong! it has to be normalized by the norm of Fx[:2]
*/
class EpipolarErrorFunctionSampsonInliers: public EpipolarErrorFunction
{
public: // public classes/enums/types etc...

public: //attributes
    double InlierThreshold;

    ///@brief store matches for nominator
    Eigen::MatrixXd MatchesStack;
    ///@brief store matches for denominator
    Eigen::MatrixXd MatchesDenominatorStack;

    double err_val;

public: // public methods

    /**
  * default constructor
  */
    EpipolarErrorFunctionSampsonInliers(double InlierThreshold):
        InlierThreshold(InlierThreshold)
    {
        MatchesStack= -1.*Eigen::MatrixXd::Ones(1,1);
        MatchesDenominatorStack= -1.*Eigen::MatrixXd::Ones(1,1);
        err_val=1000.;
    }

    inline double squ(const double& a) const { return a*a;}

    virtual void precompute(const CameraInput& In)
    {
        auto CurMatches=In.get_matches();

        if(CurMatches.size() == 0){ std::cerr<<"attention matches are empty. Functor will return "<<err_val<<std::endl;};

        Eigen::MatrixXd CurMatchesStack(9,CurMatches.size());
        Eigen::MatrixXd CurMatchesDenominatorStack(11,CurMatches.size());

        for (size_t i = 0; i < CurMatches.size(); ++i)
        {
            const auto & Tracklet=CurMatches[i];
            assert(Tracklet.size()>1);
            auto CurPoint=*(Tracklet.rbegin());
            auto LastPoint=*(Tracklet.rbegin()+1);

            Eigen::Vector3d XOld(LastPoint.first, LastPoint.second,1.);

            //assemble the matches as in Hartley and Zissermann for eight point algorithm
            CurMatchesStack.col(i) << XOld*CurPoint.first, XOld*CurPoint.second, XOld;

            
            //assemble matches for the denominator of Sampson Dist(Hartley)
            const double& xc=CurPoint.first;    const double& yc=CurPoint.second;
            const double& xp=LastPoint.first;   const double& yp=LastPoint.second;
            CurMatchesDenominatorStack.col(i) <<    squ(xc),   squ(yc),    xc*yc,  xc, yc,
                                                    squ(xp),   squ(yp),    xp*yp,  xp, yp, 
                                                    1.;
        }

        MatchesDenominatorStack=CurMatchesDenominatorStack;
        MatchesStack=CurMatchesStack;
//        std::cout<<"precomputed matches and denom"<<std::endl;
    }

    Eigen::Matrix<double,1,11> transformed_fundamental_matrix_denominator(const Eigen::Matrix<double,1,9>& f) const
    {
        //f is in form f11,f12,f13, f21,f22,f23, f31,f32,f33
        const double& f11=f(0,0); const double& f12=f(0,1); const double& f13=f(0,2);
        const double& f21=f(0,3); const double& f22=f(0,4); const double& f23=f(0,5);
        const double& f31=f(0,6); const double& f32=f(0,7); const double& f33=f(0,8);

        //this is JJt=squ(Fx[0])+squ(Fx[1])+squ(Ftxp[0])+squ(Ftxp[1]); ->denominator for samspon distance
        Eigen::Matrix<double,1,11> Out;
        Out<<   squ(f11)+squ(f21), squ(f12)+squ(f22), 2.*(f11*f12+f21*f22), 2.*(f11*f13+f21*f23), 2.*(f12*f13+f22*f23), 
                squ(f11)+squ(f12), squ(f21)+squ(f22), 2.*(f11*f21+f12*f22), 2.*(f11*f31+f12*f32), 2.*(f21*f31+f22*f32), 
                squ(f13)+squ(f23)+squ(f31)+squ(f32) ;
        return Out;
    }

    Eigen::VectorXd get_error(const Eigen::Matrix<double,1,9>& FundamentalVec)
    {
        if(MatchesStack.cols()==0) {
            Eigen::VectorXd a(1);
            a<<err_val;
            return a;
        }
        assert(MatchesStack.cols()>0);
        assert(MatchesStack.rows()==9);
        assert(MatchesStack(0,0)!=-1.&& "Matches need to be precomputed");
        assert(MatchesDenominatorStack.cols() > 0);
        assert(MatchesDenominatorStack.rows() == 11);
        assert(MatchesDenominatorStack(0,0)!=1.&& "Matches need to be precomputed");

        assert(FundamentalVec.allFinite());

        //calc epipolar error
        Eigen::VectorXd NominSqrt=(FundamentalVec*MatchesStack).transpose();
        Eigen::VectorXd Denomin=(transformed_fundamental_matrix_denominator(FundamentalVec)*MatchesDenominatorStack).transpose();
        assert((Denomin.array()>1e-12).all());
        assert(Denomin.rows()==NominSqrt.rows());
        
        //elementwise square and then division
        Eigen::VectorXd ErrorVec=NominSqrt.array().square() / Denomin.array();

        //threshold error and return it
        std::for_each(ErrorVec.data(),ErrorVec.data()+ErrorVec.size(),[&](double& a){
                if(a<InlierThreshold) a=0.;
                else a=1.;
                });
        return ErrorVec;
    }
}


}
}
