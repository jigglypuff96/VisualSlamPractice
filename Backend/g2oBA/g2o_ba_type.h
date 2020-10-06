#ifndef G2O_BA_TYPE_H
#define G2O_BA_TYPE_H

#include <iostream>
using namespace std;

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>

#include "ceres/ceres.h"

typedef Eigen::Matrix<double, 9, 1> Vector9d;
using Eigen::Vector3d;
using Eigen::Vector2d;

#include "tools.h"

class VertexCameraBAL : public g2o::BaseVertex<9, Vector9d> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	VertexCameraBAL() {}
	virtual bool read ( istream& is ) {return false;}
	virtual bool write ( ostream& os ) const {return false;}
	virtual void setToOriginImpl() {}
	
	virtual void oplusImpl ( const double* update ) {
		Vector9d::ConstMapType v ( update );
		_estimate += v;
	}
};

class VertexPointBAL : public g2o::BaseVertex<3, Vector3d> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	VertexPointBAL() {}
	virtual bool read ( istream& is ) {return false;}
	virtual bool write ( ostream& os ) const {return false;}
	virtual void setToOriginImpl() {}
	
	virtual void oplusImpl ( const double* update ) {
		Vector3d::ConstMapType v ( update );
		_estimate += v;
	}
};

class EdgeObservationBAL : public g2o::BaseBinaryEdge<2, Vector2d, VertexCameraBAL, VertexPointBAL> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	EdgeObservationBAL() {}
	virtual bool read ( istream& is ) {return false;}
	virtual bool write ( ostream& os ) const {return false;}
	
	virtual void computeError() override {
		const VertexCameraBAL* cam = static_cast<const VertexCameraBAL*> ( vertex( 0 ) );
		const VertexPointBAL* point = static_cast<const VertexPointBAL*> ( vertex( 1 ) );
		( *this ) ( cam->estimate().data(), point->estimate().data(), _error.data() );
    }

    template<typename T>
    bool operator() ( const T* camera, const T* point, T* residuals ) const {
        T predictions[2];
        CamProjectionWithDistortion ( camera, point, predictions );
        residuals[0] = predictions[0] - T ( measurement() ( 0 ) );
        residuals[1] = predictions[1] - T ( measurement() ( 1 ) );

        return true;
    }
    
    virtual void linearizeOplus() override {
        const VertexCameraBAL* cam = static_cast<const VertexCameraBAL*> ( vertex ( 0 ) );
        const VertexPointBAL* point = static_cast<const VertexPointBAL*> ( vertex ( 1 ) );
        typedef ceres::internal::AutoDiff<EdgeObservationBAL, double, VertexCameraBAL::Dimension, VertexPointBAL::Dimension> BalAutoDiff;

        Eigen::Matrix<double, Dimension, VertexCameraBAL::Dimension, Eigen::RowMajor> dError_dCamera;
        Eigen::Matrix<double, Dimension, VertexPointBAL::Dimension, Eigen::RowMajor> dError_dPoint;
        double *parameters[] = { const_cast<double*> ( cam->estimate().data() ), const_cast<double*> ( point->estimate().data() ) };
        double *jacobians[] = { dError_dCamera.data(), dError_dPoint.data() };
        double value[Dimension];
        bool diffState = BalAutoDiff::Differentiate ( *this, parameters, Dimension, value, jacobians );

        // copy over the Jacobians (convert row-major -> column-major)
        if ( diffState )
        {
            _jacobianOplusXi = dError_dCamera;
            _jacobianOplusXj = dError_dPoint;
        }
        else
        {
            assert ( 0 && "Error while differentiating" );
            _jacobianOplusXi.setZero();
            _jacobianOplusXj.setZero();
        }
    }
};
#endif 
