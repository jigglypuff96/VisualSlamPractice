#include "BALProblem.h"
#include "g2o_ba_type.h"
#include "tools.h"
#include <iostream>

#include "g2o/stuff/sampler.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/core/batch_stats.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/optimization_algorithm_dogleg.h"

#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/solvers/pcg/linear_solver_pcg.h"
#include "g2o/types/sba/types_six_dof_expmap.h"

#include "g2o/solvers/structure_only/structure_only_solver.h"

double Median(std::vector<double>* data){
  int n = data->size();
  std::vector<double>::iterator mid_point = data->begin() + n/2;
  std::nth_element(data->begin(),mid_point,data->end());
  return *mid_point;
}

BALProblem::BALProblem(const string& filename) {
	filename_ = filename;
	ifstream fin;
	fin.open(filename_);
	
	if (!fin) {
		cout << "Error opening the file!\n";
		return;
	}
	
	fin >> num_cameras_;
	fin >> num_points_;
	fin >> num_observations_;
	
	cout << "Header: he number of cameras is: " << num_cameras_ << ", ";
	cout << "the number of points is: " << num_points_ << ", ";
	cout << "the number of observations is: " << num_observations_ << ".\n";
	
	point_index_ = new int[num_observations_];
    camera_index_ = new int[num_observations_];
    observations_ = new double[2 * num_observations_];
    
    num_parameters_ = camera_block_size()*num_cameras_ + point_block_size()*num_points_;
    parameters_ = new double[num_parameters_];

	for (int i=0; i<num_observations_; i++) {
		fin >> camera_index_[i];
		fin >> point_index_[i];
		fin >> observations_[2*i];
		fin >> observations_[2*i + 1];	
	}
	
	for (int j=0; j<num_parameters_; j++) {
		fin >> parameters_[j];
	}
	fin.close();
}

// Write the problem to a PLY file for inspection in Meshlab or CloudCompare
void BALProblem::WriteToPLYFile(const std::string& filename)const{
  std::ofstream of(filename.c_str());

  of<< "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex " << num_cameras_ + num_points_
    << '\n' << "property float x"
    << '\n' << "property float y"
    << '\n' << "property float z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue"
    << '\n' << "end_header" << std::endl;

    // Export extrinsic data (i.e. camera centers) as green points.
    double angle_axis[3];
    double center[3];
    for(int i = 0; i < num_cameras(); ++i){
      const double* camera = cameras() + camera_block_size() * i;
      CameraToAngelAxisAndCenter(camera, angle_axis, center);
      of << center[0] << ' ' << center[1] << ' ' << center[2]
         << "0 255 0" << '\n';
    }

    // Export the structure (i.e. 3D Points) as white points.
    const double* points = parameters_ + camera_block_size() * num_cameras_;
    for(int i = 0; i < num_points(); ++i){
      const double* point = points + i * point_block_size();
      for(int j = 0; j < point_block_size(); ++j){
        of << point[j] << ' ';
      }
      of << "255 255 255\n";
    }
    of.close();
}

void BALProblem::WriteToFile(const string& filename) const {
    ofstream fout;
    fout.open(filename);
    fout << num_cameras_ << " " << num_points_ << " " << num_observations_;
    fout << endl;
    for (int i=0; i<num_observations_; i++) {
        fout << camera_index_[i] << " " << point_index_[i] << " " << observations_[2*i] << " " << observations_[2*i+1] << endl;
    }
    
    for (int i=0; i<num_parameters_; i++) {
        fout << parameters_[i];
    }
    fout.close();
}

void BALProblem::CameraToAngelAxisAndCenter(const double* camera, 
                                            double* angle_axis,
                                            double* center) const{
    Eigen::Map<Vector3d> angle_axis_ref(angle_axis);
    angle_axis_ref = Eigen::Map<const Vector3d>(camera);
      
    // c = -R't
    Eigen::Vector3d inverse_rotation = -angle_axis_ref;
    AngleAxisRotatePoint(inverse_rotation.data(),
                         camera + camera_block_size() - 6,
                         center);
    Eigen::Map<Vector3d>(center) *= -1.0;
}

void BALProblem::AngleAxisAndCenterToCamera(const double* angle_axis,
                                            const double* center,
                                            double* camera) const{
    Eigen::Map<const Vector3d> angle_axis_ref(angle_axis);
    Eigen::Map<Vector3d> cam(camera); cam = angle_axis_ref;
    // t = -R * c 
    AngleAxisRotatePoint(angle_axis,center,camera + camera_block_size() - 6);
    Eigen::Map<Vector3d>(camera + camera_block_size() - 6) *= -1.0;
}

void BALProblem::Normalize(){
  // Compute the marginal median of the geometry
  std::vector<double> tmp(num_points_);
  Eigen::Vector3d median;
  double* points = mutable_points();
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < num_points_; ++j){
      tmp[j] = points[3 * j + i];      
    }
    median(i) = Median(&tmp);
  }

  for(int i = 0; i < num_points_; ++i){
    Eigen::Map<Vector3d> point(points + 3 * i, 3);
    tmp[i] = (point - median).lpNorm<1>();
  }

  const double median_absolute_deviation = Median(&tmp);

  // Scale so that the median absolute deviation of the resulting
  // reconstruction is 100

  const double scale = 100.0 / median_absolute_deviation;

  // X = scale * (X - median)
  for(int i = 0; i < num_points_; ++i){
    Eigen::Map<Vector3d> point(points + 3 * i, 3);
    point = scale * (point - median);
  }
  double* cameras = mutable_cameras();
  double angle_axis[3];
  double center[3];
  for(int i = 0; i < num_cameras_ ; ++i){
    double* camera = cameras + camera_block_size() * i;
    CameraToAngelAxisAndCenter(camera, angle_axis, center);
    // center = scale * (center - median)
    Eigen::Map<Vector3d> cen(center);
    cen = scale * (cen - median);
    AngleAxisAndCenterToCamera(angle_axis, center,camera);
  }
}

class G2O_TYPES_SBA_API VertexSE3Expmap : public BaseVertex<6, SE3Quat>{
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VertexSE3Expmap();

    bool read(std::istream& is);

    bool write(std::ostream& os) const;

    virtual void setToOriginImpl() {
        _estimate = SE3Quat();
    }

    virtual void oplusImpl(const number_t* update_)  {
        Eigen::Map<const Vector6> update(update_);
        setEstimate(SE3Quat::exp(update)*estimate());
    }
};

class G2O_TYPES_SBA_API EdgeProjectXYZ2UV : public  BaseBinaryEdge<2, Vector2, VertexSBAPointXYZ, VertexSE3Expmap>{
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    EdgeProjectXYZ2UV();

    bool read(std::istream& is);

    bool write(std::ostream& os) const;

    void computeError()  {
        const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
        const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
        const CameraParameters * cam
                = static_cast<const CameraParameters *>(parameter(0));
        Vector2 obs(_measurement);
        _error = obs-cam->cam_map(v1->estimate().map(v2->estimate()));
    }

    virtual void linearizeOplus();

    CameraParameters * _cam;
};
void EdgeProjectXYZ2UV::linearizeOplus() {
    VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
    SE3Quat T(vj->estimate());
    VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
    Vector3 xyz = vi->estimate();
    Vector3 xyz_trans = T.map(xyz);

    number_t x = xyz_trans[0];
    number_t y = xyz_trans[1];
    number_t z = xyz_trans[2];
    number_t z_2 = z*z;

    const CameraParameters * cam = static_cast<const CameraParameters *>(parameter(0));

    Matrix<number_t,2,3,Eigen::ColMajor> dudq;
    auto fx = cam->focal_length;
    auto fy = cam->focal_length;
    dudq(0,0) = fx/z;
    dudq(0,1) = 0;
    dudq(0,2) = -x*fx/z_2;

    dudq(1,0) = 0;
    dudq(1,1) = fy/z;
    dudq(1,2) = -fy*y/z_2;

    _jacobianOplusXi =  -1 * dudq * T.rotation().toRotationMatrix();


    _jacobianOplusXj(0,0) = fx/z;
    _jacobianOplusXj(0,1) = 0;
    _jacobianOplusXj(0,2) = fx*x/z_2;
    _jacobianOplusXj(0,3) = -fx*x*y/z_2;
    _jacobianOplusXj(0,4) = fx*(1+(x*x/z_2));
    _jacobianOplusXj(0,5) = -fx*y/z;

    _jacobianOplusXj(1,0) = 0;
    _jacobianOplusXj(1,1) = fy/z;
    _jacobianOplusXj(1,2) = -fy*y/z_2;
    _jacobianOplusXj(1,3) = -fy*(1+y*y/z_2);
    _jacobianOplusXj(1,4) = fy*x*y/z_2;
    _jacobianOplusXj(1,5) = fy*x/z;

    _jacobianOplusXj = -_jacobianOplusXj;
}