#ifndef BALPROBLEM_H
#define BALPROBLEM_H
#include <iostream>
#include <cstdio>
#include <string>
#include <fstream>
#include "g2o_ba_type.h"
#include "tools.h"
using namespace std;


class BALProblem {
public:
    BALProblem(const string& filename);
    ~BALProblem() {
        delete [] point_index_;
        delete [] camera_index_;
        delete [] observations_;
        delete [] parameters_;
    }
    
    void WriteToPLYFile(const string& filename) const;
    void WriteToFile(const string& filename) const;
    
    int camera_block_size() const {return 9;}
    int point_block_size() const {return 3;}
    
    int num_cameras() const {return num_cameras_;}
    int num_points() const {return num_points_;}
    int num_observations() const {return num_observations_;}
    int num_parameters() const {return num_parameters_;}
    
    const int* point_index() const {return point_index_;}
    const int* camera_index() const {return camera_index_;}
    const double* observations() const {return observations_;}
    const double* parameters() const {return parameters_;}
    
    const double* cameras() const {return parameters_;}
    const double* points() const {return parameters_ + camera_block_size() * num_cameras_;}
    
    const double* camera_for_observation(int i) const {
        return cameras() + camera_index_[i] * camera_block_size();
    }
    
    const double* point_for_observation(int i) const {
        return points() + point_index_[i] * point_block_size();
    }
    
    
    double* mutable_cameras() {return parameters_;}
    double* mutable_points() {return parameters_ + camera_block_size() * num_cameras_;}
    double* mutable_camera_for_observation(int i) {
        return mutable_cameras() + camera_index_[i] * camera_block_size();
    }
    double* mutable_point_for_observation(int i) {
        return mutable_points() + point_index_[i] * point_block_size();
    } 
       
    void Normalize();
private:
    void CameraToAngelAxisAndCenter(const double* camera,
                                    double* angle_axis,
                                    double* center)const;
    void AngleAxisAndCenterToCamera(const double* angle_axis,
                                    const double* center,
                                    double* camera)const;
    int num_cameras_;
    int num_points_;
    int num_observations_;
    int num_parameters_;
    string filename_;
    int* point_index_;
    int* camera_index_;
    double* observations_;
    double* parameters_;   
};
#endif