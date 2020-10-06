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

typedef Eigen::Map<const Vector3d> ConstVecRef;
typedef Eigen::Map<Vector3d> VecRef;
typedef g2o::BlockSolver<g2o::BlockSolverTraits<9,3> > BalBlockSolver;

void BuildProblem(const BALProblem* bal_problem, g2o::SparseOptimizer* optimizer) {
    const int num_cameras = bal_problem->num_cameras();
    const int num_points = bal_problem->num_points();
    const int num_observations = bal_problem->num_observations();
    const int camera_block_size = bal_problem->camera_block_size();
    const int point_block_size = bal_problem->point_block_size();
    
    const double* cameras = bal_problem->cameras();
    for (int i=0; i<num_cameras; i++) {
        
        Eigen::Map<const Eigen::Matrix<double, 9, 1>> tempcam(cameras + i*camera_block_size);
        VertexCameraBAL* pcam = new VertexCameraBAL();
        pcam->setId(i);
        pcam->setEstimate(tempcam);
        
        optimizer->addVertex(pcam);
    }
    
    const double* points = bal_problem->points();
    for (int i=0; i<num_points; i++) {
        ConstVecRef temppoi(points + i*point_block_size);
        VertexPointBAL* ppoint = new VertexPointBAL();
        ppoint->setId(num_cameras + i);
        ppoint->setEstimate(temppoi);
        ppoint->setMarginalized(true);
        optimizer->addVertex(ppoint);
    }
    
    const double* observations = bal_problem->observations();
    for (int i=0; i<num_observations; i++) {
        const int cameraid = bal_problem->camera_index()[i];
        const int pointid = bal_problem->point_index()[i] + num_cameras;
        
        EdgeObservationBAL* pedge = new EdgeObservationBAL();
        
        g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
        rk->setDelta(1.0);
        pedge->setRobustKernel(rk);
        
        pedge->setVertex(0, dynamic_cast<VertexCameraBAL*>(optimizer->vertex(cameraid)));
        pedge->setVertex(1, dynamic_cast<VertexPointBAL*>(optimizer->vertex(pointid)));
        
        pedge->setInformation(Eigen::Matrix2d::Identity());
        pedge->setMeasurement(Eigen::Vector2d(observations[2*i],observations[2*i+1]));
        
        optimizer->addEdge(pedge);        
    }   
    
}

void WriteToBALProblem(BALProblem* bal_problem, g2o::SparseOptimizer* optimizer) {
    const int num_cameras = bal_problem->num_cameras();
    const int num_points = bal_problem->num_points();
    const int num_observations = bal_problem->num_observations();
    const int camera_block_size = bal_problem->camera_block_size();
    const int point_block_size = bal_problem->point_block_size();
    
    double* cameras = bal_problem->mutable_cameras();
    for (int i=0; i<num_cameras; i++) {
        VertexCameraBAL* pcam = dynamic_cast<VertexCameraBAL*> (optimizer->vertex(i));
        Vector9d newcam = pcam->estimate();
        memcpy(cameras + i*camera_block_size, newcam.data(), sizeof(double)*camera_block_size);
    }
    
    double* points = bal_problem->mutable_points();
    for (int i=0; i<num_points; i++) {
        VertexPointBAL* ppoint = dynamic_cast<VertexPointBAL*> (optimizer->vertex(i+num_cameras));
        Vector3d newpoint = ppoint->estimate();
        memcpy(points + i*point_block_size, newpoint.data(), sizeof(double)*point_block_size);
    }
}

// 最新版的g2o的写法在注释中
void setSolverOptions(g2o::SparseOptimizer* optimizer) {
    
    BalBlockSolver* solver_ptr;
    g2o::LinearSolver<BalBlockSolver::PoseMatrixType>* linearSolver = new g2o::LinearSolverCholmod<BalBlockSolver::PoseMatrixType>();
    dynamic_cast<g2o::LinearSolverCholmod<BalBlockSolver::PoseMatrixType>*>(linearSolver)->setBlockOrdering(true);
    
    solver_ptr = new BalBlockSolver(linearSolver);
    
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
	
	/*
    std::unique+ptr<g2o::LinearSolver<BalBlockSolver::PoseMatrixType>> linearSolver(new g2o::LinearSolverCholmod<BalBlockSolver::PoseMatrixType>());
    dynamic_cast<g2o::LinearSolverCholmod<BalBlockSolver::PoseMatrixType>*>(linearSolver.get())->setBlockOrdering(true);
    
    std::unique_ptr<BalBlockSolver> solver_ptr( new BalBlockSolver(std::move(linearSolver)) );
    
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(std::move(solver_ptr));
    */ 
    optimizer->setAlgorithm(solver);
}

void SolveBALProblem(const string& filename) {
    BALProblem bal_problem(filename);
    
    // show some information here ...
    std::cout << "bal problem file loaded..." << std::endl;
    std::cout << "bal problem have " << bal_problem.num_cameras() << " cameras and " << bal_problem.num_points() << " points. " << std::endl;
    std::cout << "Forming " << bal_problem.num_observations() << " observatoins. " << std::endl;
	
    std::cout << "beginning problem..." << std::endl;
    
    bal_problem.Normalize();
    std::cout << "Normalization complete..." << std::endl;
    
    bal_problem.WriteToPLYFile("initial.ply");
    
    g2o::SparseOptimizer optimizer;
    setSolverOptions(&optimizer);
    
    BuildProblem(&bal_problem, &optimizer);
    
    std::cout << "begin optimizaiton .."<< std::endl;
    // perform the optimizaiton 
    optimizer.initializeOptimization();
    optimizer.setVerbose(true);
    optimizer.optimize(20);

    std::cout << "optimization complete.. "<< std::endl;
    // write the optimized data into BALProblem class
    WriteToBALProblem(&bal_problem, &optimizer);
    bal_problem.WriteToPLYFile("final.ply");
    bal_problem.WriteToFile("final_data.txt");
    
}

int main(int argc, char** argv) {    
    if (argc != 2) {
        cout << "usage: bundle_adjustment <path for datasest>";
        return 1;
    }  
    string filename = argv[1];
    SolveBALProblem(filename);
    return 0;
}
