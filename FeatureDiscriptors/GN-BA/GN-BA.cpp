//
// Created by xiang on 12/21/17.
//

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "sophus/se3.hpp"

using namespace std;

typedef vector<Vector3d, Eigen::aligned_allocator<Vector3d>> VecVector3d;
typedef vector<Vector2d, Eigen::aligned_allocator<Vector3d>> VecVector2d;
typedef Matrix<double, 6, 1> Vector6d;

string p3d_file = "../p3d.txt";
string p2d_file = "../p2d.txt";

int main(int argc, char **argv) {

    VecVector2d p2d;
    VecVector3d p3d;
    Matrix3d K;
    double fx = 520.9, fy = 521.0, cx = 325.1, cy = 249.7;
    K << fx, 0, cx, 0, fy, cy, 0, 0, 1;

    // load points in to p3d and p2d 
    // START YOUR CODE HERE
    ifstream file3d(p3d_file);
    ifstream file2d(p2d_file);

    double pd1a, pd2a, pd3a, pd1b, pd2b;
    while (!file3d.eof()){
        file3d >> pd1a >> pd2a >> pd3a;
        Vector3d temp(pd1a,pd2a,pd3a);
        p3d.push_back(temp);
    }
    while (!file2d.eof()){
        file2d >> pd1b >> pd2b;
        Vector2d temp(pd1b,pd2b);
        p2d.push_back(temp);
    }

    // END YOUR CODE HERE
    assert(p3d.size() == p2d.size());

    int iterations = 100;
    double cost = 0, lastCost = 0;
    int nPoints = p3d.size();
    cout << "points: " << nPoints << endl;

    Sophus::SE3d T_esti; // estimated pose

    for (int iter = 0; iter < iterations; iter++) {

        Matrix<double, 6, 6> H = Matrix<double, 6, 6>::Zero();
        Vector6d b = Vector6d::Zero();

        cost = 0;
        // compute cost
        for (int i = 0; i < nPoints; i++) {
            // compute cost for p3d[I] and p2d[I]
            // START YOUR CODE HERE 
            Vector3d P_prime = T_esti*p3d[i];
            Vector3d U = K*P_prime;
            //the coordinates of the projection point
            Vector2d e = p2d[i] - Vector2d(U(0)/U(2),U(1)/U(2)); // projection
            cost += e.squaredNorm()/2;
	    // END YOUR CODE HERE

	    // compute jacobian
            Matrix<double, 2, 6> J;
            // START YOUR CODE HERE 
            double x_prime = P_prime(0);
            double y_prime = P_prime(1);
            double z_prime = P_prime(2);
            double x_prime_sqr = x_prime*x_prime;
            double y_prime_sqr = y_prime*y_prime;
            double z_prime_sqr = z_prime*z_prime;
            // u = fx*x_prime/z_prime + cx
            // v = fy*y_prime/z_prime + cy;

            J(0,0) = fx/z_prime;
            J(0,1) = 0;
            J(0,2) = -fx*x_prime/z_prime_sqr;
            J(0,3) = -fx*x_prime*y_prime/z_prime_sqr;
            J(0,4) = fx + fx*x_prime_sqr/z_prime_sqr;
            J(0,5) = -fx*y_prime/z_prime;

            J(1,0) = 0;
            J(1,1) = fy/z_prime;
            J(1,2) = -fy/z_prime_sqr;
            J(1,3) = -fy-fy*y_prime_sqr/z_prime_sqr;
            J(1,4) = fy*x_prime*y_prime/z_prime_sqr;
            J(1,5) = fy*x_prime/z_prime;

            J = -J;


	    // END YOUR CODE HERE

            H += J.transpose() * J;
            b += -J.transpose() * e;
        }

	// solve dx 
        Vector6d dx;

        // START YOUR CODE HERE 
        dx = H.lu().solve(b);
        // END YOUR CODE HERE

        if (isnan(dx[0])) {
            cout << "result is nan!" << endl;
            break;
        }

        if (iter > 0 && cost >= lastCost) {
            // cost increase, update is not good
            cout << "cost: " << cost << ", last cost: " << lastCost << endl;
            break;
        }

        // update your estimation
        // START YOUR CODE HERE 
        T_esti = Sophus::SE3d::exp(dx)*T_esti;
        // END YOUR CODE HERE
        
        lastCost = cost;

        cout << "iteration " << iter << " cost=" << cout.precision(12) << cost << endl;
    }

    cout << "estimated pose: \n" << T_esti.matrix() << endl;
    return 0;
}
