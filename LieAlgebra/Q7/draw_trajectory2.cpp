#include <sophus/se3.hpp>
#include <unistd.h>
#include <string>
#include <iostream>
#include <fstream>

// need pangolin for plotting trajectory
#include <pangolin/pangolin.h>

using namespace std;

// path to trajectory file
//string trajectory_file = "../trajectory.txt";
string groundtruth_file = "../groundtruth.txt";
string estimated_file = "../estimated.txt";

// function for plotting trajectory, don't edit this code
// start point is red and end point is blue
void DrawTrajectory(vector<Sophus::SE3d, Eigen::aligned_allocator<Sophus::SE3d>>);

typedef vector<Sophus::SE3d, Eigen::aligned_allocator<Sophus::SE3d>> myTrajectory;
myTrajectory ReadTrajectory (const string& path);
int main(int argc, char **argv) {
    // bonus question starts here:
    double rmse = 0;

    myTrajectory gt = ReadTrajectory(groundtruth_file);
    myTrajectory est = ReadTrajectory(estimated_file);
    int upperlimit = gt.size();
    for (int i = 0; i < upperlimit;  i++){
        rmse += ((gt[i].inverse()*est[i]).log()).squaredNorm();
    }
    rmse = rmse/double (gt.size());
    rmse = sqrt(rmse);
    cout <<"RMSE = " << rmse << endl;


    // draw trajectory in pangolin
    vector<Sophus::SE3d, Eigen::aligned_allocator<Sophus::SE3d>> twoPoses;
    for (int i = 0; i < upperlimit;  i++){
        twoPoses.push_back(gt[i]);
    }
    for (int i = 0; i < upperlimit;  i++){
        twoPoses.push_back(est[i]);
    }
    DrawTrajectory(twoPoses);
    return 0;
}


myTrajectory ReadTrajectory(const string& path){
    myTrajectory poses;
    double t,tx,ty,tz,qx,qy,qz,qw;
    ifstream myfile (path);
    if (myfile.is_open()) {
        while (!myfile.eof()) {
            myfile >> t >> tx >> ty >> tz >> qx >> qy >> qz >> qw;
            Eigen::Quaterniond q = Eigen::Quaterniond(qw, qx, qy, qz).normalized();
            Eigen::Vector3d t(tx, ty, tz);
            Sophus::SE3d SE3_qt(q, t);
            poses.push_back(SE3_qt);
        }
        myfile.close();
    }
    return poses;
}
/*******************************************************************************************/
void DrawTrajectory(vector<Sophus::SE3d, Eigen::aligned_allocator<Sophus::SE3d>> poses) {
    if (poses.empty()) {
        cerr << "Trajectory is empty!" << endl;
        return;
    }

    // create pangolin window and plot the trajectory
    pangolin::CreateWindowAndBind("Trajectory Viewer", 1024, 768);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    pangolin::OpenGlRenderState s_cam(
            pangolin::ProjectionMatrix(1024, 768, 500, 500, 512, 389, 0.1, 1000),
            pangolin::ModelViewLookAt(0, -0.1, -1.8, 0, 0, 0, 0.0, -1.0, 0.0)
    );

    pangolin::View &d_cam = pangolin::CreateDisplay()
            .SetBounds(0.0, 1.0, pangolin::Attach::Pix(175), 1.0, -1024.0f / 768.0f)
            .SetHandler(new pangolin::Handler3D(s_cam));


    while (pangolin::ShouldQuit() == false) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        d_cam.Activate(s_cam);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

        glLineWidth(2);
        for (size_t i = 0; i < poses.size() - 1; i++) {
            glColor3f(1 - (float) i / poses.size(), 0.0f, (float) i / poses.size());
            glBegin(GL_LINES);
            auto p1 = poses[i], p2 = poses[i + 1];
            glVertex3d(p1.translation()[0], p1.translation()[1], p1.translation()[2]);
            glVertex3d(p2.translation()[0], p2.translation()[1], p2.translation()[2]);
            glEnd();
        }
        pangolin::FinishFrame();
        usleep(5000);   // sleep 5 ms
    }

}