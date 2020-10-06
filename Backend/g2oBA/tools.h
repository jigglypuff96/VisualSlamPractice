#ifndef TOOLS_H
#define TOOLS_H
#include <algorithm>
#include <cmath>
#include <limits>

template<typename T> 
inline T DotProduct(const T x[3], const T y[3]) {
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

template<typename T>
inline void CrossProduct(const T x[3], const T y[3], T result[3]){
  result[0] = x[1] * y[2] - x[2] * y[1];
  result[1] = x[2] * y[0] - x[0] * y[2];
  result[2] = x[0] * y[1] - x[1] * y[0];
}

template<typename T>
inline void AngleAxisRotatePoint(const T angle_axis[3], const T pt[3], T result[3]) {
  const T theta2 = DotProduct(angle_axis, angle_axis);
  if (theta2 > T(std::numeric_limits<double>::epsilon())) {
  
    const T theta = sqrt(theta2);
    const T costheta = cos(theta);
    const T sintheta = sin(theta);
    const T theta_inverse = 1.0 / theta;

    const T w[3] = { angle_axis[0] * theta_inverse,
                     angle_axis[1] * theta_inverse,
                     angle_axis[2] * theta_inverse };

    T w_cross_pt[3];
    CrossProduct(w, pt, w_cross_pt);                          


    const T tmp = DotProduct(w, pt) * (T(1.0) - costheta);

    result[0] = pt[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp;
    result[1] = pt[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp;
    result[2] = pt[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp;
  } else {
    T w_cross_pt[3];
    CrossProduct(angle_axis, pt, w_cross_pt); 

    result[0] = pt[0] + w_cross_pt[0];
    result[1] = pt[1] + w_cross_pt[1];
    result[2] = pt[2] + w_cross_pt[2];
  }
}

template<typename T>
inline bool CamProjectionWithDistortion(const T* camera, const T* point, T* predictions){
    T p[3];
    AngleAxisRotatePoint(camera, point, p);
    // camera[3,4,5] are the translation
    p[0] += camera[3]; p[1] += camera[4]; p[2] += camera[5];

    // Compute the center fo distortion
    T xp = -p[0]/p[2];
    T yp = -p[1]/p[2];

    // Apply second and fourth order radial distortion
    const T& l1 = camera[7];
    const T& l2 = camera[8];

    T r2 = xp*xp + yp*yp;
    T distortion = T(1.0) + r2 * (l1 + l2 * r2);

    const T& focal = camera[6];
    predictions[0] = focal * distortion * xp;
    predictions[1] = focal * distortion * yp;

    return true;
}
#endif
