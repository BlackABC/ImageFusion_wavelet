#include <opencv2/core.hpp>

cv::Rect getRoi(cv::Mat infrared, cv::Mat visible);

void round_and_clip_to_0_255(float **im, int N1, int N2);

void wavelet_transform(float **original, int Ni, int Nj, int inverse);