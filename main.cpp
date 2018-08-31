#include <opencv2\opencv.hpp>
#include <opencv2\xfeatures2d\nonfree.hpp>
#include <opencv2\features2d\features2d.hpp>
#include <opencv2\highgui\highgui.hpp>

#include <iostream>
#include <vector>
#include "alloc.h"
#include "util.h"

using namespace cv;
using namespace cv::xfeatures2d;
using namespace std;

int main()
{
	////test RGB to YUV
	//cv::Mat infrared = cv::imread("Tank_IR.tif", 0);
	//cv::Mat visible_origin = cv::imread("Tank_Vis.tif");
	cv::Mat infrared = cv::imread("output_up_ir.jpg", 0);
	cv::Mat visible_origin = cv::imread("output_up_vis.jpg");

	cv::Rect roi = getRoi(infrared, visible_origin);
	cv::Mat visible = visible_origin(roi);

	cv::imshow("infrared", infrared);
	cv::imshow("visible", visible);

	cv::Mat yuv;
	cv::cvtColor(visible, yuv, CV_BGR2YUV);

	std::vector<cv::Mat> yuvArr;
	cv::split(yuv, yuvArr);
	//yuvArr[0] = infrared;
	//cv::imshow("h1", yuvArr[0]);
	//cv::imshow("h2", yuvArr[1]);
	//cv::imshow("h3", yuvArr[2]);

	double start = cv::getTickCount();

	float **inf;
	float **vis;
	inf = allocate_2d_float(roi.height, roi.width, 0);
	vis = allocate_2d_float(roi.height, roi.width, 0);

	for (int i = 0; i < roi.height; i++)
		for (int j = 0; j < roi.width; j++)
		{
			inf[i][j] = infrared.at<uchar>(i, j);
			vis[i][j] = yuvArr[0].at<uchar>(i, j);
		}
	
	wavelet_transform(inf, roi.height, roi.width, 1);
	wavelet_transform(vis, roi.height, roi.width, 1);

	//cv::Mat aaa(roi.height, roi.width, CV_8U);
	//cv::Mat bbb(roi.height, roi.width, CV_8U);
	//for (int i = 0; i < roi.height; i++)
	//	for (int j = 0; j < roi.width; j++)
	//	{
	//		aaa.at<uchar>(i, j) = (int)inf[i][j];
	//		bbb.at<uchar>(i, j) = (int)vis[i][j];
	//	}
	//cv::imshow("aaa", aaa);
	//cv::imshow("bbb", bbb);

	for (int i = 0; i < roi.height; i++)
		for (int j = 0; j < roi.width; j++)
		{
			//yuvArr[0].at<uchar>(i, j) = (yuvArr[0].at<uchar>(i, j) + infrared.at<uchar>(i, j)) / 2;
			////线性组合
			vis[i][j] = inf[i][j] * 0.2 + vis[i][j] * 0.8;
			////取极大值
			//if (inf[i][j]>vis[i][j])
			//{
			//	vis[i][j] = inf[i][j];
			//}
		}
	
	wavelet_transform(vis, roi.height, roi.width, 0);

	round_and_clip_to_0_255(vis, roi.height, roi.width);

	for (int i = 0; i < roi.height; i++)
		for (int j = 0; j < roi.width; j++)
		{
			yuvArr[0].at<uchar>(i, j) = vis[i][j];
			//取大情况下会出现像素值大于255的情况
		}

	double end = cv::getTickCount();
	double t = (end - start) / cv::getTickFrequency();
	std::cout << t << std::endl;

	cv::Mat fusionImage;
	cv::merge(yuvArr, fusionImage);

	cv::cvtColor(fusionImage, fusionImage, CV_YUV2BGR);
	cv::imshow("fusionImage", fusionImage);

	cv::imwrite("fusion.jpg", fusionImage);

	cv::waitKey(0);
}