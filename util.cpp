#include "util.h"

cv::Rect getRoi(cv::Mat infrared, cv::Mat visible)
{
	int infrared_rows = infrared.rows;
	int infrared_cols = infrared.cols;
	int visible_rows = visible.rows;
	int visible_cols = visible.cols;

	int x = (visible_cols - infrared_cols) / 2;
	int y = (visible_rows - infrared_rows) / 2;

	cv::Rect roi(x, y, infrared_cols, infrared_rows);

	return roi;
}