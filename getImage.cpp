#include <opencv.hpp>

using namespace cv;
using namespace std;

int getImage_main()
{
	VideoCapture cap("visible.avi");
	Mat image;
	cap.read(image);
	imwrite("visible.jpg", image);
	return 0;
}