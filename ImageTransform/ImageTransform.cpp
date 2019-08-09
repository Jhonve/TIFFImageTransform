#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <tiffio.h>

#define PI 3.141592653

using namespace std;
using namespace cv;

void rotate(Vec2d* ray_index, int position_x, int position_y, int max_length, float angle)
{
	for (int i = 0; i < max_length; i++)
	{
		ray_index[i][0] = round(((position_x + i - position_x)*cos(angle)) - ((position_y - position_y) * sin(angle)) + position_x);
		ray_index[i][1] = round(((position_x + i - position_x)*sin(angle)) + ((position_y - position_y) * cos(angle)) + position_y);
	}
}

bool ImageTransform(const char* input_file, int position_x, int position_y, int shape_w, int shape_h) {
	Mat input_image = imread(input_file, IMREAD_UNCHANGED);
	// Mat input_image = imread(input_file);
	// Mat input_image = imread(input_file, IMREAD_GRAYSCALE);
	
	int input_height = input_image.rows;
	int input_width = input_image.cols;
	// cout << input_image.dims << " "<< input_image.channels() << " "<< input_height << " " << input_width << endl;

	int max_length = ceil(sqrt(pow(max(input_width - position_x, position_x), 2) + \
		pow(max(input_height - position_y, position_y), 2)));
	// cout << "max_length: " << max_length << endl;

	Mat output_image = Mat(max_length, 360, CV_8UC1, Scalar(0));
	// imwrite("OutputTempShape.tif", output_image);

	Vec2d* ray_index;
	ray_index = (Vec2d*)malloc(sizeof(Vec2d)*max_length);

	/*for (int i = 135; i <= 135; i++)
	{
		rotate(ray_index, position_x, position_y, max_length, -1 * float(i) / 180.f*PI);
	}

	int count = 0;
	for (int i = 0; i < max_length; i++)
	{
		int draw_pos_x = ray_index[i][0];
		int draw_pos_y = ray_index[i][1];
		if (draw_pos_x < input_width && draw_pos_y < input_height && draw_pos_x > -1 && draw_pos_y > -1)
		{
			Point point(draw_pos_x, draw_pos_y);
			circle(input_image, point, 1, Scalar(255));
			count++;
		}
	}*/

	for (int i = 1; i <= 360; i++)
	{
		rotate(ray_index, position_x, position_y, max_length, -1 * float(i) / 180.f*PI);

		for (int j = 0; j < max_length; j++)
		{
			int select_pos_x = ray_index[j][0];
			int select_pos_y = ray_index[j][1];
			
			if (select_pos_x < input_width && select_pos_y < input_height && select_pos_x > -1 && select_pos_y > -1)
			{
				int value = input_image.at<uchar>(select_pos_x, select_pos_y);
				// cout << "x: " << select_pos_x << "  y:" << select_pos_y << "  value: " << value << endl;

				output_image.at<uchar>(max_length - j - 1, i-1) = value;
			}
		}
	}

	Size output_size = Size(shape_w, shape_h);
	resize(output_image, output_image, output_size, cv::INTER_AREA);
	imwrite("OutputImage.tif", output_image);

	return true;
}

bool TiffTransform(const char* input_file, int position_x, int position_y, int shape_w, int shape_h) {
	TIFF* input_tif = TIFFOpen(input_file, "r");

	if (input_tif)
	{
		cout << "Read success!!" << endl;
	}
	else
	{
		cout << "Load input failed. Please input right .tif" << endl;
		return false;
	}

	int dir_count = 0;
	do
	{
		dir_count++;
	} while (TIFFReadDirectory(input_tif));

	uint32 input_width, input_height = 0;
	int input_channels = 0;
	int input_bits = 0;
	uint16 input_data_type = 0;
	int input_photo_metric = 0;
	float input_x_resolution, input_y_resolution = 0.f;
	uint16 input_planar_config = 0;
	uint16 input_res_unit = 0;
	
	TIFFGetField(input_tif, TIFFTAG_IMAGEWIDTH, &input_width);
	TIFFGetField(input_tif, TIFFTAG_IMAGELENGTH, &input_height);
	TIFFGetField(input_tif, TIFFTAG_SAMPLESPERPIXEL, &input_channels);
	TIFFGetField(input_tif, TIFFTAG_BITSPERSAMPLE, &input_bits);
	TIFFGetField(input_tif, TIFFTAG_DATATYPE, &input_data_type);
	TIFFGetField(input_tif, TIFFTAG_PHOTOMETRIC, &input_photo_metric);
	TIFFGetField(input_tif, TIFFTAG_XRESOLUTION, &input_x_resolution);
	TIFFGetField(input_tif, TIFFTAG_YRESOLUTION, &input_y_resolution);
	TIFFGetField(input_tif, TIFFTAG_RESOLUTIONUNIT, &input_res_unit);
	TIFFGetField(input_tif, TIFFTAG_PLANARCONFIG, &input_planar_config);
	
	int tag_c = TIFFGetTagListCount(input_tif);

	cout << "Tags show: (number of tags:" << tag_c << ")" << endl;
	cout << "numOffiles * width * height * channels * bitsPerChannel" << endl;
	cout << dir_count << " * " << input_width << " * " << input_height << " * "  << input_channels << " * " << input_bits << endl;
	cout << "data type: " << input_data_type << endl;
	cout << "photo metric: " << input_photo_metric << endl; // 1: gary image, 0 means the darkest.
	cout << "xresolution: " << input_x_resolution << " , y resolution " << input_y_resolution << ", resolution unit " << input_res_unit << endl;
	cout << "planar config: " << input_planar_config << endl;

	TIFF* output_tif = TIFFOpen("OutputTif.tif", "w");
	TIFFSetField(output_tif, TIFFTAG_IMAGEWIDTH, input_width);
	TIFFSetField(output_tif, TIFFTAG_IMAGELENGTH, input_height);
	TIFFSetField(output_tif, TIFFTAG_SAMPLESPERPIXEL, input_channels);
	TIFFSetField(output_tif, TIFFTAG_BITSPERSAMPLE, input_bits);
	TIFFSetField(output_tif, TIFFTAG_PHOTOMETRIC, input_photo_metric);
	TIFFSetField(output_tif, TIFFTAG_XRESOLUTION, input_x_resolution);
	TIFFSetField(output_tif, TIFFTAG_YRESOLUTION, input_y_resolution);
	TIFFSetField(output_tif, TIFFTAG_RESOLUTIONUNIT, input_res_unit);
	TIFFSetField(output_tif, TIFFTAG_PLANARCONFIG, input_planar_config);
	TIFFSetField(output_tif, TIFFTAG_DATATYPE, input_data_type);

	// Scanline-based Image I/O
	tdata_t input_buf;
	uint32 row;

	uint32** middle_buffer;
	middle_buffer = new uint32*[input_height];
	for (int i = 0; i < input_height; i++)
	{
		middle_buffer[i] = new uint32[input_width];
	}

	/*uint32** middle_buffer = (uint32**)malloc(sizeof(uint32*)*input_height);
	for (int i = 0; i < input_height; i++)
	{
		middle_buffer[i] = (uint32*)malloc(sizeof(uint32)*input_width);
	}*/

	input_buf = _TIFFmalloc(TIFFStripSize(input_tif));
	if (input_planar_config == PLANARCONFIG_CONTIG)
	{
		cout << "write correct" << endl;
		for (row = 0; row < input_height; row++)
		{
			TIFFReadScanline(input_tif, input_buf, row);
			// TIFFWriteScanline(output_tif, input_buf, row);

			middle_buffer[row] = (uint32*)input_buf;
			TIFFWriteScanline(output_tif, middle_buffer[row], row);

			//for (int i = 0; i < 1024; i++)
			//{
			//	/*if (middle_buffer[row][i] > 0 && middle_buffer[row][i] < 2147483648)
			//	{
			//		cout << middle_buffer[row][i] << endl;
			//	}*/
			//	cout << middle_buffer[row][i] << endl;
			//}
		}
	}
	else if (input_planar_config == PLANARCONFIG_SEPARATE)
	{
		uint16 s, input_nsamples;

		TIFFGetField(input_tif, TIFFTAG_SAMPLESPERPIXEL, &input_nsamples);
		TIFFSetField(output_tif, TIFFTAG_SAMPLESPERPIXEL, input_nsamples);

		for (s = 0; s < input_nsamples; s++)
		{
			for (row = 0; row < input_height; row++)
			{
				TIFFReadScanline(input_tif, input_buf, row, s);
				TIFFWriteScanline(output_tif, input_buf, row, s);
			}
		}
	}

	/*for (int i = 250; i >= 0; i--)
	{
		cout << i << " " << middle_buffer[i][i] << endl;
	}*/

	_TIFFfree(input_buf);
	TIFFClose(input_tif);
	TIFFClose(output_tif);
	free(middle_buffer);

	return true;
}

struct myVec2d
{
	int x;
	int y;
};

void rotatemy(myVec2d ray_index[], int position_x, int position_y, int max_length, float angle)
{
	for (int i = 0; i < max_length; i++)
	{
		ray_index[i].x = round(((position_x + i - position_x)*cos(angle)) - ((position_y - position_y) * sin(angle)) + position_x - 1);
		ray_index[i].y = round(((position_x + i - position_x)*sin(angle)) + ((position_y - position_y) * cos(angle)) + position_y - 1);
	}
}

bool TiffTransformReal(const char* input_file, int position_x, int position_y, int shape_w, int shape_h) {
	TIFF* input_tif = TIFFOpen(input_file, "r");

	if (input_tif)
	{
		cout << "Read success!!" << endl;
	}
	else
	{
		cout << "Load input failed. Please input right .tif" << endl;
		return false;
	}

	int dir_count = 0;
	do
	{
		dir_count++;
	} while (TIFFReadDirectory(input_tif));

	uint32 input_width, input_height = 0;
	int input_channels = 0;
	int input_bits = 0;
	uint16 input_data_type = 0;
	int input_photo_metric = 0;
	float input_x_resolution, input_y_resolution = 0.f;
	uint16 input_planar_config = 0;
	uint16 input_res_unit = 0;
	uint16 input_nsamples = 0;

	TIFFGetField(input_tif, TIFFTAG_IMAGEWIDTH, &input_width);
	TIFFGetField(input_tif, TIFFTAG_IMAGELENGTH, &input_height);
	TIFFGetField(input_tif, TIFFTAG_SAMPLESPERPIXEL, &input_channels);
	TIFFGetField(input_tif, TIFFTAG_BITSPERSAMPLE, &input_bits);
	TIFFGetField(input_tif, TIFFTAG_DATATYPE, &input_data_type);
	TIFFGetField(input_tif, TIFFTAG_PHOTOMETRIC, &input_photo_metric);
	TIFFGetField(input_tif, TIFFTAG_XRESOLUTION, &input_x_resolution);
	TIFFGetField(input_tif, TIFFTAG_YRESOLUTION, &input_y_resolution);
	TIFFGetField(input_tif, TIFFTAG_RESOLUTIONUNIT, &input_res_unit);
	TIFFGetField(input_tif, TIFFTAG_PLANARCONFIG, &input_planar_config);

	int tag_c = TIFFGetTagListCount(input_tif);

	cout << "Tags show: (number of tags:" << tag_c << ")" << endl;
	cout << "numOffiles * width * height * channels * bitsPerChannel" << endl;
	cout << dir_count << " * " << input_width << " * " << input_height << " * " << input_channels << " * " << input_bits << endl;
	cout << "data type: " << input_data_type << " , photo metric: " << input_photo_metric << endl; // 1: gary image, 0 means the darkest.
	cout << "xresolution: " << input_x_resolution << " , y resolution " << input_y_resolution << ", resolution unit " << input_res_unit << endl;
	cout << "planar config: " << input_planar_config << endl;
	cout << endl;

	int max_length = ceil(sqrt(pow(max(int(input_width) - position_x, position_x), 2) + \
		pow(max(int(input_height) - position_y, position_y), 2)));
	cout << "Max length of output Tif is: " << max_length << endl;

	// Scanline-based Image I/O
	tdata_t input_buf;
	uint32 row;

	uint32** middle_buffer;
	middle_buffer = new uint32*[input_height];
	for (int i = 0; i < input_height; i++)
	{
		middle_buffer[i] = new uint32[input_width];
	}

	/*uint32** middle_buffer = (uint32**)malloc(sizeof(uint32*)*input_height);
	for (int i = 0; i < input_height; i++)
	{
		middle_buffer[i] = (uint32*)malloc(sizeof(uint32)*input_width);
	}*/

	uint32* temp_buffer;
	input_buf = _TIFFmalloc(TIFFStripSize(input_tif));
	if (input_planar_config == PLANARCONFIG_CONTIG)
	{
		cout << "Read Correct" << endl;
		for (row = 0; row < input_height; row++)
		{
			TIFFReadScanline(input_tif, input_buf, row);
			temp_buffer = (uint32*)input_buf;
			
			for (int col = 0; col < input_width; col++)
			{
				middle_buffer[row][col] = temp_buffer[col];
			}
		}
	}
	else if (input_planar_config == PLANARCONFIG_SEPARATE)
	{
		uint16 s;

		TIFFGetField(input_tif, TIFFTAG_SAMPLESPERPIXEL, &input_nsamples);

		for (s = 0; s < input_nsamples; s++)
		{
			for (row = 0; row < input_height; row++)
			{
				TIFFReadScanline(input_tif, input_buf, row, s);
				temp_buffer = (uint32*)input_buf;

				for (int col = 0; col < input_width; col++)
				{
					middle_buffer[row][col] = temp_buffer[col];
				}
			}
		}
	}

	uint32** output_buffer;
	output_buffer = new uint32*[max_length];
	for (int i = 0; i < max_length; i++)
	{
		output_buffer[i] = new uint32[360];
	}
	
	/*uint32** output_buffer = (uint32**)malloc(sizeof(uint32*)*max_length);
	for (int i = 0; i < max_length; i++)
	{
		output_buffer[i] = (uint32*)malloc(sizeof(uint32)*360);
	}*/

	myVec2d* ray_index;
	ray_index = (myVec2d*)malloc(sizeof(myVec2d)*max_length);

	for (int i = 1; i <= 360; i++)
	{
		rotatemy(ray_index, position_x, position_y, max_length, -1 * float(i) / 180.f*PI);
		for (int j = 0; j < max_length; j++)
		{
			int select_pos_x = ray_index[j].y;
			int select_pos_y = ray_index[j].x;

			if (select_pos_x < input_height && select_pos_y < input_width && select_pos_x > -1 && select_pos_y > -1)
			{
				output_buffer[max_length - j - 1][i - 1] = middle_buffer[select_pos_x][select_pos_y];
				
				//if (output_buffer[max_length - j - 1][i - 1] > 0 && output_buffer[max_length - j - 1][i - 1] < 2147483648)
				//{
				//	// cout << output_buffer[max_length - j - 1][i - 1] << endl;
				//	cout << select_pos_x << " " << select_pos_y << " " << output_buffer[max_length - j - 1][i - 1] << " " << middle_buffer[select_pos_x][select_pos_y] << endl;
				//}
			}
			else
			{
				output_buffer[max_length - j - 1][i - 1] = 0;
			}
		}
	}

	//// check for rotation transformation
	//for (int i = 0; i < 511; i++)
	//{
	//	cout << middle_buffer[511][511 + i] << " " << output_buffer[max_length - i - 1][359] << endl;
	//}
	
	// Image Scale for output_buffer
	uint32** scaled_output_buffer;
	scaled_output_buffer = new uint32*[shape_h];
	for (int i = 0; i < shape_h; i++)
	{
		scaled_output_buffer[i] = new uint32[shape_w];
	}

	for (int h = 0; h < shape_h; h++)
	{
		float x = (float(h) + 0.5f) * float(max_length) / float(shape_h) - 0.5f;
		int fx = (int)x;
		x -= fx;

		int x1 = 1.f - x;
		int x2 = 1.f - x1;
		for (int w = 0; w < shape_w; w++)
		{
			float y = (float(w) + 0.5f) * 360.f / float(shape_w) - 0.5f;
			int fy = (int)y;
			y -= fy;

			int y1 = 1.f - y + 1.f;
			int y2 = 1.f - y1;
			
			scaled_output_buffer[h][w] = uint32(output_buffer[fx][fy] * x1 * y1 + output_buffer[fx + 1][fy] * x2 * y1 + \
				output_buffer[fx][fy + 1] * x1 * y2 + output_buffer[fx + 1][fy + 1] * x2 * y2);
			// cout << scaled_output_buffer[h][w] << endl;
		}
	}

	TIFF* output_tif = TIFFOpen("OutputTif.tif", "w");
	TIFFSetField(output_tif, TIFFTAG_IMAGEWIDTH, 360);
	TIFFSetField(output_tif, TIFFTAG_IMAGELENGTH, max_length);
	TIFFSetField(output_tif, TIFFTAG_SAMPLESPERPIXEL, input_channels);
	TIFFSetField(output_tif, TIFFTAG_BITSPERSAMPLE, input_bits);
	TIFFSetField(output_tif, TIFFTAG_DATATYPE, input_data_type);
	TIFFSetField(output_tif, TIFFTAG_PHOTOMETRIC, input_photo_metric);
	TIFFSetField(output_tif, TIFFTAG_XRESOLUTION, input_x_resolution);
	TIFFSetField(output_tif, TIFFTAG_YRESOLUTION, input_y_resolution);
	TIFFSetField(output_tif, TIFFTAG_RESOLUTIONUNIT, input_res_unit);
	TIFFSetField(output_tif, TIFFTAG_PLANARCONFIG, input_planar_config);
	
	if (input_planar_config == PLANARCONFIG_CONTIG)
	{
		cout << "Write Correct" << endl;
		for (row = 0; row < max_length; row++)
		{
			TIFFWriteScanline(output_tif, output_buffer[row], row);
		}
	}
	else if (input_planar_config == PLANARCONFIG_SEPARATE)
	{
		uint16 s;

		TIFFSetField(output_tif, TIFFTAG_SAMPLESPERPIXEL, input_nsamples);

		for (s = 0; s < input_nsamples; s++)
		{
			for (row = 0; row < input_height; row++)
			{
				TIFFWriteScanline(output_tif, output_buffer[row], row, s);
			}
		}
	}

	TIFF* scaled_output_tif = TIFFOpen("ScaledOutputTif.tif", "w");
	TIFFSetField(scaled_output_tif, TIFFTAG_IMAGEWIDTH, shape_w);
	TIFFSetField(scaled_output_tif, TIFFTAG_IMAGELENGTH, shape_h);
	TIFFSetField(scaled_output_tif, TIFFTAG_SAMPLESPERPIXEL, input_channels);
	TIFFSetField(scaled_output_tif, TIFFTAG_BITSPERSAMPLE, input_bits);
	TIFFSetField(scaled_output_tif, TIFFTAG_DATATYPE, input_data_type);
	TIFFSetField(scaled_output_tif, TIFFTAG_PHOTOMETRIC, input_photo_metric);
	TIFFSetField(scaled_output_tif, TIFFTAG_XRESOLUTION, input_x_resolution);
	TIFFSetField(scaled_output_tif, TIFFTAG_YRESOLUTION, input_y_resolution);
	TIFFSetField(scaled_output_tif, TIFFTAG_RESOLUTIONUNIT, input_res_unit);
	TIFFSetField(scaled_output_tif, TIFFTAG_PLANARCONFIG, input_planar_config);

	if (input_planar_config == PLANARCONFIG_CONTIG)
	{
		cout << "Write Correct" << endl;
		for (row = 0; row < shape_h; row++)
		{
			TIFFWriteScanline(scaled_output_tif, scaled_output_buffer[row], row);
		}
	}
	else if (input_planar_config == PLANARCONFIG_SEPARATE)
	{
		uint16 s;

		TIFFSetField(scaled_output_tif, TIFFTAG_SAMPLESPERPIXEL, input_nsamples);

		for (s = 0; s < input_nsamples; s++)
		{
			for (row = 0; row < shape_h; row++)
			{
				TIFFWriteScanline(scaled_output_tif, scaled_output_buffer[row], row, s);
			}
		}
	}

	_TIFFfree(input_buf);
	TIFFClose(input_tif);
	TIFFClose(output_tif);
	TIFFClose(scaled_output_tif);
	free(middle_buffer);
	free(output_buffer);
	free(scaled_output_buffer);
	return true;
}

bool TiffTransformAndScale(const char* input_file, int position_x, int position_y, int shape_w = 0, int shape_h = 0, bool is_scale = false) {
	TIFF* input_tif = TIFFOpen(input_file, "r");

	if (input_tif)
	{
		cout << "Read success!!" << endl;
	}
	else
	{
		cout << "Load input failed. Please input right .tif" << endl;
		return false;
	}

	int dir_count = 0;
	do
	{
		dir_count++;
	} while (TIFFReadDirectory(input_tif));

	uint32 input_width, input_height = 0;
	int input_channels = 0;
	int input_bits = 0;
	uint16 input_data_type = 0;
	int input_photo_metric = 0;
	float input_x_resolution, input_y_resolution = 0.f;
	uint16 input_planar_config = 0;
	uint16 input_res_unit = 0;
	uint16 input_nsamples = 0;

	TIFFGetField(input_tif, TIFFTAG_IMAGEWIDTH, &input_width);
	TIFFGetField(input_tif, TIFFTAG_IMAGELENGTH, &input_height);
	TIFFGetField(input_tif, TIFFTAG_SAMPLESPERPIXEL, &input_channels);
	TIFFGetField(input_tif, TIFFTAG_BITSPERSAMPLE, &input_bits);
	TIFFGetField(input_tif, TIFFTAG_DATATYPE, &input_data_type);
	TIFFGetField(input_tif, TIFFTAG_PHOTOMETRIC, &input_photo_metric);
	TIFFGetField(input_tif, TIFFTAG_XRESOLUTION, &input_x_resolution);
	TIFFGetField(input_tif, TIFFTAG_YRESOLUTION, &input_y_resolution);
	TIFFGetField(input_tif, TIFFTAG_RESOLUTIONUNIT, &input_res_unit);
	TIFFGetField(input_tif, TIFFTAG_PLANARCONFIG, &input_planar_config);

	int tag_c = TIFFGetTagListCount(input_tif);

	cout << "Tags show: (number of tags:" << tag_c << ")" << endl;
	cout << "numOffiles * width * height * channels * bitsPerChannel" << endl;
	cout << dir_count << " * " << input_width << " * " << input_height << " * " << input_channels << " * " << input_bits << endl;
	cout << "data type: " << input_data_type << " , photo metric: " << input_photo_metric << endl; // 1: gary image, 0 means the darkest.
	cout << "xresolution: " << input_x_resolution << " , y resolution " << input_y_resolution << ", resolution unit " << input_res_unit << endl;
	cout << "planar config: " << input_planar_config << endl;
	cout << endl;

	int max_length = ceil(sqrt(pow(max(int(input_width) - position_x, position_x), 2) + \
		pow(max(int(input_height) - position_y, position_y), 2)));
	cout << "Max length of output Tif is: " << max_length << endl;

	// Scanline-based Image I/O
	tdata_t input_buf;
	uint32 row;

	uint32** middle_buffer;
	middle_buffer = new uint32*[input_height];
	for (int i = 0; i < input_height; i++)
	{
		middle_buffer[i] = new uint32[input_width];
	}

	/*uint32** middle_buffer = (uint32**)malloc(sizeof(uint32*)*input_height);
	for (int i = 0; i < input_height; i++)
	{
		middle_buffer[i] = (uint32*)malloc(sizeof(uint32)*input_width);
	}*/

	uint32* temp_buffer;
	input_buf = _TIFFmalloc(TIFFStripSize(input_tif));
	if (input_planar_config == PLANARCONFIG_CONTIG)
	{
		cout << "Read Correct" << endl;
		for (row = 0; row < input_height; row++)
		{
			TIFFReadScanline(input_tif, input_buf, row);
			temp_buffer = (uint32*)input_buf;

			for (int col = 0; col < input_width; col++)
			{
				middle_buffer[row][col] = temp_buffer[col];
			}
		}
	}
	else if (input_planar_config == PLANARCONFIG_SEPARATE)
	{
		uint16 s;

		TIFFGetField(input_tif, TIFFTAG_SAMPLESPERPIXEL, &input_nsamples);

		for (s = 0; s < input_nsamples; s++)
		{
			for (row = 0; row < input_height; row++)
			{
				TIFFReadScanline(input_tif, input_buf, row, s);
				temp_buffer = (uint32*)input_buf;

				for (int col = 0; col < input_width; col++)
				{
					middle_buffer[row][col] = temp_buffer[col];
				}
			}
		}
	}

	uint32** output_buffer;
	output_buffer = new uint32*[max_length];
	for (int i = 0; i < max_length; i++)
	{
		output_buffer[i] = new uint32[360];
	}

	/*uint32** output_buffer = (uint32**)malloc(sizeof(uint32*)*max_length);
	for (int i = 0; i < max_length; i++)
	{
		output_buffer[i] = (uint32*)malloc(sizeof(uint32)*360);
	}*/

	myVec2d* ray_index;
	ray_index = (myVec2d*)malloc(sizeof(myVec2d)*max_length);

	for (int i = 1; i <= 360; i++)
	{
		rotatemy(ray_index, position_x, position_y, max_length, -1 * float(i) / 180.f*PI);
		for (int j = 0; j < max_length; j++)
		{
			int select_pos_x = ray_index[j].y;
			int select_pos_y = ray_index[j].x;

			if (select_pos_x < input_height && select_pos_y < input_width && select_pos_x > -1 && select_pos_y > -1)
			{
				output_buffer[max_length - j - 1][i - 1] = middle_buffer[select_pos_x][select_pos_y];

				//if (output_buffer[max_length - j - 1][i - 1] > 0 && output_buffer[max_length - j - 1][i - 1] < 2147483648)
				//{
				//	// cout << output_buffer[max_length - j - 1][i - 1] << endl;
				//	cout << select_pos_x << " " << select_pos_y << " " << output_buffer[max_length - j - 1][i - 1] << " " << middle_buffer[select_pos_x][select_pos_y] << endl;
				//}
			}
			else
			{
				output_buffer[max_length - j - 1][i - 1] = 0;
			}
		}
	}

	//// check for rotation transformation
	//for (int i = 0; i < 511; i++)
	//{
	//	cout << middle_buffer[511][511 + i] << " " << output_buffer[max_length - i - 1][359] << endl;
	//}

	TIFF* output_tif = TIFFOpen("OutputTif.tif", "w");
	TIFFSetField(output_tif, TIFFTAG_IMAGEWIDTH, 360);
	TIFFSetField(output_tif, TIFFTAG_IMAGELENGTH, max_length);
	TIFFSetField(output_tif, TIFFTAG_SAMPLESPERPIXEL, input_channels);
	TIFFSetField(output_tif, TIFFTAG_BITSPERSAMPLE, input_bits);
	TIFFSetField(output_tif, TIFFTAG_DATATYPE, input_data_type);
	TIFFSetField(output_tif, TIFFTAG_PHOTOMETRIC, input_photo_metric);
	TIFFSetField(output_tif, TIFFTAG_XRESOLUTION, input_x_resolution);
	TIFFSetField(output_tif, TIFFTAG_YRESOLUTION, input_y_resolution);
	TIFFSetField(output_tif, TIFFTAG_RESOLUTIONUNIT, input_res_unit);
	TIFFSetField(output_tif, TIFFTAG_PLANARCONFIG, input_planar_config);

	if (input_planar_config == PLANARCONFIG_CONTIG)
	{
		cout << "Write Correct" << endl;
		for (row = 0; row < max_length; row++)
		{
			TIFFWriteScanline(output_tif, output_buffer[row], row);
		}
	}
	else if (input_planar_config == PLANARCONFIG_SEPARATE)
	{
		uint16 s;

		TIFFSetField(output_tif, TIFFTAG_SAMPLESPERPIXEL, input_nsamples);

		for (s = 0; s < input_nsamples; s++)
		{
			for (row = 0; row < input_height; row++)
			{
				TIFFWriteScanline(output_tif, output_buffer[row], row, s);
			}
		}
	}

	// Image Scale for output_buffer
	if (is_scale)
	{
		uint32** scaled_output_buffer;
		scaled_output_buffer = new uint32*[shape_h];
		for (int i = 0; i < shape_h; i++)
		{
			scaled_output_buffer[i] = new uint32[shape_w];
		}

		for (int h = 0; h < shape_h; h++)
		{
			float x = (float(h) + 0.5f) * float(max_length) / float(shape_h) - 0.5f;
			int fx = (int)x;
			x -= fx;

			int x1 = 1.f - x;
			int x2 = 1.f - x1;
			for (int w = 0; w < shape_w; w++)
			{
				float y = (float(w) + 0.5f) * 360.f / float(shape_w) - 0.5f;
				int fy = (int)y;
				y -= fy;

				int y1 = 1.f - y + 1.f;
				int y2 = 1.f - y1;

				scaled_output_buffer[h][w] = uint32(output_buffer[fx][fy] * x1 * y1 + output_buffer[fx + 1][fy] * x2 * y1 + \
					output_buffer[fx][fy + 1] * x1 * y2 + output_buffer[fx + 1][fy + 1] * x2 * y2);
				// cout << scaled_output_buffer[h][w] << endl;
			}
		}

		TIFF* scaled_output_tif = TIFFOpen("ScaledOutputTif.tif", "w");
		TIFFSetField(scaled_output_tif, TIFFTAG_IMAGEWIDTH, shape_w);
		TIFFSetField(scaled_output_tif, TIFFTAG_IMAGELENGTH, shape_h);
		TIFFSetField(scaled_output_tif, TIFFTAG_SAMPLESPERPIXEL, input_channels);
		TIFFSetField(scaled_output_tif, TIFFTAG_BITSPERSAMPLE, input_bits);
		TIFFSetField(scaled_output_tif, TIFFTAG_DATATYPE, input_data_type);
		TIFFSetField(scaled_output_tif, TIFFTAG_PHOTOMETRIC, input_photo_metric);
		TIFFSetField(scaled_output_tif, TIFFTAG_XRESOLUTION, input_x_resolution);
		TIFFSetField(scaled_output_tif, TIFFTAG_YRESOLUTION, input_y_resolution);
		TIFFSetField(scaled_output_tif, TIFFTAG_RESOLUTIONUNIT, input_res_unit);
		TIFFSetField(scaled_output_tif, TIFFTAG_PLANARCONFIG, input_planar_config);

		if (input_planar_config == PLANARCONFIG_CONTIG)
		{
			cout << "Write Correct" << endl;
			for (row = 0; row < shape_h; row++)
			{
				TIFFWriteScanline(scaled_output_tif, scaled_output_buffer[row], row);
			}
		}
		else if (input_planar_config == PLANARCONFIG_SEPARATE)
		{
			uint16 s;

			TIFFSetField(scaled_output_tif, TIFFTAG_SAMPLESPERPIXEL, input_nsamples);

			for (s = 0; s < input_nsamples; s++)
			{
				for (row = 0; row < shape_h; row++)
				{
					TIFFWriteScanline(scaled_output_tif, scaled_output_buffer[row], row, s);
				}
			}
		}

		TIFFClose(scaled_output_tif);
		free(scaled_output_buffer);
	}

	_TIFFfree(input_buf);
	TIFFClose(input_tif);
	TIFFClose(output_tif);
	
	free(middle_buffer);
	free(output_buffer);
	
	return true;
}

void ReadResult(const char* input_file)
{
	TIFF* input_tif = TIFFOpen(input_file, "r");

	uint32 input_width, input_height = 0;
	TIFFGetField(input_tif, TIFFTAG_IMAGEWIDTH, &input_width);
	TIFFGetField(input_tif, TIFFTAG_IMAGELENGTH, &input_height);

	tdata_t input_buf;
	uint32 row;

	uint32** middle_buffer;
	middle_buffer = new uint32*[input_height];
	for (int i = 0; i < input_height; i++)
	{
		middle_buffer[i] = new uint32[input_width];
	}

	input_buf = _TIFFmalloc(TIFFStripSize(input_tif));
	for (row = 0; row < input_height; row++)
	{
		TIFFReadScanline(input_tif, input_buf, row);
		middle_buffer[row] = (uint32*)input_buf;
	}

	for (int i = 0; i < 511; i++)
	{
		cout << middle_buffer[input_height - i - 1][359] << endl;
	}
}

//int main() {
//	// string input_file = "InputImage.jpg";
//	string input_file = "InputImage.tif";
//	// string input_file = "OutputTif.tif";
//	int position_x = 468;
//	int position_y = 881;
//	/*int position_x = 512;
//	int position_y = 512;*/
//	int shape_w = 256;
//	int shape_h = 256;
//
//	// ImageTransform(input_file.c_str(), position_x, position_y, shape_w, shape_h);
//
//	// TiffTransform(input_file.c_str(), position_x, position_y, shape_w, shape_h);
//
//	TiffTransformReal(input_file.c_str(), position_x, position_y, shape_w, shape_h);
//
//	// ReadResult(input_file.c_str());
//
//	system("pause");
//	return 0;
//}

int main(int argc, char* argv[]) {
	if(argc == 5 && atoi(argv[4]) == 0)
	{
		char* input_file = argv[1];
		char* input_position_x = argv[2];
		char* input_position_y = argv[3];

		int position_x = atoi(input_position_x);
		int position_y = atoi(input_position_y);

		TiffTransformAndScale(input_file, position_x, position_y);
	}
	else if(argc == 7 && atoi(argv[4]) == 1)
	{
		char* input_file = argv[1];
		char* input_position_x = argv[2];
		char* input_position_y = argv[3];
		char* output_shape_w = argv[5];
		char* output_shape_h = argv[6];

		int position_x = atoi(input_position_x);
		int position_y = atoi(input_position_y);
		int shape_w = atoi(output_shape_w);
		int shape_h = atoi(output_shape_h);

		TiffTransformAndScale(input_file, position_x, position_y, shape_w, shape_h, true);
	}
	else
	{
		cout << "Wrong Input, Execute example: \'./ImageTransform.exe ./InputImage.tif 468 881 1 512 512\'(With scale) " << endl;
		cout << "Or: \'./ImageTransform.exe ./InputImage.tif 468 881 0\'(No scale) " << endl;
		cout << "It means, input image ./InputImage.tif; pixel position x=25, y=25; output shape, w=512, h=512." << endl;
	}

	system("pause");
	return 0;
}