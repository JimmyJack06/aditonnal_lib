// additional_lib.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <opencv2/opencv.hpp>
#include <ipp.h>
#include <mkl.h>

using namespace cv;
using namespace std;

#pragma region <图像批量存取相关> 
const int IMG_COUNT = 100;    //<图像序列总数
/**
* @brief <B>Destription:</B><br> 批量重命名图像（亦可用在无法打开图片数据间接使用opencv读存）
*
* @param[in] start_No 名字开始序号
* @param[in] end_No 名字结束序号
*
* @brief 使用时注意更改路径！！！
*/
void tranSaveImg(int start_No, int end_No)
{
	Mat testImg;
	char file[IMG_COUNT];
	char newfile[IMG_COUNT];
	//int file_No = 0;
	for (int i = start_No; i <= end_No; i++)
	{
		sprintf(file, "repeat_testImg/%u.bmp", i);
		sprintf(newfile, "repeat_testImg/%u(1).bmp", i);
		testImg = imread(file, 0);
		imwrite(newfile, testImg);
	}
}

/**
* @brief <B>Destription:</B><br> 批量读取图像序列
*
* @param Src 指针数组，存图像序列地址Src[0]表示首张图首地址。
* @param[in] start_No 名字开始序号
* @param[in] end_No 名字结束序号
* @param[out] rows 返回图像高
* @param[out] cols 返回图像宽
*
* @brief 使用时注意更改路径！！！
*/
void autoReadImg(unsigned char **&Src, int start_No, int end_No, int &rows, int &cols)
{
	Src = (unsigned char **)malloc(sizeof(unsigned char*) * IMG_COUNT);
	char file[IMG_COUNT];
	int img_No = 0;
	for (int i = start_No; i <= end_No; i++)
	{
		sprintf(file, "%u.bmp", i);
		Mat midMat = imread(file, 0);			//局部变量！！！！
		unsigned char *temp = midMat.data;
		Src[img_No] = (unsigned char *)malloc(sizeof(unsigned char) * midMat.rows * midMat.cols);
		for (int y = 0; y < midMat.rows; y++)
		{
			for (int x = 0; x < midMat.cols; x++)
			{
				Src[img_No][y * midMat.cols + x] = temp[y * midMat.cols + x];
			}
		}
		img_No++;
		//imshow(file, midMat);
		//waitKey(0);
		rows = midMat.rows;
		cols = midMat.cols;
	}
}

/**
* @brief <B>Destription:</B><br> 对应批量读取图像序列函数，释放指针数组
*
* @param Src 指针数组，存图像序列地址Src[0]表示首张图首地址。
*
* @brief 使用时注意更改路径！！！
*/
void freeMat(unsigned char **Src)
{
	free(Src);
}
#pragma endregion<图像批量存取相关> 

#pragma region<文件读写相关> 

//FREOPEN
#ifdef FREOPEN
	freopen("test_out.out", "w", stdout);
#endif // FREOPEN

///////////			C++			///////////////
#include <fstream>
void outputFileCpp()
{
	ofstream outfile("data.txt");
	for (int i = 0; i < 10; i++)
		outfile << i << " ";
	cout << "ok" << endl;
	getchar();
	outfile.close();
}

void inputFileCpp()
{
	int data;
	ifstream infile("data.txt");
	for (int i = 0; i < 10; i++)
	{
		infile >> data;
		cout << data << " ";
	}
	infile.close();
}

///////////		C			/////////////////////////////
void outputFileC()
{
	FILE *pOutfile = fopen("data.txt", "w+");
	int a = 1;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(pOutfile, "%d ", a++ );
		}
		fprintf(pOutfile, "\n");
	}
	fclose(pOutfile);
}

void inputFileC()
{
	FILE *pInfile = fopen("data.txt", "w+");
	int a = 1;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fscanf(pInfile, "%d ", &a);
		}
		//fprintf(pInfile, "\n");
	}
	fclose(pInfile);
}

#pragma endregion<文件读写相关> 

#pragma region<动态申请空间> 
void newSpace()
{
	int *p1 = new int;
	int *p2 = new int[1000];
	delete p1;
	delete[] p2;

	int *p3 = (int*)malloc(sizeof(int) * 1000);
	free(p3);
}
#pragma endregion<动态申请空间> 

#pragma  region<IPP related>
void threshold_OTSU()
{
	Mat srcImg = imread("2.bmp", 0);
	unsigned char *pSrc = srcImg.data;
	cout << (double)(*pSrc) << endl;
	int srcStep = srcImg.cols * sizeof(unsigned char);
	IppiSize roiSize = { srcImg.cols,srcImg.rows };
	unsigned char *pThreshold = pSrc;
	ippiComputeThreshold_Otsu_8u_C1R( pSrc, srcStep,  roiSize,  pThreshold);	//*pSrc change	
	cout << (double)(*pThreshold) << endl;
	//ofstream outImg("outImg.txt");
	//for (int i = 0; i < srcImg.rows; i++)
	//{
	//	for (int j = 0; j < srcImg.cols; j++)
	//	{
	//		outImg << (double)pSrc[i * srcImg.cols + j]<<" ";
	//	}
	//	outImg << endl;
	//}
	//outImg.close();
	unsigned char OSTUthreshold = *pThreshold;		//important
	for (int i = 0; i < srcImg.rows; i++)
	{
		for (int j = 0; j < srcImg.cols; j++)
		{
			if ((pSrc[i * srcImg.cols + j]) <= OSTUthreshold)
			{
				//cout << i << endl;
				pSrc[i * srcImg.cols + j] = 0;
			}
			else
			{
				pSrc[i * srcImg.cols + j] = 255;
			}
		}
	}
	imwrite("thresholdedImg.bmp", srcImg);
}

#pragma endregion<IPP related>

#pragma region<边缘检测相关>
#define PI 3.1415926
typedef struct
{
	signed short  x;
	signed short  y;
} IMG_COORD;

typedef struct
{
	float   x;
	float   y;
} IMG_RCOORD;
typedef struct 
{
	IMG_COORD xyInteger;
	IMG_RCOORD xyDecimal;
	int gradient;
	float angle;
}edgeInformation;

void SobelFilter_8u16s_C1_3x3or5x5(unsigned char *pSrc, IppiSize roiSize, int kernalSize,
	signed short *&pDst, float *&pAngle, int srcStep)
{
	IppiMaskSize mask = ippMskSize3x3;
	if (kernalSize == 3)
	{
		mask = ippMskSize3x3;
	}
	else
	{
		if (kernalSize == 5)
		{
			mask = ippMskSize5x5;
		}
	}

	IppiBorderType bordertype = ippBorderRepl; //Border is replicated from the edge pixels
	signed short *pHoriz, *pVert;
	//int srcStep = roiSize.width * sizeof(unsigned char);
	int dstStep = roiSize.width * sizeof(signed short);
	int angleStep = roiSize.width * sizeof(float);
	int bufferSize;
	int bufLen = roiSize.width * roiSize.height;
	IppStatus statusVert, statusHoriz, status;
	unsigned char *pBuffer;
	IppNormType normType = ippNormL2;//input gradient magnitude

	pVert = (signed short *)malloc(sizeof(Ipp16s)*bufLen);
	pHoriz = (signed short *)malloc(sizeof(signed short)*bufLen);
	pAngle = (float *)malloc(sizeof(float)*bufLen);
	pDst = (signed short *)malloc(sizeof(signed short)*bufLen);

	ippiGradientVectorGetBufferSize(roiSize, mask, ipp16s, 1, &bufferSize);
	pBuffer = (unsigned char *)malloc(bufferSize);
	ippiGradientVectorSobel_8u16s_C1R(pSrc, srcStep, pVert, dstStep, pHoriz, dstStep, pDst, dstStep, pAngle, angleStep, roiSize, mask, normType, bordertype, NULL, pBuffer);

	free(pVert);
	free(pHoriz);
	free(pBuffer);
}

void edge_detection(unsigned char *srcRoi, int roiRows, int roiCols, int threshold, int kernalSize, signed short *&dstRoi,
	Ipp16u *dstRoiE, float *&angAll, edgeInformation *&edgeArray, int &sn)//std::vector<edgeInformation> &edgeInfor)
{
	std::vector<edgeInformation> edgeInfor;
	edgeInformation edInf;

	int k = 0;//记录边缘点的个数
	int k1;//抛物线拟合的三个已知点
	int k2;
	int k3;
	float deci;//抛物线拟合顶点的小数部分，即对应的亚像素
	float sumx = 0;//边缘点的x坐标之和
	float sumy = 0;
	int numberChannels = 1; //the source image is single channel

	IppiSize dstRoiSize = { roiCols,roiRows };

	SobelFilter_8u16s_C1_3x3or5x5(srcRoi, dstRoiSize, kernalSize, dstRoi, angAll, roiCols);

	//把角度由[-pi，pi]变为[0，360]
	for (int i = 0; i < roiRows; i++)
	{
		for (int j = 0; j < roiCols; j++)
		{
			angAll[j + i * roiCols] = (float)(180 - angAll[j + i * roiCols] / PI * 180);
		}
	}

	/*FILE *sx;
	sx = fopen("E:\\project03\\learn\\a1sobel.txt", "w");
	FILE *ang;
	ang = fopen("E:\\project03\\learn\\a1ang.txt", "w");
	for (int i = 0; i<roiRows; i++)
	{
	for (int j = 0; j < roiCols; j++)
	{

	fprintf(sx, "%3d   ", dstRoi[j+i*roiCols]);
	fprintf(ang, "%f   ", angAll[j+i*roiCols]);
	}
	fprintf(sx,"\n");
	fprintf(ang, "\n");
	}*/

	//fclose(ang);

	for (int i = 1; i<roiRows - 1; i++)
	{
		for (int j = 1; j<roiCols - 1; j++)
		{
			if (dstRoi[j + i*roiCols] > threshold)
			{
				//angAll[j + i * roiCols] = (float)180 - angAll[j + i * roiCols] / PI * 180;
				if ((angAll[j + i*roiCols]>22.5) && (angAll[j + i*roiCols]<67.5))
				{
					if ((dstRoi[j + i*roiCols] > dstRoi[j - 1 + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + (i + 1)*roiCols]))
					{
						k1 = dstRoi[j - 1 + (i - 1)*roiCols];
						k2 = dstRoi[j + i*roiCols];
						k3 = dstRoi[(j + 1) + (i + 1)*roiCols];
						deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

						edInf.xyInteger.x = j;
						edInf.xyInteger.y = i;
						edInf.xyDecimal.x = j + deci;
						edInf.xyDecimal.y = i + deci;
						edInf.gradient = dstRoi[j + i*roiCols];
						edInf.angle = angAll[j + i*roiCols];
						edgeInfor.push_back(edInf);
						k++;
					}
				}
				else
				{
					if ((angAll[j + i*roiCols]>202.5) && (angAll[j + i*roiCols]<247.5))
					{
						if ((dstRoi[j + i*roiCols] > dstRoi[j - 1 + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + (i + 1)*roiCols]))
						{
							k3 = dstRoi[j - 1 + (i - 1)*roiCols];
							k2 = dstRoi[j + i*roiCols];
							k1 = dstRoi[(j + 1) + (i + 1)*roiCols];
							deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

							edInf.xyInteger.x = j;
							edInf.xyInteger.y = i;
							edInf.xyDecimal.x = j - deci;
							edInf.xyDecimal.y = i - deci;
							edInf.gradient = dstRoi[j + i*roiCols];
							edInf.angle = angAll[j + i*roiCols];
							edgeInfor.push_back(edInf);
							k++;
						}
					}
					else
					{
						if ((angAll[j + i*roiCols]>112.5) && (angAll[j + i*roiCols]<157.5))
						{

							if ((dstRoi[j + i*roiCols] > dstRoi[(j + 1) + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j - 1) + (i + 1)*roiCols]))
							{
								k1 = dstRoi[(j + 1) + (i - 1)*roiCols];
								k2 = dstRoi[j + i*roiCols];
								k3 = dstRoi[(j - 1) + (i + 1)*roiCols];
								deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

								edInf.xyInteger.x = j;
								edInf.xyInteger.y = i;
								edInf.xyDecimal.x = j - deci;
								edInf.xyDecimal.y = i + deci;
								edInf.gradient = dstRoi[j + i*roiCols];
								edInf.angle = angAll[j + i*roiCols];
								edgeInfor.push_back(edInf);
								k++;
							}
						}
						else
						{
							if ((angAll[j + i*roiCols]>292.5) && (angAll[j + i*roiCols]<337.5))
							{
								if ((dstRoi[j + i*roiCols] > dstRoi[(j + 1) + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j - 1) + (i + 1)*roiCols]))
								{
									k3 = dstRoi[(j + 1) + (i - 1)*roiCols];
									k2 = dstRoi[j + i*roiCols];
									k1 = dstRoi[(j - 1) + (i + 1)*roiCols];
									deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

									edInf.xyInteger.x = j;
									edInf.xyInteger.y = i;
									edInf.xyDecimal.x = j + deci;
									edInf.xyDecimal.y = i - deci;
									edInf.gradient = dstRoi[j + i*roiCols];
									edInf.angle = angAll[j + i*roiCols];
									edgeInfor.push_back(edInf);
									k++;
								}
							}
							else
							{
								if (((angAll[j + i*roiCols] >= -1) && (angAll[j + i*roiCols] <= 22.5)) || ((angAll[j + i*roiCols] >= 337.5) && (angAll[j + i*roiCols] <= 361)))
								{
									if ((dstRoi[j + i*roiCols] > dstRoi[(j - 1) + i*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + i*roiCols]))
									{
										k1 = dstRoi[(j - 1) + i*roiCols];
										k2 = dstRoi[j + i*roiCols];
										k3 = dstRoi[(j + 1) + i*roiCols];
										deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

										edInf.xyInteger.x = j;
										edInf.xyInteger.y = i;
										edInf.xyDecimal.x = j + deci;
										edInf.xyDecimal.y = i;
										edInf.gradient = dstRoi[j + i*roiCols];
										edInf.angle = angAll[j + i*roiCols];
										edgeInfor.push_back(edInf);
										k++;
									}
								}
								else
								{
									if ((angAll[j + i*roiCols] <= 202.5) && (angAll[j + i*roiCols] >= 157.5))
									{
										if ((dstRoi[j + i*roiCols] > dstRoi[(j - 1) + i*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[(j + 1) + i*roiCols]))
										{
											k3 = dstRoi[(j - 1) + i*roiCols];
											k2 = dstRoi[j + i*roiCols];
											k1 = dstRoi[(j + 1) + i*roiCols];
											deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

											edInf.xyInteger.x = j;
											edInf.xyInteger.y = i;
											edInf.xyDecimal.x = j - deci;
											edInf.xyDecimal.y = i;
											edInf.gradient = dstRoi[j + i*roiCols];
											edInf.angle = angAll[j + i*roiCols];
											edgeInfor.push_back(edInf);
											k++;
										}
									}
									else
									{
										if ((angAll[j + i*roiCols] >= 67.5) && (angAll[j + i*roiCols] <= 112.5))
										{

											if ((dstRoi[j + i*roiCols] > dstRoi[j + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[j + (i + 1)*roiCols]))
											{
												k1 = dstRoi[j + (i - 1)*roiCols];
												k2 = dstRoi[j + i*roiCols];
												k3 = dstRoi[j + (i + 1)*roiCols];
												deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

												edInf.xyInteger.x = j;
												edInf.xyInteger.y = i;
												edInf.xyDecimal.x = j;
												edInf.xyDecimal.y = i + deci;
												edInf.gradient = dstRoi[j + i*roiCols];
												edInf.angle = angAll[j + i*roiCols];
												edgeInfor.push_back(edInf);
												k++;
											}
										}
										else
										{
											if ((angAll[j + i*roiCols] >= 247.5) && (angAll[j + i*roiCols] <= 292.5))
											{
												if ((dstRoi[j + i*roiCols] > dstRoi[j + (i - 1)*roiCols]) && (dstRoi[j + i*roiCols] >= dstRoi[j + (i + 1)*roiCols]))
												{
													k3 = dstRoi[j + (i - 1)*roiCols];
													k2 = dstRoi[j + i*roiCols];
													k1 = dstRoi[j + (i + 1)*roiCols];
													deci = (k3 - k1) / ((float)2.0*(2.0*k2 - k1 - k3));

													edInf.xyInteger.x = j;
													edInf.xyInteger.y = i;
													edInf.xyDecimal.x = j;
													edInf.xyDecimal.y = i - deci;
													edInf.gradient = dstRoi[j + i*roiCols];
													edInf.angle = angAll[j + i*roiCols];
													edgeInfor.push_back(edInf);
													k++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	sn = k;
	for (int t = 0; t < roiCols*roiRows; t++)//二值图像，所有像素先都赋值为0，边缘点赋值255
	{
		dstRoiE[t] = 0;
	}
	for (int q = 0; q < k; q++)
	{
		sumx = sumx + edgeInfor[q].xyDecimal.x;
		sumy = sumy + edgeInfor[q].xyDecimal.y;
		dstRoiE[edgeInfor[q].xyInteger.x + edgeInfor[q].xyInteger.y * roiCols] = 255;
	}
	//printf("%d\n",k);

	/*
	FILE *db;
	db = fopen("E:\\ProjectCameraVerticality\\de.txt","w");
	for (int i = 0; i < roiRows; i++)
	{
	for (int j = 0; j < roiCols; j++)
	{
	fprintf(db,"%d   ",dstRoiE[j+i*roiCols]);
	}
	fprintf(db,"\n");
	}
	fclose(db);
	*/
	//以数组的方式传出边缘信息
	edgeArray = (edgeInformation*)malloc(k*sizeof(edgeInformation));
	for (int i = 0; i < k; i++)
	{
		edgeArray[i] = edgeInfor[i];
	}
	//printf("%d\n",k);

	/*FILE *e;
	e = fopen("E:\\project03\\e.txt", "w");
	for (int i = 0; i < k; i++)
	{
	fprintf(e,"%f   %f   \n", edgeInfor[i].xyDecimal.x, edgeInfor[i].xyDecimal.y);
	}*/

	//ecout.x = sumx / k;
	//ecout.y = sumy / k;
	//ecout.en = k;
}
#pragma endregion<边缘检测相关>


#pragma region<Fitting related>
void fitCircle(double *X, double *Y, int _num, double &cx, double &cy, double &radius)
{
	//dst vars
	bool iteration = false;
	double centerX, centerY, T;		//T = centerX * centerX + centerY * centerY - radius * radius

	//init vars
	double *matA = (double*)malloc(sizeof(double) * 3 * _num);
	memset(matA, 0, sizeof(double) * 3 * _num);
	double *matB = (double*)malloc(sizeof(double) * _num);
	memset(matB, 0, sizeof(double) *  _num);

	//set vals
	for (int i = 0; i < _num; i++)
	{
		matA[i * 3 + 0] = -2 * X[i];
		matA[i * 3 + 1] = -2 * Y[i];
		matA[i * 3 + 2] = 1;

		matB[i] = -X[i] * X[i] - Y[i] * Y[i];
	}
	
	if (!iteration)		//do not iterate
	{
		LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', _num, 3, 1, matA, 3, matB, 1);
		centerX = matB[0];
		centerY = matB[1];
		T = matB[2];

	}

	cx = centerX;
	cy = centerY;
	radius = sqrtf(centerX * centerX + centerY * centerY - T);

	//free
	free(matA);
	free(matB);

}

void fitLine(double * X, double *Y, int _num, double &lineA, double &lineB, double &lineC)
{
	bool iteration = false;
	bool slope_exist = true;
	

	//judege if slope is exist(calculate Xi 's variation)
	double veryMin = 5;
	double sum = 0;
	double variation = 0;
	for (int i = 0; i < _num; i++)
	{
		sum += X[i];
	}
	double Xmean = sum / _num;
	sum = 0;
	for (int i = 0; i < _num; i++)
	{
		sum += (X[i] - Xmean) * (X[i] - Xmean);
	}
	variation = sum / (_num - 1);
	if (variation < veryMin)		slope_exist = false;

	if (slope_exist)
	{
		double *matA = (double*)malloc(sizeof(double) * 2 * _num);
		memset(matA, 0, sizeof(double) * 2 * _num);
		double *matB = (double*)malloc(sizeof(double) * _num);
		memset(matB, 0, sizeof(double) *  _num);

		for (int i = 0; i < _num; i++)
		{
			matA[i * 2 + 0] = X[i];
			matA[i * 2 + 1] = 1;
			matB[i] = Y[i];
		}

		if (!iteration)		//do not iterate
		{
			LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', _num, 2, 1, matA, 2, matB, 1);
		}

		double k = matB[0];
		double b = matB[1];

		lineA = k;
		lineB = -1;
		lineC = b;

		free(matA);
		free(matB);
	}
	else
	{
		lineB = 0;
		lineA = 1;
		lineC = -Xmean;
	}
	
}
#pragma endregion<Fitting related>

#pragma region<异常处理>
int throwAssertTest()
{
	int a;
	cout << "please input a:" << endl;
	cin >> a;
	
	try {
		if (a == 0)
			throw("a cannot be 0!");
		if (a == 1)
			throw("a cannot be 1!");
	}
	catch (char* s)
	{
		printf("%s\n", s);
		//exit(0);
	}

	assert(a != -1);
	return 0;
}

//class ErrorResult
//{
//	int status;
//	char* errorMessage;
//};

typedef enum 
{
	noError = 0,
	Error1 = 1,
	Error2 = 2,
	Error3 = 3
}errorType;

char * getErrorMessage(errorType status)
{
	char *s = "noError";
	switch (status)
	{
	case Error1:
		s = "Error1";
		break;
	case Error2:
		s = "Error2";
		break;
	case Error3:
		s = "Error3";
		break;
	default:
		break;
	}
	return s;
}

errorType errorProRet()
{
	errorType status = Error1;
	printf("The status is %s\n", getErrorMessage(status));
	return status;
}
#pragma endregion<异常处理>


//test edgeDetect and circle fitting
int test1()
{
	Mat srcImg = imread("circle.bmp", 0);
	unsigned char *pSrc = srcImg.data;

	int gradThre = 600;
	int kernSize = 5;
	signed short *dstGrad = NULL;
	Ipp16u *dstRoiE;
	dstRoiE = (Ipp16u*)malloc(srcImg.rows * srcImg.cols *sizeof(Ipp16u));
	float *angle = NULL;

	//vector<edgeInformation> edgeInfor;
	edgeInformation *edgeArray;
	int sn = 0;
	edge_detection(pSrc, srcImg.rows, srcImg.cols, gradThre, kernSize, dstGrad, dstRoiE, angle, edgeArray, sn);

	//FILE *egA;
	//egA = fopen("egA.txt", "w");
	//for (int i = 0; i < sn; i++)
	//{
	//fprintf(egA, "%d   %d\n", edgeArray[i].xyInteger.x, edgeArray[i].xyInteger.y);
	//}
	//fclose(egA);

	double *X = (double *)malloc(sizeof(double) * sn);
	double *Y = (double *)malloc(sizeof(double) * sn);
	for (int k = 0; k < sn; k++)
	{
		X[k] = edgeArray[k].xyDecimal.x;
		Y[k] = edgeArray[k].xyDecimal.y;
	}

	double cx = 0;
	double cy = 0;
	double radius = 0;
	fitCircle(X, Y, sn, cx, cy, radius);
	
	cout << "cx = " << cx << endl << "cy = " << cy << endl << "radius = " << radius << endl;
	free(X);
	free(Y);

	return 0;
}

//test edgeDetect and line fitting
int test2()
{
	Mat srcImg = imread("line(x = 218).bmp", 0);
	unsigned char *pSrc = srcImg.data;

	int gradThre = 600;
	int kernSize = 5;
	signed short *dstGrad = NULL;
	Ipp16u *dstRoiE;
	dstRoiE = (Ipp16u*)malloc(srcImg.rows * srcImg.cols *sizeof(Ipp16u));
	float *angle = NULL;

	//vector<edgeInformation> edgeInfor;
	edgeInformation *edgeArray;
	int sn = 0;
	edge_detection(pSrc, srcImg.rows, srcImg.cols, gradThre, kernSize, dstGrad, dstRoiE, angle, edgeArray, sn);

	//FILE *egA;
	//egA = fopen("egA.txt", "w");
	//for (int i = 0; i < sn; i++)
	//{
	//fprintf(egA, "%d   %d\n", edgeArray[i].xyInteger.x, edgeArray[i].xyInteger.y);
	//}
	//fclose(egA);

	double *X = (double *)malloc(sizeof(double) * sn);
	double *Y = (double *)malloc(sizeof(double) * sn);
	for (int k = 0; k < sn; k++)
	{
		X[k] = edgeArray[k].xyDecimal.x;
		Y[k] = edgeArray[k].xyDecimal.y;
	}

	double lineA = 0;
	double lineB = 0;
	double lineC = 0;
	fitLine(X, Y, sn, lineA, lineB, lineC);

	cout << "line equation: " << lineA << "x + " << lineB << "y + " << lineC << " = 0" << endl;

	free(X);
	free(Y);
	return 0;
}

//int findCentroid(unsigned short *pSrc[], int imgNum, int srcRows, int srcCols, unsigned short threshold,float &xCentroid, float &yCentroid)
//{
//	//threshold
//	int srcStep = srcCols * sizeof(Ipp16u);
//	IppiSize roiSize = { srcCols,srcRows };
//	
//	//add image
//	Ipp32f* pSumImg = (Ipp32f*)malloc(srcRows * srcCols * sizeof(Ipp32f));
//	int tempStep = srcCols * sizeof(Ipp32f);
//
//	ippiConvert_16u32f_C1R(pSrc[0], srcStep, pSumImg, tempStep, roiSize);
//	for (int i = 1; i < imgNum; i++)
//	{
//		Ipp32f* pTemp = (Ipp32f*)malloc(srcRows * srcCols * sizeof(Ipp32f));
//		ippiConvert_16u32f_C1R(pSrc[i], srcStep, pTemp, tempStep, roiSize);
//		ippiAdd_32f_C1IR(pTemp, tempStep, pSumImg, tempStep, roiSize);
//		free(pTemp);
//	}
//	/*FILE *fTest;
//	fTest = fopen("pSumImg.txt", "w");
//	for (int i = 0; i < srcRows; i++)
//	{
//		for (int j = 0; j < srcCols; j++)
//		{
//			fprintf(fTest, "%f  ", (float)pSumImg[j + i*srcCols]);
//		}
//		fprintf(fTest, "\n");
//	}
//	fclose(fTest);*/
//
//	//threshold
//	Ipp32f *pThreshImg = (Ipp32f*)malloc(srcRows * srcCols * sizeof(Ipp32f));
//	ippiThreshold_32f_C1R(pSumImg, tempStep, pThreshImg, tempStep, roiSize, threshold, ippCmpGreater);
//	FILE *fTest;
//	fTest = fopen("pThreshImg.txt", "w");
//	for (int i = 0; i < srcRows; i++)
//	{
//		for (int j = 0; j < srcCols; j++)
//		{
//			fprintf(fTest, "%f  ", (float)pThreshImg[j + i*srcCols]);
//		}
//		fprintf(fTest, "\n");
//	}
//	fclose(fTest);
//
//	free(pSumImg);
//	return 0;
//}
//
//int testC5()
//{
//	int imgNum = 2;
//	char filename[100];
//	unsigned short **queSrc = (unsigned short **)malloc(sizeof(unsigned short*) * imgNum);
//	int srcRows, srcCols;
//	for (int _num = 0; _num < imgNum; _num++)
//	{
//		sprintf(filename, "%u.tif", _num);
//		Mat temp = imread(filename, -1);
//
//		queSrc[_num] = (unsigned short *)malloc(sizeof(unsigned short) * temp.rows * temp.cols);
//		for (int i = 0; i < temp.rows; i++)
//		{
//			for (int j = 0; j < temp.cols; j++)
//			{
//				queSrc[_num][i * temp.cols + j] = temp.at<ushort>(i, j);
//			}
//		}
//		srcCols = temp.cols;
//		srcRows = temp.rows;
//	}
//
//	unsigned short _threshold = 1000;
//	float x_cen, y_cen;
//	
//	findCentroid(queSrc, imgNum, srcRows, srcCols, _threshold, x_cen, y_cen);
//
//
//	free(queSrc);
//	return 0;
//}

//int main()
//{
//	threshold_OTSU();
//	return 0;
//}

#ifdef _DEBUG  
#define New   new(_NORMAL_BLOCK, __FILE__, __LINE__)  
#endif  

#define CRTDBG_MAP_ALLOC    
#include <stdlib.h>    
#include <crtdbg.h>    
//在入口函数中包含 _CrtDumpMemoryLeaks();    
//即可检测到内存泄露  

//以如下测试函数为例：  
int main()
{
	char* pChars = New char[10];
	_CrtDumpMemoryLeaks();
	delete[]pChars;
	return 0;
}