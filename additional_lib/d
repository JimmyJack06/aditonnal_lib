#pragma once
#include <ViType.h>
#include <memory.h>

#ifdef  DLL_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __declspec(dllimport)   
#endif

#define MAX_FILE_NAME (2560)
#define MAX_CHANNEL_NUM (100)
#define PLATTE_SIZE (1024)
#define JPG_SIZE_LIMITATION (65000)

enum IMAGE_TYPE
{
	BMP_TYPE = 0,
	TIFF_TYPE = 1,
	JPG_TYPE = 2,
	PNG_TYPE = 3,
	UNKNOWN_TYPE = 4
};

enum PIXEL_TYPE
{
	VIS_IMG_UINT8,
	VIS_IMG_SHORT16,
	VIS_IMG_REAL, //32f
};

enum CHANNEL_INFO
{
	VIS_IMG_GRAY = 1,
	VIS_IMG_RGB = 3,
	VIS_IMG_ARGB = 4
};

enum Threshold_TYPE
{
	CmpLess = 0,
	CmpLessEq = 1,
	CmpEq = 2,
	CmpGreaterEq = 3,
	CmpGreater = 4
};

enum Axis {
	AxsHorizontal = 0,
	AxsVertical = 1,
	AxsBoth = 2,
	Axs45 = 3,
	Axs135 = 4
};

enum Mesurement {
	MesureSum = 0, // IMG_REAL m = MesureSum
	MesureMean = 1,// IMG_REAL m = MesureMean
	MesureStd = 2,// IMG_REAL m = MesureStd
	MesureArea = 3,	// IMG_REAL m = MesureArea
	MesureMin = 4, // IMG_REAL m[3] = {min_value, position_x, position_y}
	MesureMax = 5, // IMG_REAL m[3] = {max_value, position_x, position_y}
	MesureCentroid = 6, // IMG_REAL m[3] = {centermass_value, position_x, position_y}
	MesureCenterpoint = 7,  // IMG_REAL m[3] = {center_value, position_x, position_y}
	MesureAll = 8   // IMG_REAL m[16] = {MesureSum, MesureMean, MesureStd, MesureArea, min_value, position_x, position_y, max_value, position_x, position_y, 
					//centermass_value, position_x, position_y,  center_value, position_x, position_y}
};


enum POSITION {

	TOP_LEFT = 11,
	TOP_CENTER = 12,
	TOP_RIGHT = 13,
	CENTER_LEFT = 21,
	CENTER = 22,
	CENTER_RIGHT = 23,
	BOTTOM_LEFT = 31,
	BOTTOM_CENTER = 32,
	BOTTOM_RIGHT = 33
};

enum COLORS {
	Black = 0,			//0,0,0
	White = 1,			//255,255,255
	Red = 2,				//255,0,0
	Green = 3,			//0,255,0
	Blue = 4,				//0,0,255
	Cyan = 5,			//0,255,255
	Fuchsia = 6,		//255,0,255
	Yellow = 7,			//255,255,0
};

struct DLLEXPORT VIS_IMAGE_INFO
{
	IMG_INT width;
	IMG_INT height;
	IMG_INT byte_per_pixel;
	IMG_INT page_num;
	CHANNEL_INFO channel;
	IMG_INT palette_size = 0;
	IMG_INT sample_per_pixel;
	enum PIXEL_TYPE pixel_type;		//necessary, pixel type, self explain
	enum IMAGE_TYPE img_type;		//necessary, state image channel, self explain
	IMG_INT linestep;

	void clear()
	{
		memset(this, 0x00, sizeof(VIS_IMAGE_INFO));
	}

	VIS_IMAGE_INFO& operator=(const VIS_IMAGE_INFO& o) // rewrite operator=
	{
		width = o.width;
		height = o.height;
		byte_per_pixel = o.byte_per_pixel;
		page_num = o.page_num;
		channel = o.channel;
		palette_size = o.palette_size;
		pixel_type = o.pixel_type;
		img_type = o.img_type;
		sample_per_pixel = o.sample_per_pixel;
		linestep = o.linestep;
		return *this;
	}

	bool operator!=(const VIS_IMAGE_INFO& o) // rewrite operator!=
	{
		bool state = true;
		if (width == o.width && height == o.height && byte_per_pixel == o.byte_per_pixel
			&& page_num == o.page_num && channel == o.channel
			&& palette_size == o.palette_size && pixel_type == o.pixel_type
			&& sample_per_pixel == o.sample_per_pixel && linestep == o.linestep)
		{
			state = false;
		}
		return state;
	}

};

struct VIS_IMG_FILE_INFO
{
	IMG_UBYTE path[MAX_FILE_NAME];
};

struct edgeInformation
{
	IMG_ICOORD xyInteger; //像素点
	IMG_RCOORD xyDecimal;//亚像素点
	int gradient;//灰度梯度
	int grayValue;//灰度值
	float angle;//角度
	bool operator == (const edgeInformation &value)
	{
		return (xyInteger.x == value.xyInteger.x&&
			xyInteger.y == value.xyInteger.y&&
			xyDecimal.x == value.xyDecimal.x&&
			xyDecimal.y == value.xyDecimal.y&&
			gradient == value.gradient&&
			grayValue == value.grayValue&&
			angle == value.angle);
	}
	edgeInformation& operator=(edgeInformation& value)
	{
		xyInteger.x = value.xyInteger.x;
		xyInteger.y = value.xyInteger.y;
		xyDecimal.x = value.xyDecimal.x;
		xyDecimal.y = value.xyDecimal.y;
		gradient = value.gradient;
		grayValue = value.grayValue;
		angle = value.angle;
		return *this;
	}
};//边缘点 
typedef struct _TIFF_FILE_INFORMATION
{
	IMG_RCOORD xyWorld;
	IMG_RSIZE szWorldResolution;
}VIS_TIFF_INFOMATION;

class  DLLEXPORT CVisImage
{
public:
	CVisImage();
	~CVisImage();
	CVisImage(IMG_SIZE szImage, IMG_UINT init_value = 0, CHANNEL_INFO channel = VIS_IMG_GRAY, PIXEL_TYPE type = VIS_IMG_UINT8);
	CVisImage(const CVisImage &img);
	operator IMG_UBBUF * ();
	operator IMG_UBBUF();
	operator IMG_WBUF * ();
	operator IMG_WBUF();
	operator IMG_RBUF * ();
	operator IMG_RBUF();
	//CVisImage(IMG_UBBUF img);			//create by reference
	//CVisImage(IMG_UBBUF img[3]);
	//CVisImage(IMG_WBUF img);
	//CVisImage(IMG_RBUF img);
	CVisImage(const IMG_UBBUF img);
	CVisImage(const IMG_UBBUF img[3]);
	CVisImage(const IMG_WBUF img);
	CVisImage(const IMG_RBUF img);
	static void ExitWriteThread(void);
	IMG_INT Width() { return m_Imginfo.width; }
	IMG_INT Height() { return m_Imginfo.height; }
	IMG_INT GetChannels() { return (IMG_INT)m_Imginfo.channel; }
	IMG_INT GetPageNum() { return m_Imginfo.page_num; }
	PIXEL_TYPE GetPixelType() { return m_Imginfo.pixel_type; }
	bool IsEmpty() { return Width()*Height() == 0; }
	bool GetImageInfo(VIS_IMAGE_INFO &imginfo);
	bool GetImageInfo(VIS_IMAGE_INFO &imginfo) const;
	bool SetImageInfo(VIS_IMAGE_INFO &imginfo);
	bool GetImagePlatte_BMP(IMG_UBYTE *img_platte, IMG_UINT &size);
	bool SetImagePlatte_BMP(IMG_UBYTE *img_platte, IMG_UINT size);
	bool SetImage(void **ppImage);
	bool SetImage(const void **ppImage);
	bool SetImage(const IMG_UBBUF &img);
	bool SetImage(const IMG_UBBUF img[3]);
	bool SetImage(const IMG_WBUF &img);
	bool SetImage(const IMG_RBUF &img);
	bool InitImage(IMG_SIZE szImage, IMG_UINT init_value = 0, CHANNEL_INFO channel = VIS_IMG_GRAY, PIXEL_TYPE type = VIS_IMG_UINT8);
	bool GetImage(const void **ppImage) const;
	bool GetImage(const IMG_UBBUF *pubbuf) const;
	bool GetImage(const IMG_RBUF *pubbuf) const;
	bool GetImage(const IMG_WBUF *pubbuf) const;
	bool GetImage(void **ppImage);
	bool GetImage(IMG_UBBUF *pubbuf);
	bool GetImage(IMG_WBUF *pubbuf);
	bool GetImage(IMG_RBUF *pubbuf);
	bool GetImageByChannel(void *&ppImage, IMG_UINT ch); //Exemple: type BMP/JPG/PNG RGB, channel 1, ch = 1; type tiff 8 bits, Page 2 , ch = 2
	bool ReadImage(const char* filename);
	bool WriteImage(const char *filename);
	bool WriteImage_Syn(const char *filename);
	bool Reset();
	bool CopyFrom(CVisImage &srcImage);
	bool CopyFrom(const CVisImage &srcImage);

	//bool Split(void *pImage);
	bool Merge(void *&pImage);
	bool Add(CVisImage &srcImage);
	bool Add_C(IMG_REAL num);
	bool Sub(CVisImage &srcImage);
	bool Sub_C(IMG_REAL num);
	bool Mul(CVisImage &srcImage);
	bool Mul_C(IMG_REAL num);
	bool Div(CVisImage &srcImage);
	bool Div_C(IMG_REAL num);
	bool Abs(); // for 32 bits
	bool Mean(IMG_LREAL& mean, IMG_INT channel = 0);
	bool MaskMean(CVisImage &maskimg, IMG_LREAL &mean, IMG_INT channel = 0);
	bool Sum(IMG_LREAL& sum, IMG_INT channel = 0);
	bool Min(IMG_REAL& min, IMG_INT channel = 0);
	bool Max(IMG_REAL& max, IMG_INT channel = 0);
	bool Min_Index(IMG_REAL& min, IMG_INT& x, IMG_INT& y, IMG_INT channel = 0);
	bool Max_Index(IMG_REAL& max, IMG_INT& x, IMG_INT& y, IMG_INT channel = 0);
	bool And(CVisImage &srcImage); //byte only
	bool Or(CVisImage &srcImage); //byte only
	bool Xor(CVisImage &srcImage);//byte only
	bool And_C(IMG_BYTE val); //byte only
	bool Or_C(IMG_BYTE val); //byte only
	bool Xor_C(IMG_BYTE val);//byte only
	bool Not();
	bool sqrt();
	bool sqr();
	bool Threshold(IMG_REAL thres, IMG_REAL thres_val, Threshold_TYPE thres_type);
	bool Threshold_double(IMG_REAL LT, IMG_REAL LT_val, IMG_REAL GT, IMG_REAL GT_val);
	bool Convert_8u_to_16s();
	bool Convert_8u_to_32f();
	bool Convert_16s_to_8u();
	bool Convert_16s_to_32f();
	bool Convert_32f_to_8u();
	bool Convert_32f_to_16s();
	bool ResizeLinear(IMG_REAL ScaleSize);
	bool Crop(IMG_COORD start, IMG_SIZE dstsize);
	bool FillPolygonByScanline(IMG_COORD points[], IMG_INT numPoints, void* value);
	bool MesurementPolygonByScanline(IMG_COORD points[], IMG_INT numPoints, IMG_INT option, IMG_REAL* value);
	bool Erosion(IMG_INT KernelSize);
	bool Dilation(IMG_INT KernelSize);
	bool Closing(IMG_INT Width, IMG_INT Height);
	bool Opening(IMG_INT Width, IMG_INT Height);
	bool Gradient(IMG_INT Width, IMG_INT Height);
	bool Mirror(Axis axis);
	bool RTS(IMG_RCOORD &rcoSrcCenter, IMG_RCOORD &rcoDstCenter, IMG_REAL rAngle, IMG_REAL rScale);
	bool RTS(CVisImage &src, IMG_RCOORD &rcoSrcCenter, IMG_RCOORD &rcoDstCenter, IMG_REAL rAngle, IMG_REAL rScale);
	bool FFT(CVisImage &fft_img);
	bool Histogram(IMG_INT nBins, IMG_INT channel, IMG_REAL *pLevels, IMG_UINT* pHistVec);
	bool FFT_Inv(CVisImage &Image_Inv, IMG_UINT Dst_Width, IMG_UINT Dst_Height);
	bool SubImage(IMG_COORD start, IMG_SIZE size, CVisImage &subimage);
	bool FilterGaussian(IMG_INT kernelSize, IMG_REAL sigma);
	bool CenterOfMass(IMG_RCOORD& outputpoint);
	bool Std(IMG_LREAL& std, IMG_LREAL& mean, IMG_INT channel = 0);
	bool Brightness(IMG_REAL BrightLevel = 0);
	bool Equalization(IMG_INT nBins);
	bool FilterPrewittHorizBorder(IMG_INT masksize);
	bool FilterPrewittVertBorder(IMG_INT masksize);
	bool Contrast(IMG_REAL Threshold, IMG_REAL ContrastLevel);
	bool FilterHiPass(IMG_INT kernelSize);
	bool FilterLowPass(IMG_INT kernelSize);
	bool FilterMax(IMG_INT masksize);
	bool FilterMin(IMG_INT masksize);
	bool FilterSobelHorizBorder(IMG_INT masksize);
	bool FilterSobelVertBorder(IMG_INT masksize);
	bool FilterSobelBorder(IMG_INT masksize);
	bool FilterSobelSubpixel(IMG_INT masksize, edgeInformation *edgeArray, IMG_INT &eNum);
	bool FilterMedian(IMG_INT masksize);
	bool FilterMean(IMG_INT masksize);
	bool CannyBorder(IMG_INT normsize, IMG_INT lowThresh, IMG_INT highThresh);
	bool FilterLaplace(IMG_INT kernelSize);
	bool FilterSharpen(IMG_INT kernelSize);
	bool FilterBorder(IMG_WORD *pKernel, IMG_SIZE kernelSize);
	bool Transpose();
	bool CanvasSize(IMG_SIZE roisize, IMG_INT position = 22);
	bool AddRandUniform(IMG_INT low, IMG_INT high, IMG_UINT pSeed);
	bool AddRandGauss(IMG_INT mean, IMG_INT stDev, IMG_UINT pSeed);
	bool CopyConstBorder(IMG_INT borderVal, IMG_INT borderWidth, IMG_INT borderHeight);
	bool CopyMirrorBorder(IMG_INT borderWidth, IMG_INT borderHeight);
	bool MaskImageMesurement(CVisImage &maskImage, Mesurement m, IMG_REAL* value);
	bool ResizeLanczos(IMG_REAL ScaleSize);
	bool ResizeNearest(IMG_REAL ScaleSize);
	bool ResizeCubic(IMG_REAL ScaleSize, IMG_REAL value_B, IMG_REAL value_C);
	bool RGBToGray();
	bool GrayTo16Colors(CVisImage & colors_img);
	bool Colors(IMG_INT R, IMG_INT G, IMG_INT B);
	bool GrayToRGB();
	bool DrawCircle(IMG_CIRCLE Circle, IMG_INT color);
	bool DrawLine(IMG_RCOORD staPoint, IMG_RCOORD endPoint, IMG_INT color,IMG_INT lineWidth = 1);
	bool DrawRect(IMG_WINDOW rectWin, IMG_INT color,IMG_INT lineWidth = 1);
	bool DrawRect(IMG_RORECT rotateRect, IMG_INT color,IMG_INT lineWidth = 1);
	bool DrawWindow(IMG_WINDOW roiWin, IMG_INT color, IMG_INT lineWidth = 1);
	bool DrawPolygon(IMG_RCOORD *Point, IMG_INT color, IMG_INT numPoint);
	bool DrawPoint(IMG_RCOORD *Point, IMG_INT numPoint, IMG_INT color);
	bool DrawArrow(IMG_RCOORD start, IMG_RCOORD end, IMG_INT color);
	bool DrawText_(const char* str, IMG_COORD leftUpCoor, IMG_INT color, int fontHeight, int fontThickness = 400, const char *fn = "Arial", bool italic = false, bool underline = false);
	bool DrawMask(CVisImage &maskimg, IMG_INT color);
	bool Min_Circumcircle(IMG_RCOORD *Point, IMG_INT numPoint, IMG_CIRCLE &Circle);
	bool Min_CircumRorect(IMG_RCOORD *Point, IMG_INT numPoint, IMG_RORECT &RoRect);
	bool ProjectVertical(IMG_VVOID *dst, IMG_INT size_dst, IMG_UWORD channel = 0);	//FOR UINT8 and WORD16, dst is LWORD, for REAL, dst is real
	bool ProjectHorizon(IMG_VVOID *dst, IMG_INT size_dst, IMG_UWORD channel = 0);		//FOR UINT8 and WORD16, dst is LWORD, for REAL, dst is real
	//bool Convolve_2D(CVisImage &dst, IMG_REAL *pKernelX, IMG_INT iKernelSizeX, IMG_REAL *pKernelY, IMG_INT iKernelSizeY);

	bool SetTiffInformation(VIS_TIFF_INFOMATION *tifinfo);
	VIS_TIFF_INFOMATION * GetTiffInformation(void) { return m_pTiffInfo; }
private:
	bool GrayToRGB(CVisImage &colors_img);
	IMG_INT Type(const char* file);
	void ExtractRGB(IMG_UBYTE* pointer, size_t size);
	void WriteData(const int channel);
	bool SetSubImage(void **ppImamge, IMG_UINT sub_start_x, IMG_UINT sub_start_y);
	bool CanvasSizeTemp(IMG_INT srcStep, IMG_INT dstStep, IMG_INT channel, IMG_UBYTE *&pCopyVec, IMG_SIZE roisize, IMG_SIZE StrPoint);
	bool getIntraClassVariance(IMG_UWORD* src, IMG_INT srcRows, IMG_INT srcCols, IMG_INT &varTh, IMG_REAL &varian);

private:
	IMG_UBYTE *m_arrPalette;
	IMG_INT m_paletteSize = 0;
	IMG_INT m_imageType;
	IMG_BOOL m_bInternalMalloc;
	VIS_IMAGE_INFO m_Imginfo;
	char *m_FileName;
	void **m_arrChannel;
	VIS_TIFF_INFOMATION *m_pTiffInfo;
	unsigned char * p_images;
	bool m_isSubImage = false;
public:
	// // free and reset the image class
	void Free();
	void LightCorrectWithRefImage(CVisImage * pRefImage);
	void GenerateLightReferenceImage(void);
};
