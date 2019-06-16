/*
 * @file ADIProcess.h FITS图像处理接口
 * @version 0.1
 * @author Xiaomeng
 * @note
 * - 图像合并: 本底, 暗场, 平场. 合并图像以float型格式存储
 * - 原始文件中, EXPTIME对应曝光时间
 */

#ifndef ADIPROCESS_H_
#define ADIPROCESS_H_

#include <boost/smart_ptr.hpp>
#include <boost/container/stable_vector.hpp>
#include <string>
#include "FitsHandler.h"

using std::string;

namespace AstroUtil {
// 声明数据类型
typedef boost::shared_array<float> fltarr;
typedef boost::container::stable_vector<float> fltvec;
typedef boost::shared_ptr<FitsHandler> FitsHPtr;
//typedef boost::container::stable_vector<FitsHPtr> FitsHPtrVec;

struct FitsInfo { // fits文件信息
	FitsHPtr hptr;	//< 访问指针
	float scale;	//< 归一化比例尺
};
typedef boost::shared_ptr<FitsInfo> FitsNFPtr;
typedef boost::container::stable_vector<FitsNFPtr> FitsNFPtrVec;

//////////////////////////////////////////////////////////////////////////////
struct param_dip {	//< 图像处理及信号提取参数
	int bkw, bkh;		//< 背景拟合窗口
	int bkfrw, bkfh;	//< 背景拟合滤波窗口
	int minarea;		//< 最小连通域面积
};

struct info_adip {
	bool valid_zero;	//< valid ZERO flag
	bool valid_dark;	//< valid DARK flag
	bool valid_flat;	//< valid FLAT flag
	int wdim, hdim;		//< image dimension

public:
	int pixels() {
		return wdim * hdim;
	}

	bool same_dimension(int w, int h) {
		return (w == wdim && h == hdim);
	}
};

class ADIProcess {
public:
	ADIProcess();
	virtual ~ADIProcess();

protected:
	info_adip info_;	//< 图像信息
	fltarr zero_;	//< 本底数据
	fltarr dark_;	//< 暗场数据
	fltarr flat_;	//< 平场数据
	fltarr back_;	//< 图像背景
	fltarr rms_;	//< 图像噪声

public:
	/*!
	 * @brief 合并本底
	 * @param pathname 文件存储路径
	 * @param prefix   文件名前缀
	 * @return
	 * 本底合并结果
	 */
	bool CombineZero(const string &pathname, const string &prefix);
	/*!
	 * @brief 设置合并后本底路径
	 * @param filepath 文件路径
	 * @return
	 * 本底加载结果
	 */
	bool SetZero(const string &filepath);
	/*!
	 * @brief 合并暗场
	 * @param pathname 文件存储路径
	 * @param prefix   文件名前缀
	 * @return
	 * 本底合并结果
	 */
	bool CombineDark(const string &pathname, const string &prefix);
	/*!
	 * @brief 设置合并后暗场路径
	 * @param filepath 文件路径
	 * @return
	 * 本底加载结果
	 */
	bool SetDark(const string &filepath);
	/*!
	 * @brief 合并平场
	 * @param pathname 文件存储路径
	 * @param prefix   文件名前缀
	 * @return
	 * 本底合并结果
	 */
	bool CombineFlat(const string &pathname, const string &prefix);
	/*!
	 * @brief 设置合并后平场路径
	 * @param filepath 文件路径
	 * @return
	 * 本底加载结果
	 */
	bool SetFlat(const string &filepath);
	/*!
	 * @brief 重置标定用图像
	 * @param type 图像类型. 0: 本底; 1: 暗场; 2: 平场
	 * @param
	 */
	void Reset(int type = 0);

protected:
	/*!
	 * @brief 扫描目录, 查找符合条件的文件
	 * @param pathname 目录名
	 * @param prefix   文件名前缀
	 * @param vec      符合条件的FITS文件操作句柄
	 * @return
	 */
	bool scan_directory(const string &pathname, const string &prefix,
			FitsNFPtrVec &vec);
	/*!
	 * @brief 输出图像为FLOAT型FITS文件
	 * @param pathname 文件路径
	 * @param cols     图像宽度
	 * @param rows     图像高度
	 * @return
	 * FITS文件指针
	 */
	FitsHPtr output_image(float *data, int cols, int rows,
			const string &pathname);
	/*!
	 * @brief 基于样本, 计算图像数据归一化比例尺
	 */
	float normal_scale(FitsHPtr fhptr);
	float normal_scale(float *data, int n);
	/*!
	 * @brief 卷积滤波
	 * @param x 卷积核
	 * @param w 卷积宽度
	 * @param h 卷积高度
	 */
	void conv_filter(float *x, int w, int h);
	/*!
	 * @brief 使用min-max计算均值
	 * @param x  待统计数据
	 * @param n  数据长度
	 * @return
	 * 统计结果
	 */
	float minmax_clip(float *x, int n);
	/*!
	 * @brief 基于信噪比的筛选统计
	 * @param x      待统计数据
	 * @param n      数据长度
	 * @param lsigma 下限信噪比
	 * @param hsigma 上限信噪比
	 * @return
	 * 统计结果
	 */
	float avsigclip(float *x, int n, float lsigma = 3.0, float hsigma = 3.0);
	/*!
	 * @brief 使用坏像素周边5*5范围内其它数据替代该值
	 * @param x    数据存储区
	 * @param cols 列数
	 * @param rows 行数
	 */
	void remove_noise(float *x, int cols, int rows);
	/*!
	 * @brief 计算临近5*5的平均值
	 */
	float av_replace(float *x, int col, int row, int cols, int rows);
	/*!
	 * @brief 加载坏像素
	 * @param filepath 坏像素记录文件
	 */
	bool load_badpixel(const string &filepath);
	/*!
	 * @brief 图像预处理
	 * @note
	 * - 减本底
	 * - 减暗场
	 * - 除平场
	 */
	void pre_process();
	/*!
	 * @brief 处理图像. 提取图像中目标
	 */
	void do_process();
};
//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */

#endif /* ADIPROCESS_H_ */
