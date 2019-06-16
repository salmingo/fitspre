/*
 * @file ADIProcess.cpp FITS图像处理接口
 */
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "ADIProcess.h"

using namespace std;
using namespace boost::filesystem;
using namespace boost::posix_time;

namespace AstroUtil {
FitsHPtr make_fits_handler() {
	return boost::make_shared<FitsHandler>();
}

FitsNFPtr make_fits_info() {
	FitsNFPtr nfptr = boost::make_shared<FitsInfo>();
	nfptr->hptr = make_fits_handler();
	return nfptr;
}

//////////////////////////////////////////////////////////////////////////////
ADIProcess::ADIProcess() {

}

ADIProcess::~ADIProcess() {
}

bool ADIProcess::CombineZero(const string &pathname, const string &prefix) {
	FitsNFPtrVec fhvec;
	if (!scan_directory(pathname, prefix, fhvec))
		return false;
	int rows, cols, pos, row, col, off1, off2;
	int nfile(fhvec.size()), ifile;
	fltarr rowbuff, pixbuff;

	fhvec[0]->hptr->GetDimension(cols, rows);
	rowbuff.reset(new float[nfile * cols]); // 所有文件的行数据
	pixbuff.reset(new float[nfile]); // 所有文件的点数据
	zero_.reset(new float[rows * cols]); // 处理结果
	info_.valid_zero = false;
	info_.wdim = cols;
	info_.hdim = rows;

	// 合并图像
	for (row = pos = 0; row < rows; ++row) { // 逐行遍历
		// 读出各文件行数据
		for (ifile = 0, off1 = 0; ifile < nfile; ++ifile, off1 += cols) {
			fhvec[ifile]->hptr->LoadPixels(rowbuff.get() + off1, cols, row);
		}
		for (col = 0; col < cols; ++col, ++pos) {
			// 加载各文件(col, row)位置数据
			for (ifile = 0, off1 = col; ifile < nfile; ++ifile, off1 += cols) {
				pixbuff[ifile] = rowbuff[off1];
			}
			zero_[pos] = minmax_clip(pixbuff.get(), nfile);
		}
	}
	// 输出合并结果
	path dst = pathname;
	FitsHPtr fhptr;

	dst /= path("ZERO.fit");
	fhptr = output_image(zero_.get(), cols, rows, dst.string());
	if (!fhptr.unique())
		return false;

	int status(0);
	string dateobs = boost::posix_time::to_iso_extended_string(
			second_clock::universal_time());
	fits_write_key((*fhptr)(), TSTRING, "DATE-OBS", (void*) dateobs.c_str(),
			"time of file genernated", &status);

	info_.valid_zero = !status;
	return info_.valid_zero;
}

bool ADIProcess::SetZero(const string &filepath) {
	FitsHandler fh;
	if (!fh.Open(filepath.c_str()))
		return false;
	int rows, cols, pixels;
	bool same;

	fh.GetDimension(cols, rows);
	if ((pixels = rows * cols) > 0 && pixels != info_.pixels()) {
		zero_.reset(new float[pixels]);
		if ((info_.valid_zero = fh.LoadImage(zero_.get()))) {
			// 检查暗场和平场是否合法
			same = rows == info_.hdim && cols == info_.wdim;
			if (!same) {
				info_.valid_dark = false;
				info_.valid_flat = false;
			}
			// 设置本底状态
			info_.wdim = cols;
			info_.hdim = rows;
		}
	}
	return info_.valid_zero;
}

bool ADIProcess::CombineDark(const string &pathname, const string &prefix) {
	return true;
}

bool ADIProcess::SetDark(const string &filepath) {

	return info_.valid_dark;
}

bool ADIProcess::CombineFlat(const string &pathname, const string &prefix) {
	FitsNFPtrVec fhvec;

	if (!scan_directory(pathname, prefix, fhvec))
		return false;
	int rows, cols, pos, row, col, off1, off2, off3;
	int nfile(fhvec.size()), ifile;
	fltarr rowbuff, pixbuff;
	float expt;

	fhvec[0]->hptr->GetDimension(cols, rows);
	rowbuff.reset(new float[nfile * cols]); // 所有文件的行数据
	pixbuff.reset(new float[nfile]); // 所有文件的点数据
	flat_.reset(new float[rows * cols]); // 处理结果

	// 统计归一化比例尺
	for (ifile = 0; ifile < nfile; ++ifile) {
		fhvec[ifile]->scale = normal_scale(fhvec[ifile]->hptr);
	}

	// 合并图像
	for (row = pos = off3 = 0; row < rows; ++row, off3 += cols) { // 逐行遍历
		// 读出各文件行数据
		for (ifile = 0, off1 = 0; ifile < nfile; ++ifile, off1 += cols) {
			fhvec[ifile]->hptr->LoadPixels(rowbuff.get() + off1, cols, row);
		}
		if (info_.valid_zero) {		// 本底
			float bias;
			// 减本底
			for (col = 0, off2 = off3; col < cols; ++col, ++off2) {
				bias = zero_[off2];
				for (ifile = 0, off1 = col; ifile < nfile;
						++ifile, off1 += cols) {
					rowbuff[off1] -= bias;
				}
			}
		}
		// 合并该行各像素
		for (col = 0; col < cols; ++col, ++pos) {
			// 加载各文件(col, row)位置数据
			for (ifile = 0, off1 = col; ifile < nfile; ++ifile, off1 += cols) {
				pixbuff[ifile] = rowbuff[off1] / fhvec[ifile]->scale;
			}
			flat_[pos] = avsigclip(pixbuff.get(), nfile);
		}
	}
	// 剔除噪声
//	remove_noise(dstbuff.get(), cols, rows); // CMOS相机 效果不明显. 2019-06-12
	// 输出合并结果
	path dst = pathname;
	FitsHPtr fhptr;

	dst /= path("FLAT.fit");
	fhptr = output_image(flat_.get(), cols, rows, dst.string());
	if (!fhptr.unique())
		return false;

	int status(0);
	string dateobs = boost::posix_time::to_iso_extended_string(
			second_clock::universal_time());
	expt = 1.0;
	fits_write_key((*fhptr)(), TSTRING, "DATE-OBS", (void*) dateobs.c_str(),
			"time of file genernated", &status);
	fits_write_key((*fhptr)(), TFLOAT, "EXPTIME", &expt, "Exposure duration",
			&status);

	if (!status)
		info_.valid_flat = true;

	return false;
}

bool ADIProcess::SetFlat(const string &filepath) {

	return info_.valid_flat;
}

void ADIProcess::Reset(int type) {
	if (type == 0) { // 本底
		info_.valid_zero = false;
		zero_.reset();
	} else if (type == 1) { // 暗场
		info_.valid_dark = false;
		dark_.reset();
	} else { // type == 2, 平场
		info_.valid_flat = false;
		flat_.reset();
	}
}

bool ADIProcess::scan_directory(const string &pathname, const string &prefix,
		FitsNFPtrVec &vec) {
	directory_iterator itend = directory_iterator();
	string filename;
	int rows1, cols1, rows2, cols2;

	// 遍历并打开文件
	for (directory_iterator x = directory_iterator(pathname); x != itend; ++x) {
		filename = x->path().filename().string();
		if (filename.find(prefix))
			continue;

		FitsNFPtr fnfptr = make_fits_info();
		if (!fnfptr->hptr->Open(x->path().c_str()))
			continue;
		vec.push_back(fnfptr);
	}
	if (vec.size() < 3)
		return false;
	// 检查图像一致性
	vec[0]->hptr->GetDimension(cols1, rows1);
	for (FitsNFPtrVec::iterator it = vec.begin() + 1; it != vec.end();) {
		(*it)->hptr->GetDimension(cols2, rows2);
		if (rows1 != rows2 || cols1 != cols2)
			it = vec.erase(it);
		else
			++it;
	}
	return vec.size() >= 3;
}

FitsHPtr ADIProcess::output_image(float *data, int cols, int rows,
		const string &pathname) {
	FitsHPtr fhptr = make_fits_handler();
	FitsHPtr fhret;

	if (exists(pathname))
		remove(pathname);
	if (fhptr->CreateImage(pathname.c_str(), FLOAT_IMG, cols, rows)
			&& fhptr->WriteImage(data, TFLOAT)) {
		fhret = fhptr;
	}
	return fhret;
}

/*
 * 基于样本中值, 对图像数据做归一化处理
 */
float ADIProcess::normal_scale(FitsHPtr fhptr) {
	int rows, cols, pixels;
	fltarr data;
	fhptr->GetDimension(cols, rows);
	pixels = rows * cols;
	data.reset(new float[pixels]);
	fhptr->LoadImage(data.get());
	return normal_scale(data.get(), pixels);
}

float ADIProcess::normal_scale(float *data, int n) {
	int ns = n > 10000 ? 10000 : n;
	int off = ns == n ? 0 : (n % ns) / 2;
	int i;
	float scale, pos(0.0);
	float step = ns == n ? 1.0 : float((n * 1.0 / ns));
	vector<float> tmp;

	for (i = 0; off < n; ++i) {
		tmp.push_back(data[off]);
		pos += step;
		off = int(pos);
	}
	nth_element(tmp.begin(), tmp.begin() + ns / 2, tmp.end());
	scale = tmp[ns / 2];
	tmp.clear();
	return scale;
}

void ADIProcess::conv_filter(float *x, int w, int h) {
	int whalf(w / 2), hhalf(h / 2);
}

float ADIProcess::minmax_clip(float *x, int n) {
	if (n < 3)
		return 0.0;

	float min(1E30), max(-1E30);
	double sum;
	for (int i = 0; i < n; ++i) {
		if (x[i] < min)
			min = x[i];
		if (x[i] > max)
			max = x[i];
		sum += x[i];
	}
	return (sum - min - max) / (n - 2);
}

float ADIProcess::avsigclip(float *x, int n, float lsigma, float hsigma) {
	double sum, sq;
	float min(1E30), max(-1E30), mean, rms, low, high, t;
	int i, n1, n2;

	sum = sq = 0.0;
	for (i = 0; i < n; ++i) {
		sum += (t = x[i]);
		sq += (x[i] * x[i]);
		if (min > t)
			min = t;
		if (max < t)
			max = t;
	}

	n2 = n;
	do {
		sq -= (min * min + max * max);
		mean = float((sum - min - max) / (n2 - 2));
		rms = float(sqrt((sq - (sum - min - max) * mean) / (n2 - 3)));
		low = mean - lsigma * rms;
		high = mean + hsigma * rms;
		n1 = n2;
		sum = sq = 0.0;
		min = 1E30;
		max = -1E30;
		for (i = 0, n2 = 0; i < n; ++i) {
			if ((low < (t = x[i])) && t < high) {
				sum += t;
				sq += (t * t);
				if (min > t)
					min = t;
				if (max < t)
					max = t;
				++n2;
			}
		}
	} while (n2 > 3 && n1 > n2);
	return n2 > 3 ? ((sum - min - max) / (n2 - 2)) : mean;
}

void ADIProcess::remove_noise(float *x, int cols, int rows) {
	int pixels = cols * rows, i;
	int row, col;
	float median, mean, rms, low, high, min(1E30), max(-1E30), t;
	double sum(0.0), sq(0.0);
	fltarr tmp;
	// 备份原始数据
	tmp.reset(new float[pixels]);
	memcpy(tmp.get(), x, pixels * sizeof(float));
	// 统计
	median = normal_scale(x, pixels);
	for (i = 0; i < pixels; ++i) {
		sum += (t = x[i]);
		sq += (t * t);
		if (min > t)
			min = t;
		if (max < t)
			max = t;
	}
	sum -= (min + max);
	sq -= (min * min + max * max);
	mean = float(sum / (pixels - 2));
	rms = float(sqrt((sq - sum * mean) / (pixels - 3)));
	low = median - 3.0 * rms;
	high = median + 3.0 * rms;
	// 遍历剔除噪声
	for (row = i = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col, ++i) {
			if (low > (t = x[i]) || high < t) {
				x[i] = av_replace(tmp.get(), col, row, cols, rows);
			}
		}
	}
}

float ADIProcess::av_replace(float *x, int col, int row, int cols, int rows) {
	int w(2);
	int r1 = row - w;
	int r2 = row + w;
	int c1 = col - w;
	int c2 = col + w;
	int n, r, c, d, offset;
	float min(1E30), max(-1E30), pix = x[row * cols + col], t;
	double sum(0.0);

	if (r1 < 0)
		r1 = 0;
	if (r2 >= rows)
		r2 = rows - 1;
	if (c1 < 0)
		c1 = 0;
	if (c2 >= cols)
		c2 = cols - 1;

	n = (r2 - r1 + 1) * (d = c2 - c1 + 1);
	for (r = r1, offset = r1 * cols + c1; r <= r2; ++r) {
		for (c = c1; c <= c2; ++c, ++offset) {
			sum += (t = x[offset]);
			if (min > t)
				min = t;
			if (max < t)
				max = t;
		}
		offset += (cols - d);
	}
	return (sum - min - max - pix) / (n - 3);
}

bool ADIProcess::load_badpixel(const string &filepath) {
	return false;
}

void ADIProcess::pre_process() {

}

void ADIProcess::do_process() {

}
//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */
