/*
 Name        : fitspre.cpp
 Author      : Xiaomeng Lu
 Copyright   : SVOM@NAOC, CAS
 Description : 预处理FITS图像
 @version 0.1
 @author Xiaomeng Lu
 @note
 - 1: 合并本底
 - 2: 合并暗场
 - 3: 合并平场
 - 4: 图像预处理: 减本底/减暗场/除平场
 @note
 - 合并后本底存储为 : ZERO.fit
 - 合并后暗场存储为 : DARK.fit
 - 合并后平场存储为 : FLAT.fit
 - 原始图像文件目录与处理结果存储路径必须不同, 处理后图像文件将以原名存至结果路径
 - 无本底不可合并平场
 @note
 - 本底和暗场合并采用 min-max算法
 - 平场合并采用 av-sigclip算法. sigma <= 2
 -
 @note
 配置文件名: fitspre.xml
 配置文件搜索路径:
 - 当前目录
 - /usr/local/etc
 配置文件参数项:
 - BIAS_DIR      : 本底文件目录
 - BIAS_PREFIX   : 本底文件前缀
 - DARK_DIR      : 暗场文件目录
 - DARK_PREFIX   : 暗场文件前缀
 - FLAT_DIR      : 平场文件目录
 - FLAT_PREFIX   : 平场文件前缀
 - IMAGE_DIR     : 图像文件目录
 - IMAGE_PREFIX  : 图像文件名前缀
 - RESULT_DIR    : 处理结果文件存储目录
 */

#include <stdio.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include "ADIProcess.h"

using namespace std;
using namespace boost::filesystem;
using namespace AstroUtil;

void print_help() {

}

/*
 * 命令行参数:
 * -m 处理模式. 0: 合并本底; 1: 合并暗场; 2: 合并平场; 3: 处理图像/提取目标. 缺省值3
 * -i 原文件目录
 * -p 原文件名前缀. 缺省时处理原文件目录下所有文件
 * -o 结果存储目录
 */
int main(int argc, char **argv) {
	// 解析命令行参数


	// 图像处理
	// 输出处理结果

	ADIProcess adip;
	string pathname = "/Users/lxm/Data/processing";
	string prefix = "Geo";

	adip.SetZero("/Users/lxm/Data/processing/ZERO.fit");
	if (!adip.CombineFlat(pathname, prefix))
		printf("failed\n");
	else
		printf("succeed\n");

//	string pathname = "/Volumes/NETAC/bias";
//	string prefix = "bias_";
//	if (!adip.CombineZero(pathname, prefix))
//		cout << "failed" << endl;
//	else
//		cout << "success" << endl;

	return 0;
}
