/**
 * @file constants.h
 * @brief Header file containing definitions of constants used in the UVLM program. UVLM程序中使用的常量定义
 */

#pragma once

#include "types.h"

namespace UVLM
{
    namespace Constants
    {
        /**
         * @brief Constants used in the UVLM program. UVLM程序中使用的常量
         */
        
        /// The number of dimensions in the program (always set to 3 for 3D space).
        /// 程序中的维度数量（始终设置为3，表示三维空间）
        const unsigned int NDIM = 3;

        /// Mathematical constant pi (π) with high precision.
        /// 数学常数π（圆周率），高精度
        const UVLM::Types::Real PI = 3.1415926535897932384626433832795028841971;

        /// Four times the value of pi (4π).
        /// 四倍的π（4π）
        const UVLM::Types::Real PI4 = 4.0 * PI;

        /// The reciprocal of four times pi (1 / 4π).
        /// 四倍π的倒数（1 / 4π）
        const UVLM::Types::Real INV_PI4 = 1.0 / PI4;

        /// Conversion factor for degrees to radians (π / 180.0).
        /// 角度转弧度的转换因子（π / 180.0）
        const UVLM::Types::Real DEGREES2RAD = PI / 180.0;

        /// A small positive value used as a numerical tolerance (10 times machine epsilon).
        /// 用作数值容差的小正值（10倍机器精度）
        const UVLM::Types::Real EPSILON = 10 * std::numeric_limits<UVLM::Types::Real>::epsilon();
    }
}
