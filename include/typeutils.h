#pragma once    // 防止头文件被多次包含

#include "EigenInclude.h"
#include "types.h"
#include "triads.h"

#include <cassert>

/**
 * @file typeutils.h
 * @brief Contains overloaded operators for vector-matrix operations. 包含向量-矩阵操作的重载运算符
 */

/**
 * @brief Subtraction operator for VecVecMatrixX types. VecVecMatrixX类型的减法运算符（逐元素减法）
 *
 * Subtracts one VecVecMatrixX from another and returns the result.
 * 从一个VecVecMatrixX中减去另一个VecVecMatrixX并返回结果。
 *
 * @param v1 The first VecVecMatrixX operand. 第一个VecVecMatrixX操作数
 * @param v2 The second VecVecMatrixX operand. 第二个VecVecMatrixX操作数
 * @return The result of the subtraction operation. 减法运算的结果
 */
UVLM::Types::VecVecMatrixX operator-
    (const UVLM::Types::VecVecMatrixX& v1,
     const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);  // 分配与v1相同大小的矩阵

    UVLM::Triads::VecVecMatrix_difference(v1, v2, vout);    // 逐元素相减
    return vout;
}

/**
 * @brief Subtraction operator for VecVecMapX and VecVecMatrixX types. VecVecMapX和VecVecMatrixX类型的减法运算符
 *
 * Subtracts a VecVecMatrixX from a VecVecMapX and returns the result.
 * 从VecVecMapX中减去VecVecMatrixX并返回结果
 *
 * @param v1 The VecVecMapX operand. VecVecMapX操作数
 * @param v2 The VecVecMatrixX operand. VecVecMatrixX操作数
 * @return The result of the subtraction operation. 减法运算的结果
 */
UVLM::Types::VecVecMatrixX operator-
    (const UVLM::Types::VecVecMapX& v1,
     const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_difference(v1, v2, vout);
    return vout;
}

/**
 * @brief Addition operator for VecVecMatrixX types. VecVecMatrixX类型的加法运算符（逐元素加法）
 *
 * Adds two VecVecMatrixX and returns the result.
 * 将两个VecVecMatrixX相加并返回结果。
 *
 * @param v1 The first VecVecMatrixX operand. 第一个VecVecMatrixX操作数
 * @param v2 The second VecVecMatrixX operand. 第二个VecVecMatrixX操作数
 * @return The result of the addition operation. 加法运算的结果
 */
UVLM::Types::VecVecMatrixX operator+
    (const UVLM::Types::VecVecMatrixX& v1,
     const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_addition(v1, v2, vout);
    return vout;
}

/**
 * @brief Compound addition operator for VecVecMatrixX types. VecVecMatrixX类型的复合加法运算符
 *
 * Adds a VecVecMatrixX to the current VecVecMatrixX and modifies the current object.
 * 将VecVecMatrixX添加到当前VecVecMatrixX并修改当前对象。
 *
 * @param own The VecVecMatrixX to be modified. 要修改的VecVecMatrixX
 * @param rhs The VecVecMatrixX to be added. 要添加的VecVecMatrixX
 * @return A reference to the modified VecVecMatrixX. 对修改后的VecVecMatrixX的引用
 */
UVLM::Types::VecVecMatrixX& operator+=
    (UVLM::Types::VecVecMatrixX& own,
     const UVLM::Types::VecVecMatrixX& rhs)
{
    unsigned int n_surf = own.size();
    for (unsigned int i_surf = 0; i_surf < n_surf; ++i_surf)
    {
        unsigned int n_dim = own[i_surf].size();
        for (unsigned int i_dim = 0; i_dim < n_dim; ++i_dim)
        {
            own[i_surf][i_dim] += rhs[i_surf][i_dim];
        }
    }
    return own;
}
