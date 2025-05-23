#pragma once

#include "types.h"

#include <iostream>

/**
 * @file triads.h
 * @brief Contains utility functions for matrix operations and transformations. 包含矩阵操作和变换的实用函数
 */

namespace UVLM
{
    namespace Triads    // UVLM::Triads命名空间
    {
        /**
         * @brief Perform a bilinear mapping on a matrix. 在一个矩阵上执行双线性映射（类似平均池化，下采样）
         *
         * This function calculates the mean of each 2x2 block of a matrix and stores the result in another matrix.
         * 这个函数计算矩阵的每个2x2块的平均值，并将结果存储在另一个矩阵中。2x2矩阵块-->单个元素
         *
         * @tparam t_mat The input matrix type. 输入矩阵类型
         * @tparam t_out The output matrix type. 输出矩阵类型
         * @param mat The input matrix to be mapped. 输入矩阵
         * @param out The output matrix to store the mapped values. 输出矩阵，用于存储映射值
         */
        template <typename t_mat, typename t_out>
        inline void BilinearMap1d(const t_mat& mat, t_out& out)
        {
            const int n_rows = out.rows();
            const int n_cols = out.cols();

            for (int i = 0; i < n_rows; ++i)
            {
                for (int j = 0; j < n_cols; ++j)
                {
                    out(i, j) = mat.template block<2, 2>(i, j).mean();  // 计算mat(i,j)的2x2块的平均值，并存入out(i,j)
                }
            }
        }

        /**
         * @brief Perform an inverse bilinear mapping on a matrix. 在一个矩阵上执行逆双线性映射（上采样）
         *
         * This function transfers information from the collocation points to the corner points in the matrix.
         * The input matrix has one less element per dimension than the output matrix.
         * 这个函数将信息从矩阵中的配点传递到角点。输入矩阵在每个维度上比输出矩阵少一个元素。
         *
         * @tparam t_mat The input matrix type. 输入矩阵类型
         * @tparam t_out The output matrix type. 输出矩阵类型
         * @param mat The input matrix to be mapped. 被映射的输入矩阵
         * @param out The output matrix to store the mapped values. 输出矩阵，用于存储映射值
         */
        template <typename t_mat, typename t_out>
        inline void InvBilinearMap1d(const t_mat& mat, t_out& out)
        {
            const int n_rows = mat.rows();
            const int n_cols = mat.cols();

            for (int i = 0; i < n_rows; ++i)
            {
                for (int j = 0; j < n_cols; ++j)
                {
                    out.template block<2, 2>(i, j) += 0.25 * mat(i, j);
                }
            }
        }

        /**
         * @brief Calculate the difference between two matrices element-wise. 逐元素计算两个矩阵的差
         *
         * @tparam t_1 The first matrix type. 第一个矩阵类型
         * @tparam t_2 The second matrix type. 第二个矩阵类型
         * @tparam t_out The output matrix type. 输出矩阵类型
         * @param mat1 The first input matrix. 第一个输入矩阵
         * @param mat2 The second input matrix. 第二个输入矩阵
         * @param mat_out The output matrix to store the element-wise difference. 输出矩阵，用于存储逐元素差异
         */
        template <typename t_1, typename t_2, typename t_out>
        void VecVecMatrix_difference(const t_1& mat1, const t_2& mat2, t_out& mat_out)
        {
            for (unsigned int i_surf = 0; i_surf < mat1.size(); ++i_surf)
            {
                for (unsigned int i_dim = 0; i_dim < mat1[i_surf].size(); ++i_dim)
                {
                    mat_out[i_surf][i_dim].noalias() = mat1[i_surf][i_dim] - mat2[i_surf][i_dim];
                }
            }
        }

        /**
         * @brief Calculate the element-wise addition of two matrices. 逐元素计算两个矩阵的和
         *
         * @tparam t_1 The first matrix type. 第一个矩阵类型
         * @tparam t_2 The second matrix type. 第二个矩阵类型
         * @tparam t_out The output matrix type. 输出矩阵类型
         * @param mat1 The first input matrix. 第一个输入矩阵
         * @param mat2 The second input matrix. 第二个输入矩阵
         * @param mat_out The output matrix to store the element-wise sum. 输出矩阵，用于存储逐元素和
         */
        template <typename t_1, typename t_2, typename t_out>
        void VecVecMatrix_addition(const t_1& mat1, const t_2& mat2, t_out& mat_out)
        {
            for (unsigned int i_surf = 0; i_surf < mat1.size(); ++i_surf)
            {
                for (unsigned int i_dim = 0; i_dim < mat1[i_surf].size(); ++i_dim)
                {
                    mat_out[i_surf][i_dim].noalias() = mat1[i_surf][i_dim] + mat2[i_surf][i_dim];   // .noalias()不创建临时对象存储中间结果，用于避免不必要的内存分配
                }
            }
        }

        /**
         * @brief Add two matrices element-wise and update the first matrix. 逐元素相加两个矩阵并更新第一个矩阵
         *
         * @tparam t_1 The first matrix type. 第一个矩阵类型
         * @tparam t_2 The second matrix type. 第二个矩阵类型
         * @param mat_in_and_out The input matrix to be updated. 要更新的输入矩阵
         * @param mat_in The second input matrix to add. 第二个输入矩阵
         */
        template <typename t_1, typename t_2>
        void VecVecMatrix_addition(t_1& mat_in_and_out, const t_2& mat_in)
        {
            for (uint i_surf = 0; i_surf < mat_in_and_out.size(); ++i_surf)
            {
                for (uint i_dim = 0; i_dim < mat_in_and_out[i_surf].size(); ++i_dim)
                {
                    mat_in_and_out[i_surf][i_dim] += mat_in[i_surf][i_dim];
                }
            }
        }

        /**
         * @brief Calculate the element-wise difference between two matrices and update the first matrix. 逐元素计算两个矩阵的差并更新第一个矩阵
         *
         * @tparam t_1 The first matrix type. 第一个矩阵类型
         * @tparam t_2 The second matrix type. 第二个矩阵类型
         * @param mat_in_and_out The input matrix to be updated. 要更新的输入矩阵
         * @param mat_in The second input matrix for subtraction. 第二个输入矩阵
         */
        template <typename t_1, typename t_2>
        void VecVecMatrix_difference(t_1& mat_in_and_out, const t_2& mat_in)
        {
            for (uint i_surf = 0; i_surf < mat_in_and_out.size(); ++i_surf)
            {
                for (uint i_dim = 0; i_dim < mat_in_and_out[i_surf].size(); ++i_dim)
                {
                    mat_in_and_out[i_surf][i_dim] -= mat_in[i_surf][i_dim];
                }
            }
        }
    }
}
