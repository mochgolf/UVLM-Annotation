/**
 * @file mapping.h
 * @brief Header file containing functions and data structures related to mapping and transformations in the UVLM framework. 
 *        UVLM框架中与坐标映射和数据类型变换相关的函数和数据结构的头文件
 */

#pragma once

#include "triads.h"
#include "types.h"

/**
 * @namespace UVLM
 * @brief Namespace for the UVLM (Unsteady Vortex Lattice Method) framework. UVLM（非定常涡格法）框架的命名空间
 */
namespace UVLM
{    
    /**
     * @namespace Mapping
     * @brief Namespace for functions and data structures related to mapping and transformations. 映射和变换相关的函数和数据结构的命名空间
     */
    namespace Mapping
    {
       /**
         * @brief Matrix containing the mapping from corner index to matrix indices. 矩阵包含从角点索引到矩阵索引的映射
         *
         * This matrix defines the mapping of corner indices to matrix indices as follows:
         * 这个矩阵定义了角点索引到矩阵索引的映射，如下所示：
         * vortex_indices = [0, 0
         *                   1, 0,
         *                   1, 1,
         *                   0, 1]
         *
         * With the corner numbering as:
         * 角点编号为（逆时针编号）：
         *
         *       N --> 展向-->
         *   0---------3
         *   |         |
         *   |         |
         *   1---------2
         *
         * So, the first element (0) has the coordinates (0, 0). 因此，第一个元素（0）具有坐标（0, 0）。
         * The second element (1) has coordinates (1, 0), and so on. 第二个元素（1）具有坐标（1, 0），依此类推。 
         * 第三个元素（2）具有坐标（1, 1），第四个元素（3）具有坐标（0, 1）。
         */
        const Eigen::Matrix<unsigned int, 4, 2> // 4行2列，如上
                vortex_indices((Eigen::Matrix<unsigned int, 4, 2>()
                                        << 0,0,1,0,1,1,0,1).finished());    // 填入值，如上


        /**
         * @brief Perform bilinear mapping on input and output vectors. 对输入和输出向量执行双线性映射（2x2平均池化）。
         *
         * This function performs bilinear mapping on the input and output vectors.
         * 这个函数对输入和输出向量执行双线性映射，等价于kernel_size=2, stride=1, padding=0的平均池化操作。
         *
         * @tparam t_in Type of the input vector. 输入向量的类型
         * @tparam t_out Type of the output vector. 输出向量的类型
         * @param in Input vector. 输入向量
         * @param out Output vector. 输出向量
         */
        template <typename t_in, typename t_out>
        void BilinearMapping(t_in& in,
                             t_out& out)
        {
            const unsigned int ndims = in.size();
            for (unsigned int idim=0; idim<ndims; ++idim)
            {
                UVLM::Triads::BilinearMap1d(in[idim],
                                            out[idim]); // 对vector中的每个子矩阵指向2x2双线性映射降采样
            }
        }
        /**
         * @brief Map double matrices to a vector of map objects. 将double类型的矩阵（指针数组）映射到一个VecVecMapX对象
         *
         * This function maps double matrices to a vector of map objects.
         * 这个函数将double类型的矩阵（指针数组）映射到一个VecVecMapX对象。
         * 通常用于从Python传入的矩阵（指针数组）到UVLM::Types::VecVecMapX对象，例如`p_zeta`、`p_normals`、`p_dynamic_forces`等。
         *
         * @param dimensions Dimensions of the matrices. 矩阵的维度
         * @param in Input matrices. 输入矩阵
         * @param map Vector of map objects. 映射对象的向量
         * @param correction Correction value for dimensions. 维度的修正值
         * @param n_dim Number of dimensions (default is UVLM::Constants::NDIM). 维度数（默认为UVLM::Constants::NDIM=3）
         */
        void map_VecVecMat(const UVLM::Types::VecDimensions& dimensions,
                           double** in,
                           UVLM::Types::VecVecMapX& map,
                           const int& correction=0,
                           const unsigned int& n_dim=UVLM::Constants::NDIM)
        {
            const unsigned int n_surf = dimensions.size();  // 表面数
            map.resize(n_surf); // 若map的长度小于n_surf，则创建新的默认构造的元素进行填充；否则删除末尾多余的元素
            unsigned int counter = 0;
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    map[i_surf].push_back(UVLM::Types::MapMatrixX (in[counter], // 指针数组的第counter个元素，指向被拉直的x或y或z坐标矩阵
                                                                   dimensions[i_surf].first + correction,
                                                                   dimensions[i_surf].second + correction));
                    counter++;
                }
            }
        }
        /**
         * @brief Map double matrices to a vector of map objects. 将double类型的矩阵（指针数组）映射到一个VecMapX对象
         *
         * This function maps double matrices to a vector of map objects.
         * 这个函数将double类型的矩阵（指针数组）映射到一个VecMapX对象。
         * 通常用于从Python传入的矩阵（指针数组）到UVLM::Types::VecMapX对象，例如`p_gamma`、`p_incidence_angle`等。
         *
         * @param dimensions Dimensions of the matrices.
         * @param in Input matrices.
         * @param map Vector of map objects.
         * @param correction Correction value for dimensions.
         */
        void map_VecMat(const UVLM::Types::VecDimensions& dimensions,
                        double** in,
                        UVLM::Types::VecMapX& map,
                        const int& correction=0)
        {
            const unsigned int n_surf = dimensions.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                map.push_back(UVLM::Types::MapMatrixX (in[i_surf],
                                                       dimensions[i_surf].first + correction,
                                                       dimensions[i_surf].second + correction));
            }
        }
        /**
         * @brief Map double matrices to a vector of vector map objects. 将double类型的矩阵（指针数组）映射到一个VecMapVX对象 std::vector<Eigen::Map<VectorX>>
         *
         * This function maps double matrices to a vector of vector map objects.
         * 这个函数将double类型的矩阵（指针数组）映射到一个VecMapVX对象。
         * 目前未被使用。
         *
         * @param dimensions Dimensions of the matrices. 矩阵的维度
         * @param in Input matrices. 输入矩阵
         * @param map Vector of vector map objects. VecMapVX, std::vector<Eigen::Map<VectorX>>
         * @param correction Correction value for dimensions.
         */
        void map_VecVec1(const UVLM::Types::VecDimensions& dimensions,
                        double** in,
                        UVLM::Types::VecMapVX& map,
                        const int& correction=0)
        {
            // Generates a variable that will be indexed as map[i_surf](i_m)
            const unsigned int n_surf = dimensions.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                map.push_back(UVLM::Types::MapVectorX (in[i_surf],
                                                       dimensions[i_surf].first + correction));
            }
        }
        /**
         * @brief Map a double array to a VectorX object. 将double数组映射到一个VectorX变长向量对象
         *
         * This function maps a double array to a VectorX object.
         * 将double数组映射到一个VectorX变长向量对象，例如`p_rbm_vel_g`
         *
         * @param N_rows Number of rows in the array. 数组中的行数
         * @param in Input array. 输入数组
         * @param out Output VectorX object. 输出VectorX对象
         * @param correction Correction value for dimensions.
        */
        void map_VecX(const uint N_rows,
                        double* in,
                        UVLM::Types::VectorX& out,
                        const int& correction=0)
        {
            // Caution: Use map VecX only for small vectors like p_rbm_vel_g
            // 仅对小向量使用map VecX，例如p_rbm_vel_g
            out.resize(N_rows);
            for (uint i_row=0; i_row < N_rows; ++i_row)
            {
                out[i_row] = in[i_row];   
            }
        }
        /**
         * @brief Transform dimensions from a double array to a vector of dimensions. 将维度从double数组转换为存储维度的std::vector<IntPair>
         *
         * This function transforms dimensions from a double array to a vector of dimensions.
         * 这个函数将维度从double数组转换为存储维度的std::vector<IntPair>。
         *
         * @param n_surf Number of surfaces. 表面数量
         * @param dimensions_in Input dimensions as a double array. 输入维度为double数组
         * @param dimensions Output vector of dimensions. 输出维度的vector
         */
        void transform_dimensions(unsigned int& n_surf,
                                  unsigned int** dimensions_in,
                                  UVLM::Types::VecDimensions& dimensions)
        {
            dimensions.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                dimensions[i_surf].first = dimensions_in[i_surf][0];
                dimensions[i_surf].second= dimensions_in[i_surf][1];
            }
        }

    }
}
