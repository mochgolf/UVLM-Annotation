/**
 * @file types.h
 * @brief Header file containing type definitions and utility functions for the UVLM (Unsteady Vortex Lattice Method) framework.
 */

#pragma once

#include "EigenInclude.h"
#include <vector>
#include <utility>
#include <iostream>

// Convenience declarations
typedef unsigned int uint;

namespace UVLM
{
    namespace Types
    {        
        /**
         * @brief Working precision type.
         */
        typedef double Real; // Real is double

        // Eigen shortcuts for matrices
        // Eigen的默认存储方式是列主序（Column Major），而Numpy和PyTorch默认是行主序（Row Major）
        // 为了正常读取内存中的ctypes数据，Eigen需要使用Eigen::RowMajor
        typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;   // 动态行列尺寸的Eigen矩阵
        typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXint;
        typedef Eigen::Map<MatrixX> MapMatrixX;         // Eigen::Map用于将已经存在的内存块（来自 Python/NumPy 的数据）包装成 Eigen 对象，而不发生内存拷贝
        typedef Eigen::Map<MatrixXint> MapMatrixXint;
        typedef Eigen::Map<const MatrixX> cMapMatrixX;
        typedef std::vector<MatrixX> VecMatrixX;;       // 用std::vector存储Eigen矩阵，实现模拟的三维数组
        typedef std::vector<VecMatrixX> VecVecMatrixX;  // 四维数组
        typedef std::vector<VecVecMatrixX> VecVecVecMatrixX;    // 五维数组
        typedef std::vector<MapMatrixX> VecMapX;        // 用std::vector存储Eigen矩阵的Map对象，实现模拟的Map对象三维数组
        typedef std::vector<MapMatrixXint> VecMapXint;
        typedef std::vector<VecMapX> VecVecMapX;        // Map对象的四维数组
        typedef std::vector<VecVecMapX> VecVecVecMapX;  // Map对象的五维数组
        
        /**
         * @brief Type for dense Eigen matrices. Eigen稠密矩阵的类型
         */
        typedef Eigen::DenseBase<Real> DenseBase;

        /**
         * @brief Type for Eigen matrix block. Eigen矩阵块类型，块是指矩阵的一个矩形子矩阵
         */
        typedef Eigen::Block<MatrixX> Block;    // 动态行列的Eigen矩阵块

        // Eigen shortcuts for vectors
        typedef Eigen::Matrix<Real, 3, 1> Vector3;  // 三维向量
        typedef Eigen::Matrix<Real, 4, 1> Vector4;  // 四维向量
        typedef Eigen::Matrix<Real, 5, 1> Vector5;  // 五维向量
        typedef Eigen::Matrix<Real, 6, 1> Vector6;  // 六维向量
        typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorX; // 动态行尺寸（向量维度）的Eigen向量
        typedef Eigen::Map<VectorX> MapVectorX; // 从内存读取动态维度的Eigen向量
        typedef std::vector<MapVectorX> VecMapVX; // 用std::vector存储Eigen向量的Map对象，实现模拟的Map对象二维数组

        /**
         * @brief Type for integer pairs (e.g., dimensions). 整数对的类型，例如表示维度 dimensions
         */
        typedef std::pair<unsigned int, unsigned int> IntPair;  // 类似Python的元组，但只含有两个元素。存储(行数，列数)或者(M, N)
        /**
         * @brief Type for a vector of integer pairs (e.g., dimensions for multiple surfaces). 包含整数对的向量类型，例如表示多个表面的维度
         */
        typedef std::vector<IntPair> VecDimensions;

        /**
         * @brief Structure representing various options for the VM (Vortex Method) solver. VMopts：代表涡格法VM求解器的各种选项的结构体
         */
        struct VMopts
        {
        	bool ImageMethod; // 对应函数内变量image_method，是否使用镜像法，大部分为False。
        	// unsigned int Mstar;
        	bool Steady;
            bool horseshoe;
        	bool KJMeth;    // 未被SHARPy使用，遗留变量。
        	bool NewAIC;    // 未被SHARPy使用，遗留变量。
        	double DelTime; // SHARPy中始终默认1.0
        	bool Rollup;    // 是否考虑尾流卷起（对流）
            bool only_lifting;
            bool only_nonlifting;
            bool phantom_wing_test; // 标识是否为幻影翼面测试，通常与consider_u_ind_by_sources_for_lifting_forces一起使用。
        	unsigned int NumCores;
        	unsigned int NumSurfaces;
        	unsigned int NumSurfacesNonlifting;
            double dt;
            unsigned int n_rollup;
            double rollup_tolerance;    // SHARPy默认1e-5
            unsigned int rollup_aic_refresh;
            bool iterative_solver;  // 未被SHARPy使用，默认为False。是否使用迭代求解器求解稀疏线性方程组。
            double iterative_tol;   // 未被SHARPy使用，默认为1e-4。
            bool iterative_precond; // 未被SHARPy使用，默认为False。是否使用预条件迭代求解器求解稀疏线性方程组。
            double vortex_radius;
            double vortex_radius_wake_ind;
            bool consider_u_ind_by_sources_for_lifting_forces;  // 是否在计算升力时考虑非升力面（源项）引起的诱导速度，默认False。针对`steady.h`情况。
            uint ignore_first_x_nodes_in_force_calculation;
        };

        /**
         * @brief Structure representing various options for the UVM (Unsteady Vortex Method) solver. 代表UVM求解器的各种选项的结构体
         */
        struct UVMopts
        {
            double dt;
            uint NumCores;
            uint NumSurfaces;
            uint NumSurfacesNonlifting;
            bool only_lifting;
            bool only_nonlifting;
            bool phantom_wing_test;
            uint convection_scheme; // 0, 2, 3。2表示与背景流场`u_ext`对流（常用）
            bool ImageMethod;
            bool iterative_solver;
            double iterative_tol;
            bool iterative_precond;
            bool convect_wake;
            bool cfl1;  // CFL数是否为1
            double vortex_radius;
            double vortex_radius_wake_ind;
            uint interp_coords;
            uint filter_method;
            uint interp_method;
            double yaw_slerp;
            bool quasi_steady;
            uint num_spanwise_panels_wo_induced_velocity;
            bool consider_u_ind_by_sources_for_lifting_forces;
            uint ignore_first_x_nodes_in_force_calculation;
        };
         
         /**
         * @brief This function creates a VMopts object from a UVMopts object. 从UVMopts对象创建VMopts对象
         */
        VMopts UVMopts2VMopts(const UVMopts& uvm)
        {
            VMopts vm;
            vm.dt = uvm.dt;
            vm.NumCores = uvm.NumCores;
            vm.NumSurfaces = uvm.NumSurfaces;
            vm.NumSurfacesNonlifting = uvm.NumSurfacesNonlifting;
            vm.ImageMethod = uvm.ImageMethod;
            vm.iterative_solver = uvm.iterative_solver;
            vm.iterative_tol = uvm.iterative_tol;
            vm.iterative_precond = uvm.iterative_precond;
            vm.horseshoe = false;
            vm.vortex_radius = uvm.vortex_radius;
            vm.vortex_radius_wake_ind = uvm.vortex_radius_wake_ind;
            vm.only_lifting = uvm.only_lifting;
            vm.only_nonlifting = uvm.only_nonlifting;
            vm.Steady = uvm.quasi_steady;
            vm.phantom_wing_test = uvm.phantom_wing_test;
            vm.consider_u_ind_by_sources_for_lifting_forces = uvm.consider_u_ind_by_sources_for_lifting_forces;
            vm.ignore_first_x_nodes_in_force_calculation = uvm.ignore_first_x_nodes_in_force_calculation;
            return vm;  // 返回创建的VMopts对象
        };

         /**
         * @brief Structure representing flight conditions. 飞行条件结构体
         * 
         * Structure includes flow speed and direction, as well as reference chord and density.
         * 该结构体包含来流速度和方向，以及参考弦长和空气密度。
         */
        struct FlightConditions
        {
            double uinf = 1.0;
            double uinf_direction[3];
            double rho = 1.225;
            double c_ref = 1.0;
        };

        /**
         * @brief Correct the dimensions by applying a correction value. 通过应用修正值来校正维度（例如表面的展向和弦向面板数量，M*N)
         *
         * This function corrects dimensions by applying a correction value. It ensures that the dimensions
         * are not negative after the correction.
         * 这个函数通过应用修正值来校正维度。它确保在修正后维度不为负数。
         *
         * @param correction The correction value to be applied. 维度修正值
         * @param dimension The original dimension. 原始维度
         * @return The corrected dimension. 修正后的维度
         */
        inline int correct_dimensions   // 内联函数，让编译器在调用时直接将函数体插入到调用处，避免函数调用的开销
        (
            const int correction,
            int dimension
        )
        {
            if (correction < 0 && dimension < abs(correction))
            {
                return 0;   // 都小于0时返回0
            }
            else
            {
                return dimension + correction;  // 任意一个大于等于0则返回修正维度和原始维度的和。
            }
        }

        /**
         * @brief Generate dimensions from a matrix. 从一个矩阵生成维度
         *
         * This function generates dimensions from a matrix and stores them in a vector of integer pairs.
         * It applies a correction to the dimensions if needed.
         * 这个函数从一个矩阵生成维度，并将它们存储在一个整数对的向量中。如果需要，它会对维度应用修正。
         *
         * @param mat The input matrix. 输入矩阵，只能处理UVLM::Types::VecVecMatrixX和UVLM::Types::VecVecMapX
         * @param dimensions The vector of dimensions to be filled. 要填充的维度向量
         * @param correction The correction value to be applied (default is 0). 维度修正值（默认值为0）
         */
        template <typename t_mat>   // 模板声明，用于泛型编程
        inline void generate_dimensions
        (
            const t_mat& mat,   // `type& name`形式的函数参数意味着引用传递，即传入的参数不会被复制，而是直接使用原始对象的引用（别名）进行in-place操作，以减少内存开销和提高性能
            UVLM::Types::VecDimensions& dimensions,
            const int& correction = 0
        )
        {
            int M, N;
            dimensions.resize(mat.size());  // size()返回矩阵的总元素数。C++在编译时会对模板参数进行检查，例如此处要求传入的t_mat类型必须有`size()`方法。
            
            for (unsigned int i_surf=0; i_surf<dimensions.size(); ++i_surf)
            {
                M = correct_dimensions(correction, mat[i_surf][0].rows());  // correct_dimensions函数接受的dimension参数没有被const修饰，因此会为其创建一个副本并返回
                N = correct_dimensions(correction, mat[i_surf][0].cols());  // 从这里推断出该函数实际上只能处理UVLM::Types::VecVecMatrixX和UVLM::Types::VecVecMapX类型的参数
                dimensions[i_surf] = UVLM::Types::IntPair(M, N);
            }
        }
        
        /**
         * @brief Allocate a vector of matrices with specified dimensions. 分配具有指定维度的矩阵向量：UVLM::Types::VecMatrixX
         *
         * This function allocates a vector of matrices with specified dimensions and initializes them with
         * an initial value.
         * 这个函数分配具有指定维度的矩阵向量，并用初始值初始化它们。
         *
         * @param mat The vector of matrices to be allocated. 要分配的矩阵向量
         * @param dimensions The dimensions for each matrix. 每个矩阵的维度
         * @param correction The correction value to be applied to dimensions (default is 0). 维度修正值（默认值为0）
         * @param initial_value The initial value to fill the matrices (default is 0.0). 矩阵的初始值（默认值为0.0）
         */
        // 分配VecMat
        inline void allocate_VecMat
        (
            UVLM::Types::VecMatrixX& mat,   // equals to std::vector<Eigen::MatrixXd>
            const UVLM::Types::VecDimensions& dimensions,
            const int& correction = 0,
            const UVLM::Types::Real& initial_value = 0.0
        )
        {
            unsigned int n_surf = dimensions.size();    // 从维度中获取表面数
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                int M = dimensions[i_surf].first + correction;  // 获取行数并修正
                int N = dimensions[i_surf].second + correction; // 获取列数并修正
                mat.push_back(UVLM::Types::MatrixX());  // 将一个空矩阵添加到std::vector向量中
                if (initial_value == 0.0)
                {
                    mat[i_surf].setZero(M,N);   // setZero函数是Eigen库的成员函数，用于将矩阵的元素设置为0
                } else if (initial_value == 1.0)
                {
                    mat[i_surf].setOnes(M,N);   // 设置为1
                } else
                {
                    mat[i_surf].setConstant(M,N,initial_value); // 设置为初始值
                }


            }
        }

        // 分配VecMat
        template <typename t_dimensions_in>
        inline void allocate_VecMat // 函数重载
        (
            UVLM::Types::VecMatrixX& mat,
            const t_dimensions_in& dimensions_in,
            const int& correction = 0,
            const double& initial_value = 0.0
        )
        {
            const unsigned int n_mats = dimensions_in.size();
            mat.resize(n_mats);
            for (unsigned int i=0; i<n_mats; ++i)
            {
                mat[i].setConstant
                (
                    dimensions_in[i].rows(),
                    dimensions_in[i].cols(),
                    initial_value

                );
            }
        }

        // 根据参考数组矩阵VecVecMat的维度分配矩阵VecMat
        inline void allocate_VecMat_from_VecVecMat
        (
            UVLM::Types::VecMatrixX& mat,
            const UVLM::Types::VecVecMatrixX& dimensions_in,    // 从这个参考矩阵VecVecMat中获取维度
            const int& correction = 0,
            const double& initial_value = 0.0
        )
        {
            const unsigned int n_mats = dimensions_in.size();
            mat.resize(n_mats);
            for (unsigned int i=0; i<n_mats; ++i)
            {
                if ((dimensions_in[i][0].rows() == 0 )|| (dimensions_in[i][0].cols() == 0))
                {
                    mat[i].resize(0,0); // 如果行数或列数为0，则将矩阵大小设置为0
                }

                else
                {
                    mat[i].setConstant
                    (
                        dimensions_in[i][0].rows() + correction,
                        dimensions_in[i][0].cols() + correction,
                        initial_value
                    );  // 根据参考矩阵维度设置目标矩阵大小
                }
            }
        }

        // 当所有待创建的矩阵都具有相同的维度时使用
        inline void allocate_VecMat
        (
             UVLM::Types::VecMatrixX& mat,
             const unsigned int& n_surf,
             const unsigned int& M,
             const unsigned int& N
        )
        {
             mat.resize(n_surf);
             for (auto& surf: mat)
             {
                 surf.resize(M, N);
                 surf.setZero(M, N);    // 设置为0
             }
        }
        /**
         * @brief Initialize a vector of matrices with a specific value. 用特定值初始化一个矩阵向量 UVLM::Types::VecMatrixX
         *
         * This function initializes a vector of matrices with a specific value.
         * 这个函数用特定值初始化一个矩阵向量。
         *
         * @param mat The vector of matrices to be initialized. 要初始化的矩阵向量
         * @param value The value to initialize the matrices with (default is 0.0). 初始化矩阵的值（默认值为0.0）
         */
        // 初始化VecMat为指定值
        inline void initialise_VecMat
        (
            UVLM::Types::VecMatrixX& mat,
            const double& value = 0.0
        )
        {
            const unsigned int n_mats = mat.size();
            for (unsigned int i=0; i<n_mats; ++i)
            {
                mat[i].setConstant(
                    mat[i].rows(),
                    mat[i].cols(),
                    value);
            }
        }

        // 初始化VecVecMat为指定值
        template <typename t_mat>
        inline void initialise_VecVecMat
        (
            t_mat& mat,
            const double& value = 0.0
        )
        {
            const unsigned int n_surf = mat.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint n_dim = mat[i_surf].size();
                for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    mat[i_surf][i_dim].setConstant(value);
                }
            }
        }
        /**
         * @brief Allocate a vector of vector of matrices with specified dimensions. 分配具有指定维度的矩阵向量向量 UVLM::Types::VecVecMatrixX
         *
         * This function allocates a vector of vector of matrices with specified dimensions and initializes
         * them with an initial value.
         * 这个函数分配具有指定维度的矩阵向量向量，并用初始值初始化它们。
         *
         * @param mat The vector of vector of matrices to be allocated. 要分配的矩阵向量向量
         * @param n_surf The number of surfaces. 表面数量
         * @param n_dim The number of dimensions for each matrix. 每个矩阵的维度数量
         * @param M The number of rows for each matrix. 每个矩阵的行数
         * @param N The number of columns for each matrix. 每个矩阵的列数
         */
        inline void allocate_VecVecMat
        (
            UVLM::Types::VecVecMatrixX& mat,
            const unsigned int& n_surf,
            const unsigned int& n_dim,
            const unsigned int& M,  // 分别传入行数和列数
            const unsigned int& N
        )
        {
            mat.resize(n_surf);
            for (auto& surf: mat)
            {
                surf.resize(n_dim);
                for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    surf.push_back(UVLM::Types::MatrixX());
                    surf[i_dim].resize(M, N);
                    surf[i_dim].setZero(M, N);
                }
            }
        }

        inline void allocate_VecVecMat
        (
            UVLM::Types::VecVecMatrixX& mat,
            const unsigned int& n_dim,
            const UVLM::Types::VecDimensions& dimensions,   // 传入期望的维度整数对Intpair
            const int& correction = 0
        )
        {
            unsigned int n_surf = dimensions.size();
            int M, N;
            mat.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                M = correct_dimensions(correction, dimensions[i_surf].first);
                N = correct_dimensions(correction, dimensions[i_surf].second);
                mat[i_surf].resize(n_dim);
                for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    mat[i_surf].push_back(UVLM::Types::MatrixX());
                    mat[i_surf][i_dim].resize(M, N);
                    mat[i_surf][i_dim].setZero(M, N);
                }
            }
        }

        template <typename t_dimensions>
        inline void allocate_VecVecMat
        (
            UVLM::Types::VecVecMatrixX& mat,
            const t_dimensions& in_dimensions,  // 传入供参考形状的矩阵
            const int& correction = 0
        )
        {
            unsigned int n_surf = in_dimensions.size();
            int M, N;
            mat.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                
                M = correct_dimensions(correction, in_dimensions[i_surf][0].rows());
                N = correct_dimensions(correction, in_dimensions[i_surf][0].cols());  
                mat[i_surf].resize(in_dimensions[i_surf].size());
                for (unsigned int i_dim=0; i_dim<in_dimensions[i_surf].size(); ++i_dim)
                {
                    mat[i_surf][i_dim].resize(M, N);
                    mat[i_surf][i_dim].setZero(M, N);
                }
            }
        }
        /**
         * @brief Copy data from a vector of vector of matrices to another. 从一个矩阵向量向量复制数据到另一个矩阵向量向量
         *
         * This function copies data from one vector of vector of matrices to another.
         * 这个函数将数据从一个矩阵向量向量`UVLM::Types::VecVecMatrixX`复制到另一个。
         *
         * @tparam t_in The input vector of vector of matrices. 输入的矩阵向量向量类型
         * @tparam t_out The output vector of vector of matrices. 输出的矩阵向量向量类型
         * @param in The input vector of vector of matrices. 输入的矩阵向量向量
         * @param out The output vector of vector of matrices. 输出的矩阵向量向量
         */
        template <typename t_in,
                  typename t_out>
        inline void copy_VecVecMat
        (
            const t_in& in,
            t_out& out
        )
        {
            uint M, N;
            uint n_surf = in.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                uint n_dim = in[i_surf].size();
                M = in[i_surf][0].rows();
                N = in[i_surf][0].cols();
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    for (uint i_m=0; i_m<M; ++i_m)
                    {
                        for (uint i_n=0; i_n<N; ++i_n)
                        {
                            out[i_surf][i_dim](i_m, i_n) = in[i_surf][i_dim](i_m, i_n);
                        }
                    }
                }
            }
        }
        /**
         * @brief Copy data from a matrix to a block within another matrix. 从一个矩阵复制数据到另一个矩阵的块中
         *
         * This function copies data from a matrix to a specified block within another matrix.
         * 这个函数将数据从一个矩阵复制到另一个矩阵的指定块中。
         *
         * @tparam mat_in The input matrix. 输入矩阵类型
         * @param in The input matrix to copy from. 要复制的输入矩阵
         * @param out The output matrix (block) to copy to. 要复制到的输出矩阵（块）
         * @param i_start The row index to start copying. 开始复制的行索引
         * @param j_start The column index to start copying. 开始复制的列索引
         */
        template<typename mat_in>
        inline void copy_Mat_to_block
        (
            mat_in& in,
            UVLM::Types::MatrixX& out,
            uint i_start,
            uint j_start
        )
        {
            uint M = in.rows();
            uint N = in.cols();
            for (uint i_m=0; i_m<M; ++i_m)
            {
                for (uint i_n=0; i_n<N; ++i_n)
                {
                    out(i_m+i_start, i_n+j_start) = in(i_m, i_n);
                }
            }
        }
        /**
         * @brief Join two vectors into a single vector. 连接两个向量为一个向量
         *
         * This function concatenates two vectors into a single vector.
         * 这个函数将两个向量连接成一个向量。
         *
         * @param vec1 The first vector. 第一个向量
         * @param vec2 The second vector. 第二个向量
         * @return The concatenated vector. 连接后的向量
         */
        UVLM::Types::VectorX join_vectors
        (
            const UVLM::Types::VectorX& vec1,
            const UVLM::Types::VectorX& vec2
        )
        {
            UVLM::Types::VectorX vec_joined(vec1.size() + vec2.size());
            vec_joined << vec1, vec2;   // <<运算符用于连接两个向量
            return vec_joined; // 返回创建的连接后的新向量
        }

        /**
         * @brief Calculate the norm of a vector of vector of matrices. 计算矩阵向量向量的范数
         *
         * This function calculates the norm of a vector of vector of matrices by summing the norms of individual matrices.
         * 这个函数通过对每个矩阵求范数并求和来计算矩阵向量向量的范数。
         *
         * @tparam t_mat The input vector of vector of matrices. 输入的矩阵向量向量类型
         * @param mat The input vector of vector of matrices. 输入的矩阵向量向量
         * @return The norm of the vector of vector of matrices. 矩阵向量向量的范数
         */
        template <typename t_mat>
        inline double norm_VecVec_mat
        (
            const t_mat& mat
        )
        {
            double norm = 0.0;  // 创建一个变量来存储范数
            uint n_surf = mat.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                uint n_dim = mat[i_surf].size();
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    norm += mat[i_surf][i_dim].norm();  // 累加每个子矩阵的范数
                }
            }
            return norm;
        }        
        /**
         * @brief Calculate the norm of a vector. 计算向量的范数
         *
         * This function calculates the norm of a 3D vector.
         * 这个函数计算3D向量（例如来流速度分量、位置分量）的范数。
         *
         * @param value_1 The first component of the vector. 第一个分量
         * @param value_2 The second component of the vector. 第二个分量
         * @param value_3 The third component of the vector. 第三个分量
         * @return The norm of the vector. 向量的范数
         */

        inline double norm_Vec
        (
            const double value_1,
            const double value_2,
            const double value_3
        )
        {
            
            UVLM::Types::Vector3 vector;
            vector << value_1, value_2, value_3;
            return vector.norm();
        }

        /**
         * @brief Calculate the squared norm of a vector. 计算向量范数的平方
         *
         * This function calculates the squared norm of a 3D vector.
         * 这个函数计算3D向量的范数的平方。
         *
         * @param value_1 The first component of the vector. 第一个分量
         * @param value_2 The second component of the vector. 第二个分量
         * @param value_3 The third component of the vector. 第三个分量
         * @return The squared norm of the vector. 向量的范数的平方
         */
        inline double norm_Vec_squared
        (
            const double value_1,
            const double value_2,
            const double value_3
        )
        {
            
            return value_1 * value_1 + 
                   value_2 * value_2 + 
                   value_3 * value_3;
        }
        /**
         * @brief Calculate the maximum value in a vector of vector of matrices. 计算矩阵向量向量元素的最大值
         *
         * This function calculates the maximum value in a vector of vector of matrices by finding the maximum
         * absolute value in all matrices.
         * 这个函数通过在所有矩阵中找到最大绝对值来计算矩阵向量向量元素的最大值。
         *
         * @tparam t_mat The input vector of vector of matrices. 输入的矩阵向量向量类型
         * @param mat The input vector of vector of matrices. 输入的矩阵向量向量
         * @return The maximum value in the vector of vector of matrices. 矩阵向量向量元素的最大值
         */
        template <typename t_mat>
        inline double max_VecVecMat
        (
            const t_mat& mat
        )
        {
            double max = 0.0;
            uint n_surf = mat.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                uint n_dim = mat[i_surf].size();
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    max = std::max(max, std::abs(mat[i_surf][i_dim].maxCoeff()));
                }
            }
            return max;
        }
        
        // 创建一个分量为0的UVLM::Types::Vector3向量
        UVLM::Types::Vector3 zeroVector3()
        {
            UVLM::Types::Vector3 vec;
            vec.setZero();
            return vec;
        }   

        // 将VecMatrix的子矩阵的值设为Vector3的分量
        template<typename t_vecmatrix>
        inline void pass_3D_Vector_to_VecMatrix
        (
            t_vecmatrix& mat,
            UVLM::Types::Vector3& vec,
            const uint row,
            const uint col

        )
        {
            for (uint i_dim=0; i_dim<3; i_dim++)
            {
                mat[i_dim](row, col) = vec(i_dim);
            }
        }

        // 移除列向量VectorX中的指定一行
        void remove_row_from_VectorX
		(
		UVLM::Types::VectorX & vector_in,
		const int rowToRemove
		)
		{
			int counter_idx = 0;
            UVLM::Types::VectorX intitial_vector_in = vector_in;
            uint vector_initial_size = vector_in.rows();
            vector_in.resize(vector_initial_size-1,1);  // 新向量大小为原向量大小减1
			for (uint i_row = 0;i_row<vector_initial_size;++i_row)
			{
				if (i_row!=rowToRemove)
				{
					vector_in[counter_idx] = intitial_vector_in[i_row];
					counter_idx++;
				}
			}
		}

        // 重新排序VectorX向量的元素，将前numb_of_push_backs行移到最后
        template<typename vec_in>
        UVLM::Types::VectorX reorder_vector_by_pushback
		(
		vec_in& vector_in,
		const int numb_of_push_backs
		)
		{
            uint vector_size = vector_in.rows();            
            UVLM::Types::VectorX vec_out(vector_size);
            uint counter = 0;
			for (uint i_row = numb_of_push_backs;i_row<vector_size;++i_row)
			{
                vec_out[counter] = vector_in[i_row];
                counter++;
			}
            for (uint i_row = 0 ;i_row<numb_of_push_backs;++i_row)
			{
                vec_out[counter] = vector_in[i_row];
                counter++;
			}
            return vec_out;
		}
    }
}



/*
Define types for C++ interface on linear UVLM routines. 定义线性UVLM程序上C++接口的类型。
Matrices size is specified whenever possible to maximise speed. 尽可能指定矩阵大小，以最大限度地提高运行速度。
*/

namespace UVLMlin{

    using namespace Eigen;

    typedef Matrix<double,4,3,RowMajor> Matrix4by3d;

    // map 1d arrays into Eigen Matrices (interface for 1D or 2D python arrays)
    typedef Map< Matrix<double,Dynamic,Dynamic,RowMajor> > map_Mat;
    typedef Map< Matrix<double,4,3,RowMajor> > map_Mat4by3;
    typedef Map< Matrix<double,3,3,RowMajor> > map_Mat3by3;
    typedef Map< Matrix<double,1,3> > map_RowVec3;

    // mapp 3D python arrays
    typedef std::vector<map_Mat3by3> Vec_map_Mat3by3;
    typedef std::vector<map_Mat> Vec_map_Mat;
}



#include "typeutils.h" // 该头文件重载了-, +, +=运算符，用于各类矩阵向量的加减法
