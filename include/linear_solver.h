#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "Eigen/IterativeLinearSolvers"

namespace UVLM
{
    /** 求解Ax=b */
    namespace LinearSolver
    {
        template <typename t_a,
                  typename t_b,
                  typename t_x,
                  typename t_options>
        void solve_system
        (
            t_a& a, // 系数矩阵
            t_b& b, // 右端向量
            t_options& options, // 求解器选项
            t_x& x  // 解向量
        )
        {
            if (options.iterative_solver)
            {
                // we use iterative solver 使用迭代法求解线性方程组
                    Eigen::BiCGSTAB<t_a> solver;    // BiCGSTAB稳定双共轭梯度法：求解大型稀疏非对称线性方程组
                    solver.compute(a);
                    x = solver.solveWithGuess(b, x);
            } else
            {
                if(a.rows() == a.cols())
                {
                    //Square Matrix 对方阵使用部分主元变换的LU分解
                    x = a.partialPivLu().solve(b);
                }
                else
                {
                    //Non-square matrix (e.g. phantom cells included in AIC) 对非方阵（例如含幻影面元的AIC）使用全主元变换的LU分解
                    x = a.fullPivLu().solve(b);
                }
            }
        }
    }
}
