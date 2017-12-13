/*

This header and the corresponding source file is external LGPL-3 source code, taken
from the alta rev. 598de8d554f repository as instructed by Gael GUENNEBAUD.
It was itself adapted from the QuadProg++ package written by Luca Di Gaspero
and Eric Moyer, as described below. The only modifications are a renaming of the files
and this header.

The quadprog_solve() function implements the algorithm of Goldfarb and Idnani
for the solution of a (convex) Quadratic Programming problem
by means of an active-set dual method.

The problem is in the form:

min 0.5 * x G x + g0 x
s.t.
    CE^T x + ce0 = 0
    CI^T x + ci0 >= 0

The matrix and vectors dimensions are as follows:
    G: n * n
    g0: n

    CE: n * p
  ce0: p

    CI: n * m
  ci0: m

    x: n

The function will return the cost of the solution written in the x vector or
std::numeric_limits::infinity() if the problem is infeasible. In the latter case
the value of the x vector is not correct.

References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
            strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

Notes:
  1. pay attention in setting up the vectors ce0 and ci0.
    If the constraints of your problem are specified in the form
    A^T x = b and C^T x >= d, then you should set ce0 = -b and ci0 = -d.
  2. The matrix G is modified within the function since it is used to compute
    the G = L^T L cholesky factorization for further computations inside the function.
    If you need the original matrix G you should make a copy of it and pass the copy
    to the function.

Author: Luca Di Gaspero
        DIEGM - University of Udine, Italy
        l.digaspero@uniud.it
        http://www.diegm.uniud.it/digaspero/

The author will be grateful if the researchers using this software will
acknowledge the contribution of this function in their research papers.

LICENSE

This file is part of QuadProg++: a C++ library implementing
the algorithm of Goldfarb and Idnani for the solution of a (convex)
Quadratic Programming problem by means of an active-set dual method.
Copyright (C) 2007-2009 Luca Di Gaspero.
Copyright (C) 2009 Eric Moyer.

QuadProg++ is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

QuadProg++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with QuadProg++. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _QUADPROGPP
#define _QUADPROGPP

#include <Eigen/Core>

namespace QuadProgPP{
  using namespace Eigen;

  struct Scheduling
  {
    enum ModeType {
      WorstFirst,     // Default strategy: check all inequalities, pick the worst one (if any), solve it
      WorstSetFirst,  // chek all inequalities, and pick about "initial_chunk_size" worst ones, and then satisfy them all using the default strategy
      SlidingWindows, // check the current chunk of initial_chunk_size inequalities, pick and solve the worst one (if any), shift the next chunk and repeat until all chunks are satisfied
    };

    ModeType mode;
    int initial_chunk_size;
    double grow_factor;

    Scheduling()
      : mode(WorstFirst), initial_chunk_size(-1), grow_factor(1)
    {}

    void initSlidingWindows(int n, int m)
    {
      mode = SlidingWindows;
      initial_chunk_size = std::min(std::max(n/2,m/64+1),m);
      grow_factor = 1;
    }

    void initWorstSetFirst(int n, int m)
    {
      mode = WorstSetFirst;
      initial_chunk_size = std::min(2*n,m);
      grow_factor = 1;
    }

    static void shift_window(int &chunk_start, int &chunk_size, int m)
    {
      if(chunk_start+chunk_size==m)
        chunk_start = 0;
      else
      {
        chunk_start += chunk_size;
        if(chunk_start+chunk_size>m)
          chunk_start = m-chunk_size;
      }
    }
  };

  void init_qp(Ref<MatrixXd> G);

  // Shortcut for init_qp(G) followed by a solve of the quadratic energy, and a call to solve_quadprog_with_guess().
  double solve_quadprog(Ref<MatrixXd> G, Ref<const VectorXd> g0,
                        Ref<const MatrixXd> CE, Ref<const VectorXd> ce0,
                        Ref<const MatrixXd> CI, Ref<const VectorXd> ci0,
                        Ref<VectorXd> x);

  /**
  * \param L Cholesky factor of the quadratic objective as computed by init_qp
  * \param g0 linear part of the quadratic objective
  * \param x on input: initial solution, on output the optimal solution (if any)
  * \param scheduling strategy used to prioritize inequalities
  * \param active_set if not null, treats the indexed inequality first, and on output active_set is filled with the remaining active-set
  */
  double solve_quadprog_with_guess(Ref<const MatrixXd> L, Ref<const VectorXd> g0,
                                   Ref<const MatrixXd> CE, Ref<const VectorXd> ce0,
                                   Ref<const MatrixXd> CI, Ref<const VectorXd> ci0,
                                   Ref<VectorXd> x,
                                   Scheduling scheduling = Scheduling(),
                                   VectorXi *active_set = 0);
}

#endif // #define _QUADPROGPP
