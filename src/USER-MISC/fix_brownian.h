/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(brownian,FixBrownian)

#else

#ifndef LMP_FIX_BROWNIAN_H
#define LMP_FIX_BROWNIAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBrownian : public Fix {
 public:
  FixBrownian(class LAMMPS *, int, char **);
  virtual ~FixBrownian();
  virtual int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void reset_dt();

 protected:
  int seed;
  class AtomVecEllipsoid *avec;
  double dtb;
  double t_start,t_target,t_stop,tsqrt;
  double gamma,Gamma,ascale;
  double dpref,spref;
  boolean do_orientational_dynamics=false;
  class RanMars *random;
  void compute_target();

 private:
  int atoms_have_quaternion();
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
