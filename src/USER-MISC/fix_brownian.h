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
  double temp,gamma,ang_gamma;
  double dtpref,drpref,stpref,srpref;
  bool do_orientational_dynamics=false;
  bool do_exact_rotation=false;
  bool is_2d=false;
  class RanMars *random;

 private:
  int atoms_have_quaternion();
  // conversion from rotation vector to quaternion
  inline void rvec_to_quat(const double * r, double * q)
    {
      double norm_r=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
      if(norm_r > 0.0){
        double pref=sin(0.5*norm_r)/norm_r;
        q[0]=cos(0.5*norm_r);
        q[1]=pref*r[0];
        q[2]=pref*r[1];
        q[3]=pref*r[2];
      }
      else{
        q[0]=1.0;
        q[1]=0.0;
        q[2]=0.0;
        q[3]=0.0;
      }
    }
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
