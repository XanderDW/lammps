/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Xander de Wit (Eindhoven University of Technology, Netherlands)
------------------------------------------------------------------------- */

#include "fix_brownian.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "math_extra.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;


/* ---------------------------------------------------------------------- */

FixBrownian::FixBrownian(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if ((narg != 7) && (narg != 9)) error->all(FLERR,"Illegal fix brownian command");

  t_start = utils::numeric(FLERR,arg[3],false,lmp);
  t_target = t_start;
  t_stop = utils::numeric(FLERR,arg[4],false,lmp);
  gamma = utils::numeric(FLERR,arg[5],false,lmp);
  if (gamma <= 0.0) error->all(FLERR,"Fix brownian gamma must be > 0.0");
  seed = utils::inumeric(FLERR,arg[6],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal fix brownian command");

  if (strcmp(arg[7],"angmom") == 0) {
    if (strcmp(arg[8],"no") == 0) {
      ang_diff = 0.0;
    }
    else {
      if (!atoms_have_quaternion()) {
        error->all(FLERR, "All fix atoms need to be extended particles");
      }
      ang_diff = utils::numeric(FLERR,arg[8],false,lmp);
      do_orientational_dynamics=true;
    }

  }

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixBrownian::~FixBrownian()
{

  delete random;

}

/* ---------------------------------------------------------------------- */

int FixPropelSelf::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}


void FixNVE::init()
{
  if(do_orientational_dynamics) avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");

  reset_dt();
}

/* ----------------------------------------------------------------------
   set current t_target and tsqrt
------------------------------------------------------------------------- */

void FixBrownian::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // Only homogeneous temperature supported
  t_target = t_start + delta * (t_stop-t_start);
  tsqrt = sqrt(t_target);

}


/* ---------------------------------------------------------------------- */

void FixBrownian::initial_integrate(int /*vflag*/)
{
  double mdpref;
  double mspref;

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;

  double *quat;
  AtomVecEllipsoid::Bonus *bonus;
  int *ellipsoid;
  double ome[3];

  if(do_orientational_dynamics){
    bonus = avec->bonus;
    ellipsoid = atom->ellipsoid;
  }

  int *type = atom->type
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set square root of temperature
  compute_target();

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if(rmass){
        mdpref=dpref/rmass[i];
        mspref=spref*tsqrt/sqrt(rmass[i]);
      }
      else{
        mdpref=dpref/mass[type[i]];
        mspref=spref*tsqrt/mass[type[i]];
      }

      if(do_orientational_dynamics) quat = bonus[ellipsoid[i]].quat;

      // update position by 1 step
      x[i][0] += mdpref * f[i][0]+mspref*random->gaussian();
      x[i][1] += mdpref * f[i][1]+mspref*random->gaussian();
      x[i][2] += mdpref * f[i][2]+mspref*random->gaussian();

      //update quaternion by 1 step
      if(do_orientational_dynamics){
        //calculate random angular displacement
        ome[0]=apref*random->gaussian(); //TODO could add +constant*dt*torque;
        ome[1]=apref*random->gaussian();
        ome[2]=apref*random->gaussian();
        //rotate quaternion
        std::array<double,4> quat0 = quat;
        MathExtra::vecquat(ome,quat0,quat);

        qnormalize(quat);
      }

    }

}

/* ---------------------------------------------------------------------- */

void FixBrownian::reset_dt()
{
  dtb = update->dt;

  dpref=dtb/gamma;
  spref=sqrt(2*force->boltz*dtb/gamma);
  apref=sqrt(2*force-boltz*dtb*ang_diff);
}

/* ---------------------------------------------------------------------- */

int FixBrownian::atoms_have_quaternion()
{
  if (!atom->ellipsoid_flag) {
    error->all(FLERR, "Mode 'quat' requires atom style ellipsoid");
    return 0;
  }

  int *mask = atom->mask;
  int flag=0,flagall=0;

  // Make sure all atoms have ellipsoid data:

  for (int i = 0; i < atom->nlocal; ++i)
    if (mask[i] & groupbit)
      if (atom->ellipsoid[i] < 0) ++flag;

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall > 0) return 0;

  return 1;
}
