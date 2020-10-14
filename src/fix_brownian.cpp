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

#include "domain.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "math_extra.h"
#include "random_mars.h"
#include "update.h"
#include "force.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;


/* ---------------------------------------------------------------------- */

FixBrownian::FixBrownian(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if ((narg != 6) && (narg != 8) && (narg != 9)) error->all(FLERR,"Illegal fix brownian command");

  temp = utils::numeric(FLERR,arg[3],false,lmp);
  if (temp < 0.0) error->all(FLERR,"Fix brownian temperature must be >= 0.0");
  gamma = utils::numeric(FLERR,arg[4],false,lmp);
  if (gamma <= 0.0) error->all(FLERR,"Fix brownian gamma must be > 0.0");
  seed = utils::inumeric(FLERR,arg[5],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal fix brownian command");

  if(narg>6){
    if(strcmp(arg[6],"angdiff") == 0){
      if (strcmp(arg[7],"no") == 0) {
        ang_gamma = 0.0;
        do_orientational_dynamics=false;
      }
      else {
        if (atoms_have_quaternion()){
          ang_gamma = utils::numeric(FLERR,arg[7],false,lmp);
          do_orientational_dynamics=true;
        }
        else{
          error->all(FLERR, "All fix atoms need to be extended particles");
        }
      }
    }
    else{
      error->all(FLERR,"Illegal fix brownian command");
    }
  }

  if(narg>8){
    if(strcmp(arg[8],"exact") == 0){
      do_exact_rotation=true;
    }
    else{
      error->all(FLERR,"Illegal fix brownian command");
    }
  }

  if(domain->dimension == 2) is_2d=true;

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixBrownian::~FixBrownian()
{

  delete random;

}

/* ---------------------------------------------------------------------- */

int FixBrownian::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}


void FixBrownian::init()
{
  if(do_orientational_dynamics) avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");

  reset_dt();
}

/* ---------------------------------------------------------------------- */

void FixBrownian::initial_integrate(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;

  double *quat;
  AtomVecEllipsoid::Bonus *bonus;
  int *ellipsoid;
  double rot[3];
  double qrot[4]; //for ALG1
  double dquat[4]; //for ALG2

  if(do_orientational_dynamics){
    bonus = avec->bonus;
    ellipsoid = atom->ellipsoid;
  }

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      if(do_orientational_dynamics) quat = bonus[ellipsoid[i]].quat;

      // update position by 1 step
      x[i][0] += dtpref * f[i][0]+stpref*random->gaussian();
      x[i][1] += dtpref * f[i][1]+stpref*random->gaussian();
      if(!is_2d) x[i][2] += dtpref * f[i][2]+stpref*random->gaussian();

      //update quaternion by 1 step
      if(do_orientational_dynamics){
        //calculate random angular displacement
        //TODO test nonzero torque
        if(is_2d){
          rot[0]=0.0;
          rot[1]=0.0;
        }
        else{
          rot[0]=drpref*torque[i][0]+srpref*random->gaussian();
          rot[1]=drpref*torque[i][1]+srpref*random->gaussian();
        }
        rot[2]=drpref*torque[i][2]+srpref*random->gaussian();

        //rotate quaternion
        if(do_exact_rotation){
          //ALG1 (exact)
          rvec_to_quat(rot,qrot);
          double quat0[4];
          std::copy(quat,quat+4,quat0);
          MathExtra::quatquat(qrot,quat0,quat);
        }
        else{
          //ALG2 (approximate)
          MathExtra::vecquat(rot,quat,dquat);
          quat[0]=quat[0]+0.5*dquat[0];
          quat[1]=quat[1]+0.5*dquat[1];
          quat[2]=quat[2]+0.5*dquat[2];
          quat[3]=quat[3]+0.5*dquat[3];
        }

        //normalize quaternion to avoid error propagation
        qnormalize(quat);
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixBrownian::reset_dt()
{
  dtb = update->dt;

  dtpref=dtb/gamma;
  drpref=dtb/ang_gamma;
  stpref=sqrt(2*force->boltz*temp*dtb/gamma);
  srpref=sqrt(2*force->boltz*temp*dtb/ang_gamma);
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
