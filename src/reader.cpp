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

#include "reader.h"

#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

// only proc 0 calls methods of this class, except for constructor/destructor

/* ---------------------------------------------------------------------- */

Reader::Reader(LAMMPS *lmp) : Pointers(lmp)
{
  fp = NULL;
}

/* ----------------------------------------------------------------------
   try to open given file
   generic version for ASCII files that may be compressed
------------------------------------------------------------------------- */

void Reader::open_file(const char *file)
{
  if (fp != NULL) close_file();

  compressed = 0;
  const char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file,"r");
  else {
#ifdef LAMMPS_GZIP
    std::string gunzip = fmt::format("gzip -c -d {}",file);
#ifdef _WIN32
    fp = _popen(gunzip.c_str(),"rb");
#else
    fp = popen(gunzip.c_str(),"r");
#endif

#else
    error->one(FLERR,"Cannot open gzipped file: " + utils::getsyserror());
#endif
  }

  if (fp == NULL)
    error->one(FLERR,fmt::format("Cannot open file {}: {}",
                                 file, utils::getsyserror()));
}

/* ----------------------------------------------------------------------
   close current file if open
   generic version for ASCII files that may be compressed
------------------------------------------------------------------------- */

void Reader::close_file()
{
  if (fp == NULL) return;
  if (compressed) pclose(fp);
  else fclose(fp);
  fp = NULL;
}

/* ----------------------------------------------------------------------
   detect unused arguments
------------------------------------------------------------------------- */

void Reader::settings(int narg, char** /*args*/)
{
  if (narg > 0)
    error->all(FLERR,"Illegal read_dump command");
}
