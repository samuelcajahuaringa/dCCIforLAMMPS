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

#ifdef COMMAND_CLASS

CommandStyle(dcci,DCCI)

#else

#ifndef LMP_DCCI_H
#define LMP_DCCI_H

#include "pointers.h"

namespace LAMMPS_NS {

class DCCI : protected Pointers {
 public:
  DCCI(class LAMMPS *);
  ~DCCI();
  void command(int, char **);

 private:
  int me,me_universe;          // my proc ID in world and universe
  int iworld,nworlds;          // world info
  double nktv2p;
  MPI_Comm roots;              // MPI comm with 1 root proc from each world
  bigint t_sc;                 // Total scaling steps
  int whichfix;                // index of temperature fix to use
  int fixstyle;                // what kind of temperature fix is used
  int *world2root;             // world2root[i] = root proc of world 

  double Tcoex,Pcoex;          // temperature and pressure coxistence
  double T_start,T_end;
  double P_start,P_end;
  double lambda,Pcoex_rs;
  int dcci_flag;
  int *NATOMS;                 // numbers of atoms of each phase
  double *PE;                  // potential energy of each phase
  double *VOL;                 // volume of each phase
  void print_status();
  double lambda_initial;
  double lambda_final;
  int sf;                      // scaling function option  
  double scaling_function(double, double, double);
  
  class FixAdaptDCCI * fix_adapt_dcci;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Must have more than one processor partition to temper

Cannot use the temper command with only one processor partition.  Use
the -partition command-line option.

E: Temper command before simulation box is defined

The temper command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Tempering fix ID is not defined

The fix ID specified by the temper command does not exist.

E: Invalid frequency in temper command

Nevery must be > 0.

E: Non integer # of swaps in temper command

Swap frequency in temper command must evenly divide the total # of
timesteps.

E: Tempering temperature fix is not valid

The fix specified by the temper command is not one that controls
temperature (nvt or langevin).

E: Too many timesteps

The cumulative timesteps must fit in a 64-bit integer.

E: Tempering could not find thermo_pe compute

This compute is created by the thermo command.  It must have been
explicitly deleted by a uncompute command.

*/
