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
   Contributing author: Oscar Samuel Cajahuaringa Macollunco 
                       (CEFET-MG Curvelo)
   Contact Email: samuelcajahuaringa@gmail.com 
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "dcci.h"
#include "fix_adapt_dcci.h"
#include "universe.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "output.h"
#include "thermo.h"
#include "fix.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DCCI::DCCI(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

DCCI::~DCCI()
{
  MPI_Comm_free(&roots);
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void DCCI::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to dCCI");
  if (domain->box_exist == 0)
    error->all(FLERR,"dCCI command before simulation box is defined");
  //if (narg != 6 && narg != 7)
  //  error->universe_all(FLERR,"Illegal dCCI command");
 
  // coexistence condition 
  Tcoex = force->numeric(FLERR,arg[0]);  // coexistence temperature 
  Pcoex = force->numeric(FLERR,arg[1]);  // coexistence pressure

  lambda = force->numeric(FLERR,arg[2]);
  lambda_initial = lambda;

  // Get and check if adapt/dcci fix exists
  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"fix adapt/dcci ID is not defined");
  fix_adapt_dcci = (FixAdaptDCCI*)(modify->fix[whichfix]);

  // Check input values lambdas should be equal, assign other gREM values
  if (lambda != fix_adapt_dcci->lambda)
    error->universe_all(FLERR,"Lambda from fix adapt/dcci in the same world"
        " must be the same");

  // ID of the fix that will control temperature and pressure during the run  
  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[4],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"fix ID to scaling is not defined");
 
  // Time scaling variable
  t_sc = force->bnumeric(FLERR,arg[5]);  // <- total scaling steps
  if (t_sc < 0)
    error->all(FLERR,"Invalid dCCI command");

  if (strcmp(arg[6],"temp") == 0) {
    dcci_flag = 0;
    T_start = force->numeric(FLERR,arg[7]);
    T_end = force->numeric(FLERR,arg[8]);
    if (T_start == T_end) error->all(FLERR,"Illegal dcci command");
    lambda_final = T_start/T_end;
  } else if (strcmp(arg[6],"press") == 0) {
    dcci_flag = 1;
    P_start = force->numeric(FLERR,arg[7]);
    P_end = force->numeric(FLERR,arg[8]);
    if (P_start == P_end) error->all(FLERR,"Illegal dcci command");
    sf = 1;
  } else error->all(FLERR,"Illegal dcci command");
  
  // scaling parameter
  //lambda_initial = force->numeric(FLERR,arg[4]); // <- intial lambda value
  //lambda_final = force->numeric(FLERR,arg[5]);  // <- final lambda value

  //if ((lambda_initial <= 0.0) || (lambda_final <= 0.0))
  //  error->all(FLERR,"Invalid dCCI command");

  // choose scaling function
  if (narg > 9 && dcci_flag == 0) {
    if (strcmp(arg[9], "function") == 0) sf = force->inumeric(FLERR,arg[10]);
    else error->all(FLERR,"Illegal dCCI scaling function");
    if ((sf!=1) && (sf!=2))
      error->all(FLERR,"Illegal dCCI scaling function");
  } else if (narg > 9 && dcci_flag == 1) error->all(FLERR,"Illegal dCCI command"); 

  // dCCI must be appropriate for temperature and pressure control,
  // i.e. it needs to provide a working Fix::reset_pressure() and must also
  // change the volume. This currently only applies to fix npt and
  // fix rigid/npt variants

  if ((strncmp(modify->fix[whichfix]->style,"npt",3) != 0)
      && (strncmp(modify->fix[whichfix]->style,"rigid/npt",9) != 0)
      && (strncmp(modify->fix[whichfix]->style,"nph",3) != 0))
    error->universe_all(FLERR,"controling of temperature and pressure fix is not supported");

  // setup for long tempering run

  update->whichflag = 1;
  update->nsteps = t_sc;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + t_sc;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // local storage

  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  
  nworlds = universe->nworlds;
  if (nworlds!=2) error->all(FLERR,"dcci apply only for two systems");  
  iworld = universe->iworld;
  nktv2p = force->nktv2p;

  //printf("me %i iworld %i me_universe %i nworlds %i\n",me,iworld,me_universe,nworlds);

  // pe_compute = ptr to thermo_pe compute
  // notify compute it will be called at first swap

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"dcci could not find thermo_pe compute");
  Compute *pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep+1);  // <- compute the potential energy for each step

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;   // <- father proc for each world i
  else color = 1;           // <- child proc for each world i
  MPI_Comm_split(universe->uworld,color,0,&roots); // <- roots is the communicator word between father proc of each world

  // world2root[i] = global proc that is root proc of world i

  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots); // <- send the id's of father proc for each father proc 
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world); // <- root of entired worl send the id's father proc for both worlds

  // setup dcci runs

  double lambda_k,lambda_k_1,press_k,press_k_1,press_rs_k,press_rs_k_1;
  lambda_k_1 = lambda_k = lambda; // <- step 0
  press_k_1 = press_k = Pcoex;
  press_rs_k_1 = press_rs_k = Pcoex_rs = Pcoex;
  
  int nlocal = atom->nlocal;
  int natoms;
  double pe;
  pe = pe_compute->compute_scalar();
  pe_compute->addstep(update->ntimestep + 1);

  double boxlox=domain->boxlo[0];
  double boxhix=domain->boxhi[0];
  double boxloy=domain->boxlo[1];
  double boxhiy=domain->boxhi[1];
  double boxloz=domain->boxlo[2];
  double boxhiz=domain->boxhi[2];
  double vol = (boxhix - boxlox)*(boxhiy - boxloy)*(boxhiz - boxloz);
  
  MPI_Allreduce(&nlocal,&natoms,1,MPI_INT,MPI_SUM,world);

  NATOMS = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&natoms,1,MPI_INT,NATOMS,1,MPI_INT,roots); // <- send the id's of father proc for each father proc 
  MPI_Bcast(NATOMS,nworlds,MPI_INT,0,world); // <- root of entired worl send the id's father proc for both worlds

  PE = new double[nworlds];
  if (me == 0)
    MPI_Allgather(&pe,1,MPI_DOUBLE,PE,1,MPI_DOUBLE,roots); // <- send the id's of father proc for each father proc 
  MPI_Bcast(PE,nworlds,MPI_DOUBLE,0,world); // <- root of entired worl send the id's father proc for both worlds

  VOL = new double[nworlds];
  if (me == 0)
    MPI_Allgather(&vol,1,MPI_DOUBLE,VOL,1,MPI_DOUBLE,roots); // <- send the id's of father proc for each father proc 
  MPI_Bcast(VOL,nworlds,MPI_DOUBLE,0,world); // <- root of entired worl send the id's father proc for both worlds
 
  double du,dv;

  du = (PE[0]/NATOMS[0] - PE[1]/NATOMS[1]);
  dv = (VOL[0]/NATOMS[0] - VOL[1]/NATOMS[1]);

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting dCCI ...\n");

  update->integrate->setup(1);

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step  Tcoex  Pcoex  lambda Prs pe1  pe2  vol1  vol2");
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step  Tcoex  Pcoex  lambda Prs pe1  pe2  vol1  vol2");
      fprintf(universe->ulogfile,"\n");
    }
    print_status(); //<- step Tcoex Pcoex lambda U1 U2 V1 V2
  }

  double ts = update->ntimestep + 1 - update->beginstep;
  ts /= (update->endstep - update->beginstep);

  nktv2p = force->nktv2p;

  if (dcci_flag == 0) {
    lambda_k_1 = scaling_function(lambda_initial,lambda_final,ts); 
    press_rs_k_1 = press_rs_k - (lambda_k_1 - lambda_k)*du/dv*nktv2p;   
    press_k_1 = press_rs_k_1 / lambda_k_1;
  } else if (dcci_flag == 1) {
    press_k_1 = scaling_function(P_start,P_end,ts); // <- coexistence pressure
    lambda_k_1 = lambda_k * (1.0 + (press_k/(du/dv*nktv2p))) / (1.0 + (press_k_1/(du/dv*nktv2p)));
    press_rs_k_1 = lambda_k_1 * press_k_1;
    //lambda_k_1 -= (press_k_1 - press_k)/(du/dv*nktv2p); //<- old
  }

  timer->init();
  timer->barrier_start();
 
  for (int i = 0; i < t_sc; i++) {
    
    fix_adapt_dcci->lambda = lambda_k_1;
    modify->fix[whichfix]->reset_target_pressure(press_rs_k_1);

    update->integrate->run(1);
    
    pe = pe_compute->compute_scalar();
    pe_compute->addstep(update->ntimestep + 1);

    double boxlox=domain->boxlo[0];
    double boxhix=domain->boxhi[0];
    double boxloy=domain->boxlo[1];
    double boxhiy=domain->boxhi[1];
    double boxloz=domain->boxlo[2];
    double boxhiz=domain->boxhi[2];
    double vol = (boxhix - boxlox)*(boxhiy - boxloy)*(boxhiz - boxloz);

    if (me == 0)
      MPI_Allgather(&pe,1,MPI_DOUBLE,PE,1,MPI_DOUBLE,roots); // <- send the id's of father proc for each father proc 
    MPI_Bcast(PE,nworlds,MPI_DOUBLE,0,world); // <- root of entired worl send the id's father proc for both worlds

    if (me == 0)
      MPI_Allgather(&vol,1,MPI_DOUBLE,VOL,1,MPI_DOUBLE,roots); // <- send the id's of father proc for each father proc 
    MPI_Bcast(VOL,nworlds,MPI_DOUBLE,0,world); // <- root of entired worl send the id's father proc for both worlds

    lambda = lambda_k_1;
    Pcoex = press_k_1;
    Pcoex_rs = press_rs_k_1;

    if (me_universe == 0) print_status();
   
    du = (PE[0]/NATOMS[0] - PE[1]/NATOMS[1]);
    dv = (VOL[0]/NATOMS[0] - VOL[1]/NATOMS[1]);

    double ts = update->ntimestep + 1 - update->beginstep;
    if (ts != 0.0) ts /= (update->endstep - update->beginstep);

    if (dcci_flag == 0) {
      lambda_k = lambda_k_1;
      lambda_k_1 = scaling_function(lambda_initial,lambda_final,ts);
      press_rs_k = press_rs_k_1;
      press_rs_k_1 = press_rs_k - (lambda_k_1 - lambda_k)*du/dv*nktv2p;
      press_k_1 = press_rs_k_1 / lambda_k_1;
      
    } else if (dcci_flag == 1) {
      press_k = press_k_1;
      press_k_1 = scaling_function(P_start,P_end,ts);
      lambda_k = lambda_k_1;
      lambda_k_1 = lambda_k * (1.0 + (press_k/(du/dv*nktv2p))) / (1.0 + (press_k_1/(du/dv*nktv2p)));
      press_rs_k_1 = lambda_k_1 * press_k_1;
    }
    
  }

  timer->barrier_stop();

  update->integrate->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void DCCI::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,BIGINT_FORMAT,update->ntimestep);
    fprintf(universe->uscreen,"    %.6f    %.6f    %.8f    %.8f    %.5f    %.5f    %.4f    %.4f\n",Tcoex/lambda,Pcoex,lambda,Pcoex_rs,PE[0],PE[1],VOL[0],VOL[1]);
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,BIGINT_FORMAT,update->ntimestep);
    fprintf(universe->ulogfile,"   %.6f    %.6f    %.8f    %.8f    %.5f    %.5f    %.4f    %.4f\n",Tcoex/lambda,Pcoex,lambda,Pcoex_rs,PE[0],PE[1],VOL[0],VOL[1]);
    fflush(universe->ulogfile);
  }
}

/* ----------------------------------------------------------------------
   scaling function
------------------------------------------------------------------------- */

double DCCI::scaling_function(double vi, double vf, double t)
{
  if (sf == 2) return vi / (1 + t * (vi/vf - 1));
  // Default option is sf = 1.
  return vi + (vf - vi) * t;
}

