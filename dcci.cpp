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
   Contributing author: Samuel Cajahuaringa Macollunco (CEFET-MG/BR)
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
#include "citeme.h"

using namespace LAMMPS_NS;

static const char cite_dcci[] =
  "dcci command:\n\n"
  "@article{cajahuaringa2021,\n"
  "  author={Samuel Cajahuaringa and Alex Antonelli},\n"
  "  title={Non-equilibrium free-energy calculation of phase-boundaries using LAMMPS},\n"
  "  journal={Comput. Mater. Sci.},\n"
  "  volume={207},\n"
  "  pages={111275},\n"
  "  year={2022},\n"
  "  publisher={-}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

DCCI::DCCI(LAMMPS *lmp) : Pointers(lmp) {
  if (lmp->citeme) lmp->citeme->add(cite_dcci); 
}

/* ---------------------------------------------------------------------- */

DCCI::~DCCI()
{
  MPI_Comm_free(&roots);
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform a dynamic Clausius-Clapeyron integration (dcci) method
------------------------------------------------------------------------- */

void DCCI::command(int narg, char **arg)
{  
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to dcci");
  if (domain->box_exist == 0)
    error->all(FLERR,"dcci command before simulation box is defined");
  if (narg != 9) // && narg != 7)
    error->universe_all(FLERR,"Illegal dcci command");
 
  // coexistence condition 
  Tcoex = force->numeric(FLERR,arg[0]);  // initial coexistence temperature 
  Pcoex = force->numeric(FLERR,arg[1]);  // initial coexistence pressure

  lambda = force->numeric(FLERR,arg[2]); // lambda parameter need define too in the fix adapt/dcci
  lambda_initial = lambda;               // initial lambda parameter

  // Get and check if adapt/dcci fix exists
  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"fix adapt/dcci ID is not defined");
  fix_adapt_dcci = (FixAdaptDCCI*)(modify->fix[whichfix]);

  // Check input values lambdas should be equal, assign other dcci values
  if (lambda != fix_adapt_dcci->lambda)
    error->universe_all(FLERR,"Lambda from fix adapt/dcci in the same world"
        " must be the same");

  // ID of the fix that will control temperature and pressure during the run  
  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[4],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"fix ID to scaling is not defined");
 
  // Time scaling variable
  t_sc = force->bnumeric(FLERR,arg[5]);  // total scaling steps
  if (t_sc < 0)
    error->all(FLERR,"Invalid dcci command");

  if (strcmp(arg[6],"temp") == 0) {
    dcci_flag = 0;
    T_start = force->numeric(FLERR,arg[7]);
    T_end = force->numeric(FLERR,arg[8]);
    if (T_start == T_end) error->all(FLERR,"Illegal dcci command");
    lambda_final = T_start/T_end; // end lambda value
    sf = 2;
  } else if (strcmp(arg[6],"press") == 0) {
    dcci_flag = 1;
    P_start = force->numeric(FLERR,arg[7]); // initial external pressure
    P_end = force->numeric(FLERR,arg[8]);   // end external pressure
    if (P_start == P_end) error->all(FLERR,"Illegal dcci command");
    sf = 1;
  } else error->all(FLERR,"Illegal dcci command");
 
  atomic_flag  = atom->molecular;

  if (atomic_flag != 0) 
    error->all(FLERR,"The dcci command is apply only to atomic system"); 

  // dCCI must be appropriate for temperature and pressure control,
  // i.e. it needs to provide a working Fix::reset_pressure() and must also
  // change the volume. This currently only applies to fix npt and fix nph 

  if ((strncmp(modify->fix[whichfix]->style,"npt",3) != 0)
      && (strncmp(modify->fix[whichfix]->style,"nph",3) != 0))
    error->universe_all(FLERR,"controling of pressure fix is not supported");

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

  // pe_compute = ptr to thermo_pe compute
  // notify compute it will be called at first swap

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"dcci could not find thermo_pe compute");
  Compute *pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep+1);  // compute the potential energy of the system for each step


  int color;
  if (me == 0) color = 0;  
  else color = 1;         
  MPI_Comm_split(universe->uworld,color,0,&roots); 

  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots); 
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world); 
  
  // setup dcci runs

  double lambda_k,lambda_k_1,press_k,press_k_1,press_rs_k,press_rs_k_1;
  lambda_k_1 = lambda_k = lambda; 
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
    MPI_Allgather(&natoms,1,MPI_INT,NATOMS,1,MPI_INT,roots); 
  MPI_Bcast(NATOMS,nworlds,MPI_INT,0,world); 

  PE = new double[nworlds];
  if (me == 0)
    MPI_Allgather(&pe,1,MPI_DOUBLE,PE,1,MPI_DOUBLE,roots); 
  MPI_Bcast(PE,nworlds,MPI_DOUBLE,0,world); 

  VOL = new double[nworlds];
  if (me == 0)
    MPI_Allgather(&vol,1,MPI_DOUBLE,VOL,1,MPI_DOUBLE,roots); 
  MPI_Bcast(VOL,nworlds,MPI_DOUBLE,0,world); 
 
  double du,dv;

  du = (PE[0]/NATOMS[0] - PE[1]/NATOMS[1]);
  dv = (VOL[0]/NATOMS[0] - VOL[1]/NATOMS[1]);

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting dcci ...\n");

  update->integrate->setup(1);

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step  Tcoex  Pcoex  lambda Pcoex_rs pe1  pe2  vol1  vol2");
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step  Tcoex  Pcoex  lambda Pcoex_rs pe1  pe2  vol1  vol2");
      fprintf(universe->ulogfile,"\n");
    }
    print_status(); 
  }

  double ts = update->ntimestep + 1 - update->beginstep;
  ts /= (update->endstep - update->beginstep);

  nktv2p = force->nktv2p;

  if (dcci_flag == 0) {
    lambda_k_1 = scaling_function(lambda_initial,lambda_final,ts); 
    press_rs_k_1 = press_rs_k - (lambda_k_1 - lambda_k)*du/dv*nktv2p;   
    press_k_1 = press_rs_k_1 / lambda_k_1;
  } else if (dcci_flag == 1) {
    press_k_1 = scaling_function(P_start,P_end,ts); 
    lambda_k_1 = lambda_k * (1.0 + (press_k/(du/dv*nktv2p))) / (1.0 + (press_k_1/(du/dv*nktv2p)));
    press_rs_k_1 = lambda_k_1 * press_k_1; 
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
      MPI_Allgather(&pe,1,MPI_DOUBLE,PE,1,MPI_DOUBLE,roots); 
    MPI_Bcast(PE,nworlds,MPI_DOUBLE,0,world); 

    if (me == 0)
      MPI_Allgather(&vol,1,MPI_DOUBLE,VOL,1,MPI_DOUBLE,roots); 
    MPI_Bcast(VOL,nworlds,MPI_DOUBLE,0,world); 

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
   scaling functions
------------------------------------------------------------------------- */

double DCCI::scaling_function(double vi, double vf, double t)
{
  if (sf == 2) return vi / (1 + t * (vi/vf - 1)); // <- for control of temperature
  return vi + (vf - vi) * t; // <- for control of pressure
}

