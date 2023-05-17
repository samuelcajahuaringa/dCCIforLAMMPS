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
   Contributing author:  Samuel Cajahuaringa Macollunco (CEFET-MG/BR)
   lambda scaling forces parameters works with the dynamical
   Clausius-Clapeyron Integration (dcci) method
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "fix_adapt_dcci.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{PAIR,ATOM};

/* ---------------------------------------------------------------------- */

FixAdaptDCCI::FixAdaptDCCI(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
nadapt(0), adapt(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal fix adapt/dcci command");
  lambda = force->numeric(FLERR,arg[3]);
  if (lambda < 0) error->all(FLERR,"Illegal fix adapt/dcci command");

  dynamic_group_allow = 1;
  create_attribute = 1;
  scalar_flag = 1;
  // count # of adaptations
  nadapt = 0;

  int iarg = 4;  // atribute pair or kspace
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix adapt/dcci command");
      nadapt++;
      iarg += 5;
    } else break;
  }

  if (nadapt == 0) error->all(FLERR,"Illegal fix adapt/dcci command");
  adapt = new Adapt[nadapt];

  // parse keywords

  nadapt = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix adapt/dcci command");
      adapt[nadapt].which = PAIR;
      int n = strlen(arg[iarg+1]) + 1;
      adapt[nadapt].pstyle = new char[n];
      strcpy(adapt[nadapt].pstyle,arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      adapt[nadapt].pparam = new char[n];
      adapt[nadapt].pair = NULL;
      strcpy(adapt[nadapt].pparam,arg[iarg+2]);
      force->bounds(FLERR,arg[iarg+3],atom->ntypes,
                    adapt[nadapt].ilo,adapt[nadapt].ihi);
      force->bounds(FLERR,arg[iarg+4],atom->ntypes,
                    adapt[nadapt].jlo,adapt[nadapt].jhi);
      nadapt++;
      iarg += 5;
    } else break;
  }

  // allocate pair style arrays

  int n = atom->ntypes;
  for (int m = 0; m < nadapt; m++)
    if (adapt[m].which == PAIR)
      memory->create(adapt[m].array_orig,n+1,n+1,"adapt/dcci:array_orig");

}

/* ---------------------------------------------------------------------- */

FixAdaptDCCI::~FixAdaptDCCI()
{
  for (int m = 0; m < nadapt; m++) {
    if (adapt[m].which == PAIR) {
      delete [] adapt[m].pstyle;
      delete [] adapt[m].pparam;
      memory->destroy(adapt[m].array_orig);
    }
  }
  delete [] adapt;
}

/* ---------------------------------------------------------------------- */

int FixAdaptDCCI::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_RUN;
  mask |= PRE_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::init()
{
  int i,j;

  // allow a dynamic group only if ATOM attribute not used

  if (group->dynamic[igroup])
    for (int i = 0; i < nadapt; i++)
      if (adapt[i].which == ATOM)
        error->all(FLERR,"Cannot use dynamic group with fix adapt/dcci atom");

  // setup and error checks

  anypair = 0;

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];

    /*ad->ivar = input->variable->find(ad->var);
    if (ad->ivar < 0)
      error->all(FLERR,"Variable name for fix adapt does not exist");
    if (!input->variable->equalstyle(ad->ivar))
      error->all(FLERR,"Variable for fix adapt is invalid style");*/

    if (ad->which == PAIR) {
      anypair = 1;
      ad->pair = NULL;

      // if ad->pstyle has trailing sub-style annotation ":N",
      //   strip it for pstyle arg to pair_match() and set nsub = N
      // this should work for appended suffixes as well

      int n = strlen(ad->pstyle) + 1;
      char *pstyle = new char[n];
      strcpy(pstyle,ad->pstyle);

      char *cptr;
      int nsub = 0;
      if ((cptr = strchr(pstyle,':'))) {
        *cptr = '\0';
        nsub = force->inumeric(FLERR,cptr+1);
      }

      if (lmp->suffix_enable) {
        int len = 2 + strlen(pstyle) + strlen(lmp->suffix);
        char *psuffix = new char[len];
        strcpy(psuffix,pstyle);
        strcat(psuffix,"/");
        strcat(psuffix,lmp->suffix);
        ad->pair = force->pair_match(psuffix,1,nsub);
        delete[] psuffix;
      }
      if (ad->pair == NULL) ad->pair = force->pair_match(pstyle,1,nsub);
      if (ad->pair == NULL)
        error->all(FLERR,"Fix adapt/dcci pair style does not exist");

      void *ptr = ad->pair->extract(ad->pparam,ad->pdim);
      if (ptr == NULL)
        error->all(FLERR,"Fix adapt/dcci pair style param not supported");

      // for pair styles only parameters that are 2-d arrays in atom types or
      // scalars are supported

      if (ad->pdim != 2 && ad->pdim != 0)
        error->all(FLERR,"Fix adapt/dcci pair style param is not compatible");

      if (ad->pdim == 2) ad->array = (double **) ptr;
      if (ad->pdim == 0) ad->scalar = (double *) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if (strcmp(force->pair_style,"hybrid") == 0 ||
          strcmp(force->pair_style,"hybrid/overlay") == 0) {
        PairHybrid *pair = (PairHybrid *) force->pair;
        for (i = ad->ilo; i <= ad->ihi; i++)
          for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
            if (!pair->check_ijtype(i,j,pstyle))
              error->all(FLERR,"Fix adapt/dcci type pair range is not valid for "
                         "pair hybrid sub-style");
      }

      delete [] pstyle;
    }
  }

  // make copy of original pair array values

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    if (ad->which == PAIR && ad->pdim == 2) {
      for (i = ad->ilo; i <= ad->ihi; i++)
        for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
          ad->array_orig[i][j] = ad->array[i][j];
    } else if (ad->which == PAIR && ad->pdim == 0){
      ad->scalar_orig = *ad->scalar;
    }
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::setup_pre_force(int vflag)
{
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::pre_force(int vflag)
{
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::pre_force_respa(int vflag, int ilevel, int)
{
  if (ilevel < nlevels_respa-1) return;
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdaptDCCI::post_run()
{
  //if (resetflag) restore_settings();
}

/* ----------------------------------------------------------------------
   change pair,kspace parameters based on variable evaluation
------------------------------------------------------------------------- */

void FixAdaptDCCI::change_settings()
{
  int i,j;

  // variable evaluation may invoke computes so wrap with clear/add

  //modify->clearstep_compute();

  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    double value = lambda; 

    // set global scalar or type pair array values

    if (ad->which == PAIR) {
      if (ad->pdim == 0) {
        *ad->scalar = value;
      } else if (ad->pdim == 2) {
          for (i = ad->ilo; i <= ad->ihi; i++)
            for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
              ad->array[i][j] = value;
      }
    }
  }

  //modify->addstep_compute(update->ntimestep + nevery);

  // re-initialize pair styles if any PAIR settings were changed
  // ditto for bond styles if any BOND setitings were changes
  // this resets other coeffs that may depend on changed values,
  //   and also offset and tail corrections

  if (anypair) {
    for (int m = 0; m < nadapt; m++) {
      Adapt *ad = &adapt[m];
      if (ad->which == PAIR) {
        ad->pair->reinit();
      }
    }
  }
}

/* ----------------------------------------------------------------------
   restore pair,kspace,atom parameters to original values
------------------------------------------------------------------------- */

void FixAdaptDCCI::restore_settings()
{
  for (int m = 0; m < nadapt; m++) {
    Adapt *ad = &adapt[m];
    if (ad->which == PAIR) {
      if (ad->pdim == 0) *ad->scalar = ad->scalar_orig;
      else if (ad->pdim == 2) {
        for (int i = ad->ilo; i <= ad->ihi; i++)
          for (int j = MAX(ad->jlo,i); j <= ad->jhi; j++)
            ad->array[i][j] = ad->array_orig[i][j];
      }
    }
  }

  if (anypair) force->pair->reinit();
}

double FixAdaptDCCI::compute_scalar()
{
  return lambda;
}

