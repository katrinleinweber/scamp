/*
 				catout.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for catout.c
*
*	Last modify:	06/10/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _CATOUT_H_
#define _CATOUT_H_
/*--------------------------------- constants -------------------------------*/
/*--------------------------------- typedefs --------------------------------*/

typedef enum {CAT_NONE, CAT_ASCII_HEAD, CAT_ASCII, CAT_ASCII_SKYCAT,
		CAT_ASCII_VOTABLE, CAT_FITS_LDAC} cattypenum;

/*--------------------------- structure definitions -------------------------*/
typedef struct mergedsample
  {
  double	wcspos[NAXIS];		/* Mean World Coordinate positions */
  float		wcsposerr[NAXIS];	/* Errors on mean WCS positions */
  float		wcsposdisp[NAXIS];	/* Dispersion on mean WCS positions */
  float		wcsprop[NAXIS];		/* Proper motion vectors in WCS */
  float		wcsproperr[NAXIS];	/* Proper motion vector errors in WCS */
  float		epoch;			/* Mean epoch for observations */
  float		flux[MAXPHOTINSTRU];	/* Mean flux */
  float		fluxerr[MAXPHOTINSTRU];	/* Mean flux uncertainty (1-sigma) */
  float		mag[MAXPHOTINSTRU];	/* "Mean" magnitude */
  int		nmag[MAXPHOTINSTRU];	/* Number of available measurements */
  float		magerr[MAXPHOTINSTRU];	/* Mean mag. uncertainty (1-sigma) */
  float		magdisp[MAXPHOTINSTRU];	/* Mean mag. dispersion (1-sigma) */
  int		nband;			/* Number of available bands */
  int		flags;			/* Merging flags */
  }	mergedsamplestruct;

/*-------------------------------- protos -----------------------------------*/

void		writemergedcat_fgroup(char *filename, fgroupstruct *fgroup);

#endif