/*
 *				crossid.c
 *
 * Manage source cross-identifications.
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	This file part of:	SCAMP
 *
 *	Copyright:		(C) 2002-2017 Emmanuel Bertin -- IAP/CNRS/UPMC
 *
 *	License:		GNU General Public License
 *
 *	SCAMP is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 * 	(at your option) any later version.
 *	SCAMP is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *	You should have received a copy of the GNU General Public License
 *	along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
 *
 *	Last modified:		 15/11/2017
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <stdlib.h>
#include <stdbool.h>

#include "define.h"
#include "globals.h"
#include "crossid.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "match.h"
#include "misc.h"
#include "prefs.h"
#include "samples.h"

#ifdef HAVE_CONFIG_H
#include    "config.h"
#endif

/**
 * @fn static double get_limit(fgroupstruct *fgroup, double tolerance)
 * @brief Compute the largest possible error.
 * Compute the largest possible error in pixels allowed in previous matching
 * @return largest possible error
 * @date 13/11/2017
 */
static double get_limit(fgroupstruct *fgroup, double tolerance)
{
    double lim_tmp;
    double lim_result = 0.0;
    int i;

    for (i=0; i < fgroup->naxis; i++) {
        lim_tmp = tolerance / fgroup->meanwcsscale[i];
        if (lim_tmp > lim_result)
            lim_result = lim_tmp;
    }

    return lim_result;
}

/*
 * @fn static void sort_samples(fgroupstruct *fgroup)
 * @brief Sort samples to accelerate further processing and reset pointers.
 * @date 13/11/2017
 */
static void sort_unlink_samples(fgroupstruct *fgroup)
{
    fieldstruct **field = fgroup->field;
    setstruct    *set;
    int i, j;

    for (i=0; i < fgroup->nfield; i++) {

        field[i]->prevfield = NULL;
        field[i]->nextfield = NULL;

        for (j=0; j < field[i]->nset; j++) {

            set = field[i]->set[j];

            sort_samples(set);
            unlink_samples(set);

        }
    }
}


/**
 * @fn static bool overlapping_sets(setstruct *set_A, setstruct *set_B)
 * @brief Check if frames are overlapping
 * If frames are not overlapping, no need to corss them!
 * @param set_A a pointer to the first setstruct
 * @param set_B a pointer to the second setstruct
 * @param rlim maximul distance
 * @param lng longitude
 * @param lat latitude
 * @returns boolean true or false.
 * @date 15/11/2017
 */
static bool overlapping_sets(setstruct *set_A, setstruct *set_B,
        double rlim, int lng, int lat)
{

    if (lat == lng)
        return true;

    if ((set_B->projposmax[lng] + rlim) < set_A->projposmin[lng]
     || (set_B->projposmin[lng] - rlim) > set_A->projposmax[lng]
     || (set_B->projposmax[lat] + rlim) < set_A->projposmin[lat]
     || (set_B->projposmin[lat] - rlim) > set_A->projposmax[lat]) {
        return false;
    }

    return true;
}


/**
 * @fn static bool sample_overlap(samplestruct *sample, setstruct *set
 * double rlim, int lng, int lat)
 * @brief Check if sample can match any objets in set.
 * @param set_A a pointer to the first setstruct
 * @param set_B a pointer to the second setstruct
 * @param rlim maximul distance
 * @param lng longitude
 * @param lat latitude
 * @returns boolean true or false.
 * @date 15/11/2017
 */
static bool sample_overlaps_set(
        samplestruct *sample,
        setstruct    *set,
        double        rlim,
        int           lng,
        int           lat)
{
    if (
            sample->projpos[lat] < (set->projposmin[lat] - rlim) ||
            sample->projpos[lat] > (set->projposmax[lat] + rlim) ||
            sample->projpos[lng] < (set->projposmin[lng] - rlim) ||
            sample->projpos[lng] > (set->projposmax[lng] + rlim))
        return false;

    return true;
}


/**
 * @fn void crossid_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
 *                          double tolerance)
 * @brief Perform source cross-identifications
 * Perform source cross-identifications in a group of fields.
 * Uses the global preferences
 * @param fgroup Pointer to the group of fields,
 * @param reffield Pointer to the reference field,
 * @param tolerance Tolerance (in degree. if angular coordinates).
 * @author E. Bertin (IAP)
 * @date 15/11/2017
 */
void crossid_fgroup(
        fgroupstruct *fgroup,
        fieldstruct  *reffield,
        double        tolerance) {

    fieldstruct	 *field_A, *field_B;
    setstruct    *set_A, *set_B;
    samplestruct *sample_A, *sample_B,
                 *B_match,
                 *previous_sample, *next_sample;

    double sample_A_latmin, sample_A_latmax,
           samples_A_B_dist, samples_A_B_mindist,
           lng_diff, lat_diff, dx, rlim, r2n, r2p;

    int naxis, yaxis, lng, lat;
    int i, j, k, l, m, n, o;

    naxis = fgroup->naxis;
    lng = fgroup->lng;
    lat = fgroup->lat;

    /* Compute the largest possible error in pixels allowed in previous matching */
    rlim = get_limit(fgroup, tolerance);

    /* Sort samples to accelerate further processing and reset pointers */
    sort_unlink_samples(fgroup);

    /* Now start the cross-id loop */
    /* for each fields */
    for (i=1; i<fgroup->nfield; i++) {
        field_A = fgroup->field[i];

        /* for each other fields */
        for (j=0; j<i; j++) {
            field_B = fgroup->field[j];

            /* for each sets from field_A */
            for (k=0; k < field_A->nset; k++) {
                set_A = field_A->set[k];

                /* for each sets from field_B */
                for (l=0; l < field_B->nset; l++) {
                    set_B = field_B->set[l];

                    if (!overlapping_sets(set_A, set_B, rlim, lng, lat))
                        continue;

                    /* for each samples from set_A */
                    for (m=0; m <set_A->nsample; m++) {
                        sample_A = &set_A->sample[m];

                        if (!sample_overlaps_set(sample_A, set_B, rlim, lng, lat))
                            continue;

                        /* XXX ??? */
                        if (lat != lng)
                            yaxis = lat;
                        else
                            yaxis = (naxis < 2) ? 0 : 1;

                        sample_A_latmin = sample_A->projpos[lat] - rlim;
                        sample_A_latmax = sample_A->projpos[lat] + rlim;
                        samples_A_B_mindist = rlim * rlim;
                        B_match = NULL;

                        /* for each samples from set_B */
                        for (n=0; n < set_B->nsample; n++) {
                            sample_B = &set_B->sample[n];

                            /* We did not reach interesting samples */
                            if (sample_B->projpos[yaxis] < sample_A_latmin)
                                continue;

                            /* We will not have any more matching samples */
                            if (sample_B->projpos[yaxis] > sample_A_latmax)
                                break;

                            /* XXX ??? */
                            if (sample_B->nextsamp && sample_B->nextsamp->set->field!=field_A)
                                continue;

                            /* Effectively compare objects */
                            if (lat!=lng) {
                                lng_diff = sample_A->projpos[lng] - sample_B->projpos[lng];
                                lat_diff = sample_A->projpos[lat] - sample_B->projpos[lat];
                                samples_A_B_dist = lng_diff*lng_diff + lat_diff*lat_diff;
                            } else {
                                samples_A_B_dist = 0.0;
                                for (o=0; o<naxis; o++)
                                {
                                    dx = sample_A->projpos[o] - sample_B->projpos[o];
                                    samples_A_B_dist += dx*dx;
                                }
                            }

                            /* Finally select the closest source within the search disk */
                            if (samples_A_B_dist < samples_A_B_mindist) {
                                samples_A_B_mindist = samples_A_B_dist;
                                B_match = sample_B;
                            }
                        }


                        /* Link samples if there is a match */
                        if (B_match) {

                            // TODO this part can not go parallel. Use an OpenMP lock here.
                            r2p = r2n = BIG;
                            if ((previous_sample=sample_A->prevsamp)) {
                                /* Check if it is a better match than the previous one */

                                if (lat!=lng) {
                                    lng_diff = previous_sample->projpos[lng] - sample_A->projpos[lng];
                                    lat_diff = previous_sample->projpos[lat] - sample_A->projpos[lat];
                                    r2p = lng_diff*lng_diff + lat_diff*lat_diff;
                                } else {
                                    r2p = 0.0;
                                    for (i=0; i<naxis; i++) {
                                        dx = previous_sample->projpos[i] - sample_A->projpos[i];
                                        r2p += dx*dx;
                                    }
                                }
                            }

                            /* Check if it is a better match than the previous one */
                            if ((next_sample=B_match->nextsamp)) {

                                if (lat!=lng) {
                                    lng_diff = B_match->projpos[lng] - next_sample->projpos[lng];
                                    lat_diff = B_match->projpos[lat] - next_sample->projpos[lat];
                                    r2n = lng_diff*lng_diff + lat_diff*lat_diff;
                                } else {
                                    r2n = 0.0;
                                    for (i=0; i<naxis; i++)
                                    {
                                        dx = B_match->projpos[i] - next_sample->projpos[i];
                                        r2n += dx*dx;
                                    }
                                }
                            }

                            /*------------ unlink from previous match if this is a better match */
                            if (samples_A_B_mindist<r2p && samples_A_B_mindist<r2n) {
                                if (previous_sample)
                                    previous_sample->nextsamp = NULL;
                                if (next_sample)
                                    next_sample->prevsamp = NULL;
                                sample_A->prevsamp = B_match;
                                B_match->nextsamp = sample_A;
                            }
                        }
                    }
                }
            }
        }
    }


    /* Now bring also the reference field samples to the common projection */
    /* Sort samples to accelerate further processing and reset pointers */
    if (reffield) {
        set_A = reffield->set[0];
        sort_samples(set_A);
        unlink_samples(set_A);

        for (k=0; k<fgroup->nfield; k++) {
            field_B = fgroup->field[k];

            for (l=0; l < field_B->nset; l++) {
                set_B = fgroup->field[k]->set[l];

                if (!overlapping_sets(set_A, set_B, rlim, lng, lat))
                    continue;

                /* for each samples from set_A */
                for (m=0; m <set_A->nsample; m++) {
                    sample_A = &set_A->sample[m];

                    if (!sample_overlaps_set(sample_A, set_B, rlim, lng, lat))
                        continue;
                    /* XXX ??? */
                    if (lat != lng)
                        yaxis = lat;
                    else
                        yaxis = (naxis < 2) ? 0 : 1;


                    sample_A_latmin = sample_A->projpos[lat] - rlim;
                    sample_A_latmax = sample_A->projpos[lat] + rlim;
                    samples_A_B_mindist = rlim * rlim;
                    B_match = NULL;

                    /* for each samples from set_B */
                    for (n=0; n < set_B->nsample; n++) {
                        sample_B = &set_B->sample[n];

                        /* We did not reach interesting samples */
                        if (sample_B->projpos[yaxis] < sample_A_latmin)
                            continue;

                        if (sample_B->prevsamp != NULL)
                            continue;

                        /* Effectively compare objects */
                        if (lat!=lng) {
                            lng_diff = sample_A->projpos[lng] - sample_B->projpos[lng];
                            lat_diff = sample_A->projpos[lat] - sample_B->projpos[lat];
                            samples_A_B_dist = lng_diff*lng_diff + lat_diff*lat_diff;
                        } else {
                            samples_A_B_dist = 0.0;
                            for (o=0; o<naxis; o++)
                            {
                                dx = sample_A->projpos[o] - sample_B->projpos[o];
                                samples_A_B_dist += dx*dx;
                            }
                        }

                        /* Finally select the closest source within the search disk */
                        if (samples_A_B_dist < samples_A_B_mindist) {
                            samples_A_B_mindist = samples_A_B_dist;
                            B_match = sample_B;
                        }
                    }



                    if (B_match) {
                        r2n = BIG;
                        if ((next_sample=sample_A->nextsamp)) {
                            /*------------ Check if it is a better match than the previous one */
                            if (lat!=lng) {
                                lng_diff = next_sample->projpos[lng] - sample_A->projpos[lng];
                                lat_diff = next_sample->projpos[lat] - sample_A->projpos[lat];
                                r2n = lng_diff*lng_diff + lat_diff*lat_diff;
                            } else {
                                r2n = 0.0;
                                for (m=0; m<naxis; m++) {
                                    dx = next_sample->projpos[m] - sample_A->projpos[m];
                                    r2n += dx*dx;
                                }
                            }
                        }
                        /*---------- unlink from previous match if this is a better match */
                        if (samples_A_B_mindist<r2n) {
                            if (next_sample)
                                next_sample->prevsamp = NULL;
                            sample_A->nextsamp = B_match;
                            B_match->prevsamp = sample_A;
                        }
                    }
                }
            }
        }
    }

    return;
}


/**
 * @fn void recenter_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
 * @brief Perform field re centering.
 * Perform field re centering with respect to a reference catalog in a
 * group of fields. Uses the global preferences.
 * @param ptr to the group of fields,
 * @param ptr to the reference field.
 * @author E. Bertin (IAP)
 * @date 09/06/2011
 */
void recenter_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
{
    fieldstruct	*field;
    setstruct	**sets, *set;
    samplestruct	*samp, *samp2;
    double	*offsetbuf[NAXIS],
    offset[NAXIS], rawpos[NAXIS], wcspos[NAXIS], dwcspos[NAXIS];
    int		d,f,s, naxis, nsamp, o,omax;

    NFPRINTF(OUTPUT, "Re-centering fields...");

    set = reffield->set[0];
    naxis = fgroup->naxis;
    omax = 0;
    for (f=0; f<fgroup->nfield; f++)
    {
        o = 0;
        field = fgroup->field[f];
        samp = set->sample;
        for (nsamp=set->nsample; nsamp--; samp++)
        {
            for (samp2=samp; (samp2=samp2->nextsamp);)
            {
                if (samp2->set && samp2->set->field == field)
                {
                    if (o>=omax)
                    {
                        omax += 8192;
                        if (o)
                            for (d=0; d<naxis; d++)
                            {
                                QREALLOC(offsetbuf[d], double, omax);
                            }
                        else
                            for (d=0; d<naxis; d++)
                            {
                                QMALLOC(offsetbuf[d], double, omax);
                            }
                    }
                    for (d=0; d<naxis; d++)
                        offsetbuf[d][o] = samp2->projpos[d] - samp->projpos[d];
                    o++;
                }
            }
        }
        /*-- Compute the median reprojected shift in each dimension */
        for (d=0; d<naxis; d++)
            offset[d] = fast_median(offsetbuf[d], o);
        /*-- Convert it to a shift in world coordinates */
        for (d=0; d<naxis; d++)
            rawpos[d] = field->set[0]->wcs->crpix[d] - offset[d];
        raw_to_wcs(field->set[0]->wcs, rawpos, wcspos);
        for (d=0; d<naxis; d++)
            dwcspos[d] = wcspos[d] - field->set[0]->wcs->crval[d];
        sets = field->set;
        for (s=0; s<field->nset; s++)
            update_wcsll(sets[s]->wcs, dwcspos[set->lng], dwcspos[set->lat]);
    }

    for (d=0; d<naxis; d++)
        free(offsetbuf[d]);

    return;
}


/**
 * @fn int check_fieldoverlap(fieldstruct *field1, fieldstruct *field2)
 * @brief Check if two fields overlap or not.
 * @param field1 Pointer to the first field,
 * @param field2 Pointer to the second field.
 * @return 1 if they overlap, 0 otherwise.
 * @author E. Bertin (IAP)
 * @date 07/02/2005
 ***/
int check_fieldoverlap(fieldstruct *field1, fieldstruct *field2)

{
    setstruct	**pset,
    *set;
    samplestruct	*samp,*samp2;
    int		n,s;

    pset = field1->set;
    set = *(pset++);
    for (s=field1->nset; s--; set=*(pset++))
    {
        samp = set->sample;
        for (n=set->nsample; n--; samp++)
        {
            if (samp->nextsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->nextsamp))
                    if (samp2->set->field == field2)
                        return 1;
            }
            if (samp->prevsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->prevsamp))
                    if (samp2->set->field == field2)
                        return 1;
            }
        }
    }

    /* No link found between both fields */
    return 0;
}


/**
 * @fn int check_fieldphotomoverlap(fieldstruct *field, int instru)
 * @brief Check if a field overlaps a photometric field or not.
 * @param field Pointer to the field to check,
 * @param instru Photometric instrument index.
 * @return Photometric code (1 for genuine, 2 for dummy) if it overlaps,
 * 0 otherwise.
 * @author E. Bertin (IAP)
 * @date 25/02/2005
 */
int check_fieldphotomoverlap(fieldstruct *field, int instru)

{
    setstruct	**pset,
    *set;
    samplestruct	*samp,*samp2;
    int		n,s;

    pset = field->set;
    set = *(pset++);
    for (s=field->nset; s--; set=*(pset++))
    {
        samp = set->sample;
        for (n=set->nsample; n--; samp++)
        {
            if (samp->nextsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->nextsamp))
                    if (samp2->set->field->photomflag
                            && samp2->set->field->photomlabel==instru)
                        return samp2->set->field->photomflag;
            }
            if (samp->prevsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->prevsamp))
                    if (samp2->set->field->photomflag
                            && samp2->set->field->photomlabel==instru)
                        return samp2->set->field->photomflag;
            }
        }
    }

    /* No photometric field found */
    return 0;
}

