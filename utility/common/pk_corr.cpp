#include "Functions.h"

double secondOrder_Pk(int nrun, double pkgen, int mcPartType){
  /* 4 Feb 2010 "crude" quadratic weights */
  double slope = 0.0;
  if (nrun>=20154 && nrun<=20256) slope = -0.05;
  if (nrun>=20268 && nrun<=20303) slope = -0.08;
  if (nrun>=20304 && nrun<=20324) slope = -0.02;
  if (nrun>=20332 && nrun<=20371) slope = +0.08;
  if (nrun>=20387 && nrun<=20486) slope = +0.07;
  if (nrun>=20487 && nrun<=20521) slope = +0.06;
  if (nrun>=20522 && nrun<=20612) slope = +0.01;
  if (nrun>=20613 && nrun<=20695) slope = +0.10;

  double wt_pk = 1 + slope*pow(pkgen-74.0, 2);

  /* Dec 2010 fine corrections to the weights (K+) */
  if (mcPartType) {

    // period 1
    if (nrun>=20114 && nrun<=20154) {
      wt_pk *= (1.0 - 0.005*(pkgen-74.0));
      if (pkgen>77.5) wt_pk *= (1.0+4.0*(pkgen-77.5));
    }
    if (nrun>=20155 && nrun<=20174) {
      wt_pk *= (1.0 + 0.065*(pkgen-74.0));
      if (pkgen>77.5) wt_pk *= (1.0+1.0*(pkgen-77.5));
    }
    if (nrun>=20175 && nrun<=20203) {
      wt_pk *= (1.0 + 0.10*(pkgen-74.0));
      if (pkgen>77.0) wt_pk *= (1.0+0.5*(pkgen-77.0));
      if (pkgen>70.0 && pkgen<71.5) wt_pk *= (1.0-0.5*(71.5-pkgen));
    }

    // periods 2+3
    if (nrun>=20209 && nrun<=20226) {
      wt_pk *= (1.0 + 0.05*(pkgen-74.0));
      if (pkgen>77.0) wt_pk *= (1.0+0.5*(pkgen-77.0));
      if (pkgen>70.0 && pkgen<71.5) wt_pk *= (1.0-0.5*(71.5-pkgen));
    }
    if (nrun>=20228 && nrun<=20256) {
      wt_pk *= (1.0 + 0.035*(pkgen-74.0));
      if (pkgen>76.8) wt_pk *= (1.0+0.7*(pkgen-76.8));
      if (pkgen>70.0 && pkgen<71.6) wt_pk *= (1.0-0.5*(71.6-pkgen));
    }
    if (nrun>=20268 && nrun<=20291) {
      wt_pk *= (1.0 + 0.04*(pkgen-74.0));
      if (pkgen>77.0) wt_pk *= (1.0+6.0*(pkgen-77.0));
    }
    if (nrun>=20296 && nrun<=20303) {
      wt_pk *= (1.0 + 0.015*(pkgen-74.0));
      if (pkgen>76.9) wt_pk *= (1.0+5.0*(pkgen-76.9));
    }
    if (nrun>=20304 && nrun<=20324) {
      wt_pk *= (1.0 + 0.02*(pkgen-74.0));
      if (pkgen>77.2) wt_pk *= (1.0+1.3*(pkgen-77.2));
      if (pkgen>70.0 && pkgen<71.5) wt_pk *= (1.0-0.25*(71.5-pkgen));
    }

    // period 4
    if (nrun>=20332 && nrun<=20351) {
      if (pkgen>77.3) wt_pk *= (1.0+4.0*(pkgen-77.3));
      if (pkgen>70.0 && pkgen<71.0) wt_pk *= (1.0+1.0*(71.0-pkgen));
    }
    if (nrun>=20352 && nrun<=20371) {
      wt_pk *= (1.0 - 0.01*(pkgen-74.0));
      if (pkgen>77.3) wt_pk *= (1.0+6.0*(pkgen-77.3));
      if (pkgen>70.0 && pkgen<71.0) wt_pk *= (1.0+2.0*(71.0-pkgen));
    }
    if (nrun>=20387 && nrun<=20404) {
      wt_pk *= (1.0 - 0.015*(pkgen-74.0));
      if (pkgen>77.3) wt_pk *= (1.0+4.5*(pkgen-77.3));
      if (pkgen>70.0 && pkgen<71.0) wt_pk *= (1.0+4.0*(71.0-pkgen));
    }

    // period 5
    if (nrun>=20410 && nrun<=20424) {
      wt_pk *= (1.0 - 0.01*(pkgen-74.0));
      if (pkgen>77.5) wt_pk *= (1.0+12.0*(pkgen-77.5));
      if (pkgen<71.2) wt_pk *= (1.0+ 2.5*(71.2-pkgen));
    }
    if (nrun>=20438 && nrun<=20453) {
      wt_pk *= (1.0 - 0.015*(pkgen-74.0));
      if (pkgen>77.5) wt_pk *= (1.0+16.0*(pkgen-77.5));
      if (pkgen<71.1) wt_pk *= (1.0+ 5.0*(71.1-pkgen));
    }
    if (nrun>=20459 && nrun<=20478) {
      wt_pk *= (1.0 - 0.01*(pkgen-74.0));
      if (pkgen>77.4) wt_pk *= (1.0+10.0*(pkgen-77.4));
      if (pkgen<71.0) wt_pk *= (1.0+ 3.0*(71.0-pkgen));
    }
    if (nrun>=20482 && nrun<=20485) {
      wt_pk *= (1.0 + 0.00*(pkgen-74.0));
      if (pkgen>77.3) wt_pk *= (1.0+6.0*(pkgen-77.3));
      if (pkgen<71.0) wt_pk *= (1.0+5.0*(71.0-pkgen));
    }
  }

  /* Jan-Feb 2011 fine corrections to the weights (K-) */
  else {

    // period 1
    if (nrun>=20114 && nrun<=20154) {
      wt_pk *= (1.0 - 0.005*(pkgen-74.0));
      if (pkgen>77.5) wt_pk *= (1.0+4.0*(pkgen-77.5));
    }
    if (nrun>=20155 && nrun<=20174) {
      if (pkgen>77.5) wt_pk *= (1.0+1.0*(pkgen-77.5));
    }
    if (nrun>=20175 && nrun<=20203) {
      wt_pk *= (1.0 + 0.05*(pkgen-74.0));
      if (pkgen>77.0) wt_pk *= (1.0+0.5*(pkgen-77.0));
      if (pkgen>70.5 && pkgen<71.5) wt_pk *= (1.0-0.5*(71.5-pkgen));
    }

    // periods 2+3
    if (nrun>=20209 && nrun<=20226) {
      wt_pk *= (1.0 + 0.03*(pkgen-74.0));
      if (pkgen>77.0) wt_pk *= (1.0+0.5*(pkgen-77.0));
      if (pkgen>70.0 && pkgen<71.5) wt_pk *= (1.0-0.5*(71.5-pkgen));
    }
    if (nrun>=20228 && nrun<=20256) {
      wt_pk *= (1.0 + 0.035*(pkgen-74.0));
      if (pkgen>76.8) wt_pk *= (1.0+0.7*(pkgen-76.8));
      if (pkgen>70.0 && pkgen<71.6) wt_pk *= (1.0-0.5*(71.6-pkgen));
    }
    if (nrun>=20268 && nrun<=20291) {
      wt_pk *= (1.0 + 0.04*(pkgen-74.0));
      if (pkgen>77.0) wt_pk *= (1.0+6.0*(pkgen-77.0));
    }
    if (nrun>=20296 && nrun<=20303) {
      wt_pk *= (1.0 + 0.05*(pkgen-74.0));
      if (pkgen>76.9) wt_pk *= (1.0+5.0*(pkgen-76.9));
    }

    // period 6
    if (nrun>=20487 && nrun<=20531) {
      wt_pk *= (1.0 + 0.02*(pkgen-74.0));
      if (pkgen>77.0) wt_pk *= (1.0+4.0*(pkgen-77.0));
      if (pkgen>70.0 && pkgen<71.5) wt_pk *= (1.0+0.4*(71.5-pkgen));
    }
  }

  if (wt_pk<0) wt_pk = 0.01;

  return wt_pk;
}
