#ifndef __PNORM_H
#define __PNORM_H

#pragma once
extern void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
extern double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p);
#endif
