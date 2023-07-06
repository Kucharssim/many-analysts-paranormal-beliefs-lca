data {
  int<lower=0> N_obs;
  int<lower=0> N_itm;
  int<lower=0> N_ppt;
  int<lower=1> N_cat;
  int<lower=0> N_cls;
  int<lower=1,upper=N_cat> y[N_obs];
  int<lower=1,upper=N_ppt> ppt[N_obs];
  int<lower=1,upper=N_itm> itm[N_obs];
  real<lower=0> weight[N_ppt];

  int<lower=1,upper=4> N_fac;
  int<lower=1,upper=N_fac> fac[N_obs];
  int<lower=1> N_lvl;
  int<lower=1,upper=N_lvl> lvl[N_obs,N_cls];
}

parameters {
  simplex[N_cls]    class_prob;
  ordered[N_lvl]    eta;
  positive_ordered[N_cat-2] cutpoints_raw[N_itm];
}

transformed parameters {
  vector[N_cls] class_membership[N_ppt];
  ordered[N_cat-1] cutpoints[N_itm];
  real logLik = 0;
  real bic;

  for (p in 1:N_ppt) {
    class_membership[p] = log(class_prob);
  }

  for (i in 1:N_itm) {
    cutpoints[i] = append_row(0.0, cutpoints_raw[i]);
    cutpoints[i] = cutpoints[i] - mean(cutpoints[i]);
  }

  for (o in 1:N_obs) {
    int p = ppt[o];
    int i = itm[o];
    int f = fac[o];

    for (c in 1:N_cls) {
      int l = lvl[o][c];
      real eta_ = eta[l];
      class_membership[p][c] += ordered_logistic_lpmf(y[o] | eta_, cutpoints[i]);
    }
  }

  for (p in 1:N_ppt) {
    logLik += weight[p] * log_sum_exp(class_membership[p]);
    class_membership[p] = softmax(class_membership[p]);
  }
}

model {
  target += logLik;
}

