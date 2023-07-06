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
}

parameters {
  simplex[N_cls]   class_prob;
  vector[N_cls]    eta[N_itm];
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

    for (c in 1:N_cls) {
      class_membership[p][c] += ordered_logistic_lpmf(y[o] | eta[i][c], cutpoints[i]);
    }
  }

  for (p in 1:N_ppt) {
    logLik += weight[p] * log_sum_exp(class_membership[p]);
    class_membership[p] = softmax(class_membership[p]);
  }
  bic = ((N_cls - 1) + (N_cls * N_itm) + N_itm * (N_cls - 2)) * log(sum(weight)) - 2 * logLik;

  // for (ppt in 1:N_ppt) {
  //   vector[N_cls] ll;
  //
  //   ll = log(class_prob);
  //   for (cls in 1:N_cls) {
  //     for (itm in 1:N_itm) {
  //       ll[cls] += ordered_logistic_lpmf(y[ppt, itm] | eta[itm][cls], append_row(0.0, cutpoints[itm]));
  //     }
  //   }
  //
  //   logLik += log_sum_exp(ll);
  //   class_membership[ppt] = softmax(ll);
  // }

  // for (itm in 1:N_itm) {
  //   for (cls in 1:N_cls) {
  //     // probs[itm, cls][1] = 1 - inv_logit(eta[itm, cls] - cutpoints[itm, 1      ]);
  //     // probs[itm, cls][N_cat] = inv_logit(eta[itm, cls] - cutpoints[itm, N_cat-1]);
  //     // for(cat in 2:(N_cat-1)) {
  //     //   probs[itm, cls][cat] = inv_logit(eta[itm, cls] - cutpoints[itm, cat-1  ]) - inv_logit(eta[itm, cls] - cutpoints[itm, cat]);
  //     // }
  //     probs[itm, cls][1] = 1 - inv_logit(eta[itm, cls]);
  //     probs[itm, cls][N_cat] = inv_logit(eta[itm, cls] - cutpoints[itm, N_cat-2]);
  //     for(cat in 2:(N_cat-2)) {
  //       probs[itm, cls][cat] = inv_logit(eta[itm, cls] - cutpoints[itm, cat-1  ]) - inv_logit(eta[itm, cls] - cutpoints[itm, cat]);
  //     }
  //   }
  // }
  //
}

model {
  target += logLik;
  // for (cls in 1:N_cls) {
  //   target += std_normal_lpdf(cutpoints[cls]);
  // }
  // for (itm in 1:N_itm) {
  //   target += std_normal_lpdf(eta[itm]);
  // }
  // target += dirichlet_lpdf(class_prob | rep_vector(1.0, N_cls));
}


