#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include "Rmath.h"
#include <Eigen/Dense>
using namespace cpp11;

// ------------------------------------------------
Eigen::MatrixXd get_trans_mat_cpp(double s, double sigma, double dx) {

  int nx = round(1 / dx) + 1;
  Eigen::MatrixXd trans_mat(nx, nx);
  trans_mat.fill(0.0);
  trans_mat(0,0) = 1;
  trans_mat(nx-1, nx-1) = 1;

  for (int i = 1; i < (nx - 1); ++i) {

    // mean and scaled SD of diffusion
    double x_i = i*dx;
    double mu = x_i + s*x_i*(1 - x_i);
    double tau = sqrt(sigma*sigma*x_i*(1 - x_i));

    // transition prob for intermediate frequencies
    for (int j = 1; j < (nx - 1); ++j) {
      double x_j = j*dx;
      trans_mat(i,j) = Rf_pnorm5(x_j + dx/2, mu, tau, 1, 0) - Rf_pnorm5(x_j - dx/2, mu, tau, 1, 0);
    }

    // special case at boundaries
    trans_mat(i,0) = Rf_pnorm5(dx/2, mu, tau, 1, 0);
    trans_mat(i,nx-1) = Rf_pnorm5(1 - dx/2, mu, tau, 0, 0);
  }

  return trans_mat;
}

// ------------------------------------------------
double get_forward_ll_cpp(Eigen::MatrixXd trans_mat,
                          integers t_samp,
                          integers n_samp,
                          integers n_pos,
                          int t_min,
                          Eigen::VectorXd state_vec) {

  // get basic quantities
  int nx = trans_mat.rows();
  double dx = 1/double(nx - 1);
  int n_data = t_samp.size();
  int t_max = t_samp[n_data - 1];

  // apply transitions
  double ret = 0.0;
  int j = 0;
  for (int t = t_min; t < (t_max + 1); ++t) {

    if (t > t_min) {
      state_vec = trans_mat.transpose() * state_vec;
    }

    if (t == t_samp[j] && false) {
      for (int i = 0; i < nx; ++i) {
        state_vec(i) *= Rf_dbinom(n_pos[j], n_samp[j], i*dx, 0);
      }
      double vs = state_vec.sum();
      if (vs == 0) {
        return -1e300;
      }
      ret += log(vs);
      state_vec /= vs;
      j++;
    }
  }

  return ret;
}

// ------------------------------------------------
// calculate and return log-likelihood
[[cpp11::register]]
double loglike_cpp11(doubles params, list data, list misc) {

  // extract inputs
  double s = params["s"];
  double sigma = params["sigma"];
  double dx = doubles(misc["dx"])[0];
  int nx = round(1/dx) + 1;
  doubles x_init_raw = misc["x_init"];

  Eigen::VectorXd x_init(nx);
  for (int i = 0; i < nx; ++i) {
    x_init(i) = x_init_raw[i];
  }

  // get transition matrix
  Eigen::MatrixXd trans_mat = get_trans_mat_cpp(s, sigma, dx);

  // sum log-likelihood over all data
  double ret = 0.0;
  for (int i = 0; i < data.size(); ++i) {
    list data_pop = data[i];
    integers t_samp = data_pop["t"];
    integers n_samp = data_pop["n_samp"];
    integers n_pos = data_pop["n_pos"];
    int t_min = 1;
    ret += get_forward_ll_cpp(trans_mat, t_samp, n_samp, n_pos, t_min, x_init);
  }

  // return as SEXP
  return(ret);
}

// ------------------------------------------------
// calculate and return log-prior
[[cpp11::register]]
double logprior_cpp11(doubles params, list misc) {

  // extract parameters
  double s = params["s"];
  double sigma = params["sigma"];

  // calculate logprior
  double ret = Rf_dnorm4(s, 0.0, 1.0, true) + Rf_dnorm4(sigma, 0.0, 1.0, true);

  // return as SEXP
  return(ret);
}
