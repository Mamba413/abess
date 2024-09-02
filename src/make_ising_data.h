#ifndef SRC_MAKE_ISING_DATA_H
#define SRC_MAKE_ISING_DATA_H

#include <Eigen/Eigen>

Eigen::MatrixXd comp_conf(int num_conf, int p);
Eigen::MatrixXd sample_by_conf(long long n, Eigen::MatrixXd theta, int seed);
void iteration(Eigen::VectorXd &sample, Eigen::MatrixXd &theta,
               Eigen::VectorXd &value, int set_seed, int iter_time);
Eigen::MatrixXd Ising_Gibbs(Eigen::MatrixXd theta, int n_sample, int burn, int skip,
                            Eigen::VectorXd value, bool using_seed = false, int set_seed = 1);
#endif //SRC_MAKE_ISING_DATA_H