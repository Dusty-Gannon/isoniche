// Generated by rstantools.  Do not edit by hand.

/*
    isoniche is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    isoniche is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with isoniche.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_isoniche_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_isoniche");
    reader.add_event(77, 75, "end", "model_isoniche");
    return reader;
}
#include <stan_meta_header.hpp>
class model_isoniche
  : public stan::model::model_base_crtp<model_isoniche> {
private:
        int N;
        std::vector<int> P;
        std::vector<int> K;
        int J;
        matrix_d X;
        matrix_d Z;
        matrix_d G;
        matrix_d y;
public:
    model_isoniche(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_isoniche(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_isoniche_namespace::model_isoniche";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 6;
            validate_non_negative_index("P", "2", 2);
            context__.validate_dims("data initialization", "P", "int", context__.to_vec(2));
            P = std::vector<int>(2, int(0));
            vals_i__ = context__.vals_i("P");
            pos__ = 0;
            size_t P_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < P_k_0_max__; ++k_0__) {
                P[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("K", "2", 2);
            context__.validate_dims("data initialization", "K", "int", context__.to_vec(2));
            K = std::vector<int>(2, int(0));
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            size_t K_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < K_k_0_max__; ++k_0__) {
                K[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "J", "int", context__.to_vec());
            J = int(0);
            vals_i__ = context__.vals_i("J");
            pos__ = 0;
            J = vals_i__[pos__++];
            current_statement_begin__ = 9;
            validate_non_negative_index("X", "N", N);
            validate_non_negative_index("X", "(get_base1(P, 1, \"P\", 1) + get_base1(P, 2, \"P\", 1))", (get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1)));
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(N,(get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1))));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, (get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1)));
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = (get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1));
            size_t X_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("Z", "N", N);
            validate_non_negative_index("Z", "(get_base1(K, 1, \"K\", 1) + get_base1(K, 2, \"K\", 1))", (get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1)));
            context__.validate_dims("data initialization", "Z", "matrix_d", context__.to_vec(N,(get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1))));
            Z = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, (get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1)));
            vals_r__ = context__.vals_r("Z");
            pos__ = 0;
            size_t Z_j_2_max__ = (get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1));
            size_t Z_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < Z_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < Z_j_1_max__; ++j_1__) {
                    Z(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("G", "N", N);
            validate_non_negative_index("G", "J", J);
            context__.validate_dims("data initialization", "G", "matrix_d", context__.to_vec(N,J));
            G = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, J);
            vals_r__ = context__.vals_r("G");
            pos__ = 0;
            size_t G_j_2_max__ = J;
            size_t G_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < G_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < G_j_1_max__; ++j_1__) {
                    G(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 12;
            validate_non_negative_index("y", "N", N);
            validate_non_negative_index("y", "2", 2);
            context__.validate_dims("data initialization", "y", "matrix_d", context__.to_vec(N,2));
            y = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, 2);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_2_max__ = 2;
            size_t y_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < y_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                    y(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 18;
            validate_non_negative_index("beta_1", "get_base1(P, 1, \"P\", 1)", get_base1(P, 1, "P", 1));
            num_params_r__ += get_base1(P, 1, "P", 1);
            current_statement_begin__ = 19;
            validate_non_negative_index("beta_2", "get_base1(P, 2, \"P\", 1)", get_base1(P, 2, "P", 1));
            num_params_r__ += get_base1(P, 2, "P", 1);
            current_statement_begin__ = 20;
            validate_non_negative_index("zeta_1", "get_base1(K, 1, \"K\", 1)", get_base1(K, 1, "K", 1));
            num_params_r__ += get_base1(K, 1, "K", 1);
            current_statement_begin__ = 21;
            validate_non_negative_index("zeta_2", "get_base1(K, 2, \"K\", 1)", get_base1(K, 2, "K", 1));
            num_params_r__ += get_base1(K, 2, "K", 1);
            current_statement_begin__ = 22;
            validate_non_negative_index("gamma", "J", J);
            num_params_r__ += J;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_isoniche() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 18;
        if (!(context__.contains_r("beta_1")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_1 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_1");
        pos__ = 0U;
        validate_non_negative_index("beta_1", "get_base1(P, 1, \"P\", 1)", get_base1(P, 1, "P", 1));
        context__.validate_dims("parameter initialization", "beta_1", "vector_d", context__.to_vec(get_base1(P, 1, "P", 1)));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_1(get_base1(P, 1, "P", 1));
        size_t beta_1_j_1_max__ = get_base1(P, 1, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_1_j_1_max__; ++j_1__) {
            beta_1(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_1: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 19;
        if (!(context__.contains_r("beta_2")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_2 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_2");
        pos__ = 0U;
        validate_non_negative_index("beta_2", "get_base1(P, 2, \"P\", 1)", get_base1(P, 2, "P", 1));
        context__.validate_dims("parameter initialization", "beta_2", "vector_d", context__.to_vec(get_base1(P, 2, "P", 1)));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_2(get_base1(P, 2, "P", 1));
        size_t beta_2_j_1_max__ = get_base1(P, 2, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_2_j_1_max__; ++j_1__) {
            beta_2(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_2: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 20;
        if (!(context__.contains_r("zeta_1")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable zeta_1 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("zeta_1");
        pos__ = 0U;
        validate_non_negative_index("zeta_1", "get_base1(K, 1, \"K\", 1)", get_base1(K, 1, "K", 1));
        context__.validate_dims("parameter initialization", "zeta_1", "vector_d", context__.to_vec(get_base1(K, 1, "K", 1)));
        Eigen::Matrix<double, Eigen::Dynamic, 1> zeta_1(get_base1(K, 1, "K", 1));
        size_t zeta_1_j_1_max__ = get_base1(K, 1, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_1_j_1_max__; ++j_1__) {
            zeta_1(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(zeta_1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable zeta_1: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 21;
        if (!(context__.contains_r("zeta_2")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable zeta_2 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("zeta_2");
        pos__ = 0U;
        validate_non_negative_index("zeta_2", "get_base1(K, 2, \"K\", 1)", get_base1(K, 2, "K", 1));
        context__.validate_dims("parameter initialization", "zeta_2", "vector_d", context__.to_vec(get_base1(K, 2, "K", 1)));
        Eigen::Matrix<double, Eigen::Dynamic, 1> zeta_2(get_base1(K, 2, "K", 1));
        size_t zeta_2_j_1_max__ = get_base1(K, 2, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_2_j_1_max__; ++j_1__) {
            zeta_2(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(zeta_2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable zeta_2: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 22;
        if (!(context__.contains_r("gamma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable gamma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("gamma");
        pos__ = 0U;
        validate_non_negative_index("gamma", "J", J);
        context__.validate_dims("parameter initialization", "gamma", "vector_d", context__.to_vec(J));
        Eigen::Matrix<double, Eigen::Dynamic, 1> gamma(J);
        size_t gamma_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < gamma_j_1_max__; ++j_1__) {
            gamma(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(gamma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable gamma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 18;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_1;
            (void) beta_1;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_1 = in__.vector_constrain(get_base1(P, 1, "P", 1), lp__);
            else
                beta_1 = in__.vector_constrain(get_base1(P, 1, "P", 1));
            current_statement_begin__ = 19;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_2;
            (void) beta_2;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_2 = in__.vector_constrain(get_base1(P, 2, "P", 1), lp__);
            else
                beta_2 = in__.vector_constrain(get_base1(P, 2, "P", 1));
            current_statement_begin__ = 20;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> zeta_1;
            (void) zeta_1;  // dummy to suppress unused var warning
            if (jacobian__)
                zeta_1 = in__.vector_constrain(get_base1(K, 1, "K", 1), lp__);
            else
                zeta_1 = in__.vector_constrain(get_base1(K, 1, "K", 1));
            current_statement_begin__ = 21;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> zeta_2;
            (void) zeta_2;  // dummy to suppress unused var warning
            if (jacobian__)
                zeta_2 = in__.vector_constrain(get_base1(K, 2, "K", 1), lp__);
            else
                zeta_2 = in__.vector_constrain(get_base1(K, 2, "K", 1));
            current_statement_begin__ = 22;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> gamma;
            (void) gamma;  // dummy to suppress unused var warning
            if (jacobian__)
                gamma = in__.vector_constrain(J, lp__);
            else
                gamma = in__.vector_constrain(J);
            // model body
            {
            current_statement_begin__ = 29;
            validate_non_negative_index("B", "(get_base1(P, 1, \"P\", 1) + get_base1(P, 2, \"P\", 1))", (get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1)));
            validate_non_negative_index("B", "2", 2);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> B((get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1)), 2);
            stan::math::initialize(B, DUMMY_VAR__);
            stan::math::fill(B, DUMMY_VAR__);
            current_statement_begin__ = 30;
            validate_non_negative_index("Zeta", "(get_base1(K, 1, \"K\", 1) + get_base1(K, 2, \"K\", 1))", (get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1)));
            validate_non_negative_index("Zeta", "2", 2);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Zeta((get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1)), 2);
            stan::math::initialize(Zeta, DUMMY_VAR__);
            stan::math::fill(Zeta, DUMMY_VAR__);
            current_statement_begin__ = 31;
            validate_non_negative_index("mu", "N", N);
            validate_non_negative_index("mu", "2", 2);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> mu(N, 2);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 32;
            validate_non_negative_index("sigma", "N", N);
            validate_non_negative_index("sigma", "2", 2);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> sigma(N, 2);
            stan::math::initialize(sigma, DUMMY_VAR__);
            stan::math::fill(sigma, DUMMY_VAR__);
            current_statement_begin__ = 33;
            validate_non_negative_index("rho", "N", N);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> rho(N);
            stan::math::initialize(rho, DUMMY_VAR__);
            stan::math::fill(rho, DUMMY_VAR__);
            current_statement_begin__ = 34;
            validate_non_negative_index("Omega", "2", 2);
            validate_non_negative_index("Omega", "2", 2);
            validate_non_negative_index("Omega", "N", N);
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic>  > Omega(N, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic>(2, 2));
            stan::math::initialize(Omega, DUMMY_VAR__);
            stan::math::fill(Omega, DUMMY_VAR__);
            current_statement_begin__ = 35;
            validate_non_negative_index("Sigma", "2", 2);
            validate_non_negative_index("Sigma", "2", 2);
            validate_non_negative_index("Sigma", "N", N);
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic>  > Sigma(N, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic>(2, 2));
            stan::math::initialize(Sigma, DUMMY_VAR__);
            stan::math::fill(Sigma, DUMMY_VAR__);
            current_statement_begin__ = 38;
            stan::model::assign(B, 
                        stan::model::cons_list(stan::model::index_min_max(1, get_base1(P, 1, "P", 1)), stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list())), 
                        beta_1, 
                        "assigning variable B");
            current_statement_begin__ = 39;
            stan::model::assign(B, 
                        stan::model::cons_list(stan::model::index_min_max((get_base1(P, 1, "P", 1) + 1), (get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1))), stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list())), 
                        rep_vector(0, get_base1(P, 2, "P", 1)), 
                        "assigning variable B");
            current_statement_begin__ = 40;
            stan::model::assign(B, 
                        stan::model::cons_list(stan::model::index_min_max(1, get_base1(P, 1, "P", 1)), stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list())), 
                        rep_vector(0, get_base1(P, 1, "P", 1)), 
                        "assigning variable B");
            current_statement_begin__ = 41;
            stan::model::assign(B, 
                        stan::model::cons_list(stan::model::index_min_max((get_base1(P, 1, "P", 1) + 1), (get_base1(P, 1, "P", 1) + get_base1(P, 2, "P", 1))), stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list())), 
                        beta_2, 
                        "assigning variable B");
            current_statement_begin__ = 44;
            stan::model::assign(Zeta, 
                        stan::model::cons_list(stan::model::index_min_max(1, get_base1(K, 1, "K", 1)), stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list())), 
                        zeta_1, 
                        "assigning variable Zeta");
            current_statement_begin__ = 45;
            stan::model::assign(Zeta, 
                        stan::model::cons_list(stan::model::index_min_max((get_base1(K, 1, "K", 1) + 1), (get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1))), stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list())), 
                        rep_vector(0, get_base1(K, 2, "K", 1)), 
                        "assigning variable Zeta");
            current_statement_begin__ = 46;
            stan::model::assign(Zeta, 
                        stan::model::cons_list(stan::model::index_min_max(1, get_base1(K, 1, "K", 1)), stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list())), 
                        rep_vector(0, get_base1(K, 1, "K", 1)), 
                        "assigning variable Zeta");
            current_statement_begin__ = 47;
            stan::model::assign(Zeta, 
                        stan::model::cons_list(stan::model::index_min_max((get_base1(K, 1, "K", 1) + 1), (get_base1(K, 1, "K", 1) + get_base1(K, 2, "K", 1))), stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list())), 
                        zeta_2, 
                        "assigning variable Zeta");
            current_statement_begin__ = 49;
            stan::math::assign(mu, multiply(X, B));
            current_statement_begin__ = 50;
            stan::math::assign(sigma, stan::math::exp(multiply(Z, Zeta)));
            current_statement_begin__ = 51;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 52;
                stan::model::assign(rho, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            ((2 * inv_logit(multiply(stan::model::rvalue(G, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "G"), gamma))) - 1), 
                            "assigning variable rho");
                current_statement_begin__ = 53;
                stan::model::assign(Omega, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            diag_matrix(rep_vector(1, 2)), 
                            "assigning variable Omega");
                current_statement_begin__ = 54;
                stan::model::assign(Omega, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_uni(1), stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()))), 
                            get_base1(rho, i, "rho", 1), 
                            "assigning variable Omega");
                current_statement_begin__ = 55;
                stan::model::assign(Omega, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_uni(2), stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()))), 
                            get_base1(rho, i, "rho", 1), 
                            "assigning variable Omega");
            }
            current_statement_begin__ = 59;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 60;
                stan::model::assign(Sigma, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            multiply(multiply(diag_matrix(transpose(stan::model::rvalue(sigma, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "sigma"))), get_base1(Omega, i, "Omega", 1)), diag_matrix(transpose(stan::model::rvalue(sigma, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "sigma")))), 
                            "assigning variable Sigma");
            }
            current_statement_begin__ = 64;
            lp_accum__.add(normal_log<propto__>(beta_1, 0, 5));
            current_statement_begin__ = 65;
            lp_accum__.add(normal_log<propto__>(beta_2, 0, 5));
            current_statement_begin__ = 66;
            lp_accum__.add(std_normal_log<propto__>(zeta_1));
            current_statement_begin__ = 67;
            lp_accum__.add(std_normal_log<propto__>(zeta_2));
            current_statement_begin__ = 68;
            lp_accum__.add(std_normal_log<propto__>(gamma));
            current_statement_begin__ = 71;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 72;
                lp_accum__.add(multi_normal_log<propto__>(transpose(stan::model::rvalue(y, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "y")), transpose(stan::model::rvalue(mu, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "mu")), get_base1(Sigma, i, "Sigma", 1)));
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta_1");
        names__.push_back("beta_2");
        names__.push_back("zeta_1");
        names__.push_back("zeta_2");
        names__.push_back("gamma");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(get_base1(P, 1, "P", 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(get_base1(P, 2, "P", 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(get_base1(K, 1, "K", 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(get_base1(K, 2, "K", 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(J);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_isoniche_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_1 = in__.vector_constrain(get_base1(P, 1, "P", 1));
        size_t beta_1_j_1_max__ = get_base1(P, 1, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_1_j_1_max__; ++j_1__) {
            vars__.push_back(beta_1(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_2 = in__.vector_constrain(get_base1(P, 2, "P", 1));
        size_t beta_2_j_1_max__ = get_base1(P, 2, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_2_j_1_max__; ++j_1__) {
            vars__.push_back(beta_2(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> zeta_1 = in__.vector_constrain(get_base1(K, 1, "K", 1));
        size_t zeta_1_j_1_max__ = get_base1(K, 1, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_1_j_1_max__; ++j_1__) {
            vars__.push_back(zeta_1(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> zeta_2 = in__.vector_constrain(get_base1(K, 2, "K", 1));
        size_t zeta_2_j_1_max__ = get_base1(K, 2, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_2_j_1_max__; ++j_1__) {
            vars__.push_back(zeta_2(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> gamma = in__.vector_constrain(J);
        size_t gamma_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < gamma_j_1_max__; ++j_1__) {
            vars__.push_back(gamma(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_isoniche";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_1_j_1_max__ = get_base1(P, 1, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_1_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_1" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_2_j_1_max__ = get_base1(P, 2, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_2_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_2" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t zeta_1_j_1_max__ = get_base1(K, 1, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_1_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "zeta_1" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t zeta_2_j_1_max__ = get_base1(K, 2, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_2_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "zeta_2" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t gamma_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < gamma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gamma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_1_j_1_max__ = get_base1(P, 1, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_1_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_1" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_2_j_1_max__ = get_base1(P, 2, "P", 1);
        for (size_t j_1__ = 0; j_1__ < beta_2_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_2" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t zeta_1_j_1_max__ = get_base1(K, 1, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_1_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "zeta_1" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t zeta_2_j_1_max__ = get_base1(K, 2, "K", 1);
        for (size_t j_1__ = 0; j_1__ < zeta_2_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "zeta_2" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t gamma_j_1_max__ = J;
        for (size_t j_1__ = 0; j_1__ < gamma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gamma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_isoniche_namespace::model_isoniche stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
