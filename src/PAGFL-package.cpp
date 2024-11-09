#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tuple>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

arma::vec bspline_basis(arma::vec &x, const arma::vec &knot_vec, const unsigned int &d, const unsigned int &indx)
{
    arma::vec B_i;
    arma::vec zero_vec = arma::zeros<arma::vec>(x.n_elem);
    if (d == 0)
    {
        // Edge case, i.e. either ones or zeros
        B_i = zero_vec;
        arma::uvec indx_vec = arma::find((x >= knot_vec[indx]) && (x < knot_vec[indx + 1]));
        B_i.elem(indx_vec).fill(1);
    }
    else
    {
        // Convex combinations of the two lower degree basis functions
        arma::vec term_1 = (knot_vec[d + indx] - knot_vec[indx]) == 0 ? zero_vec : (x - knot_vec[indx]) / (knot_vec[d + indx] - knot_vec[indx]);
        arma::vec term_2 = (knot_vec[d + indx + 1] - knot_vec[indx + 1]) == 0 ? zero_vec : (knot_vec[indx + d + 1] - x) / (knot_vec[indx + d + 1] - knot_vec[indx + 1]);
        // Recursive definition of the final basis
        arma::vec lower_term_1 = bspline_basis(x, knot_vec, (d - 1), indx);
        arma::vec lower_term_2 = bspline_basis(x, knot_vec, (d - 1), indx + 1);
        B_i = term_1 % lower_term_1 + term_2 % lower_term_2;
    }
    return B_i;
}

// [[Rcpp::export]]
arma::mat bspline_system(arma::vec &x, const unsigned int &d, const arma::vec &knots, bool intercept)
{
    // Define the complete vector of boundary and interior knots
    unsigned int n = knots.n_elem;
    arma::vec boundary_knots = {knots[0], knots[n - 1]};
    arma::vec interior_knots = knots.subvec(1, n - 2);
    arma::vec knot_vec = arma::join_cols(arma::join_cols(arma::ones<arma::vec>(d + 1) * boundary_knots[0], interior_knots), arma::ones<arma::vec>(d + 1) * boundary_knots[1]);
    // Initialize the B-spline system
    unsigned int M = interior_knots.n_elem + d + 1;
    arma::mat B(x.n_elem, M, arma::fill::zeros);
    // Compute the B-spline basis
    for (unsigned int l = 0; l < M; l++)
    {
        B.col(l) = bspline_basis(x, knot_vec, d, l);
    }
    // Fix edge cases
    arma::uvec boundary_indices = arma::find(x == boundary_knots[1]);
    if (boundary_indices.n_elem > 0)
    {
        for (unsigned int i = 0; i < boundary_indices.n_elem; i++)
        {
            B(boundary_indices[i], M - 1) = 1;
        }
    }
    // Remove the intercept if required
    if (!intercept)
    {
        B.shed_col(0);
    }
    return B;
}

// Project the regressors on the spline bases
arma::mat buildZ(const arma::mat &X, const arma::mat &B, const arma::uvec &t_index, const unsigned int &p)
{
    arma::mat Z(X.n_rows, B.n_cols * p, arma::fill::zeros);
    unsigned int t;
    for (unsigned int i = 0; i < t_index.n_elem; ++i)
    {
        t = t_index[i] - 1;
        Z.row(i) = arma::kron(X.row(i), B.row(t));
    }
    return Z;
}

// Ols via inverse chol decomposition
arma::vec ols_chol(arma::mat &XtX, arma::vec &Xty)
{
    arma::mat chol_decomp = arma::chol(XtX);
    arma::mat chol_inv = arma::inv(arma::trimatu(chol_decomp));
    arma::vec estim = chol_inv * chol_inv.t() * Xty;
    return estim;
}

arma::vec ols_naive(arma::mat &XtX, arma::vec &Xty)
{
    return arma::pinv(XtX) * Xty;
}

// Delete one observation per individual unit of the index vector
arma::uvec deleteOneObsperI(const arma::uvec &vec)
{
    arma::uvec unique_i = arma::unique(vec);
    arma::uvec result = vec;
    for (unsigned int i = 0; i < unique_i.n_elem; i++)
    {
        arma::uvec indices = arma::find(result == unique_i[i], 1);
        result.shed_row(indices[0]);
    }
    return result;
}

// Add one observation per individual unit to the index vector
arma::uvec addOneObsperI(const arma::uvec &vec)
{
    arma::uvec unique_values = arma::unique(vec);
    arma::uvec result = vec;
    arma::uvec indices;
    for (unsigned int i = 0; i < unique_values.n_elem; ++i)
    {
        indices = {unique_values[i]};
        result.insert_rows(result.n_elem, indices);
    }
    return arma::sort(result);
}

// Demean a matrix per individual
arma::mat demeanIndMat(arma::mat x, unsigned int N, arma::uvec i_index)
{
    arma::mat x_tilde = x;
    arma::uvec ind_seq;
    arma::mat x_red;
    arma::rowvec means;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        x_red = x.rows(ind_seq);
        means = arma::mean(x_red, 0);
        x_tilde.rows(ind_seq) = x_red - arma::repmat(means, x_red.n_rows, 1);
    }
    return x_tilde;
}

// Demean a vector per individual
// [[Rcpp::export]]
arma::vec demeanIndVec(arma::vec x, unsigned int N, arma::uvec i_index)
{
    arma::vec x_tilde = x;
    arma::uvec ind_seq, finite_inds;
    arma::vec x_red;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        x_red = x.elem(ind_seq);
        finite_inds = arma::find_finite(x_red);
        if (!finite_inds.is_empty())
        {
            x_tilde.elem(ind_seq) = x_red - arma::mean(x_red.elem(finite_inds));
        }
    }
    return x_tilde;
}

// Take column-wise differences of a matrix per individual
arma::mat fdIndMat(arma::mat x, unsigned int N, arma::uvec i_index)
{
    arma::mat x_tilde;
    arma::uvec ind_seq;
    arma::mat x_red, delta;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        x_red = x.rows(ind_seq);
        delta = arma::diff(x_red);
        x_tilde = arma::join_vert(x_tilde, delta);
    }
    return x_tilde;
}

// Take column-wise differences of a vector per individual
arma::vec fdIndVec(arma::vec x, unsigned int N, arma::uvec i_index)
{
    arma::mat x_tilde;
    arma::uvec ind_seq;
    arma::vec x_red, delta;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        x_red = x.rows(ind_seq);
        delta = arma::diff(x_red);
        x_tilde = arma::join_vert(x_tilde, delta);
    }
    return x_tilde;
}

// Obtain the penalty matrix
arma::sp_mat buildLambda(const unsigned int &p, const unsigned int &N)
{
    std::vector<arma::mat> D_vec;
    // Create a complete fusion penalty matrix D
    for (unsigned int x = 0; x < (N - 1); x++)
    {
        arma::mat temp = arma::zeros<arma::mat>(N - 1 - x, x);
        arma::vec ones = arma::ones<arma::vec>(N - 1 - x);
        arma::mat diag = -arma::eye<arma::mat>(N - 1 - x, N - 1 - x);
        temp = arma::join_horiz(temp, ones);
        temp = arma::join_horiz(temp, diag);
        D_vec.push_back(temp);
    }
    // Initialize D with the first matrix in D_vec
    arma::mat D = D_vec[0];
    // Join the rest of the matrices in D_vec vertically
    for (unsigned int i = 1; i < D_vec.size(); i++)
    {
        D = arma::join_vert(D, D_vec[i]);
    }
    // Create a complete group penalty matrix Lambda
    arma::sp_mat Lambda = sp_mat(arma::kron(D, arma::eye<arma::mat>(p, p)));
    return Lambda;
}

// Obtain adaptive penalty weights
arma::vec getOmega(const arma::vec &delta, const double &kappa, const unsigned int &N, const unsigned int &p)
{
    int n = N * (N - 1) * p / 2;
    // Build a matrix that holds all the indices of delta that belong to the same individual per row
    arma::mat ind_mat = arma::linspace<arma::mat>(1, n, n);
    ind_mat.reshape(p, n / p);
    ind_mat = ind_mat.t();
    arma::vec omega_single(n / p);
    // Apply the (to the power of -kappa) group-wise L2 norm to obtain the penalty weights
    for (unsigned int i = 0; i < n / p; i++)
    {
        arma::uvec ind = arma::conv_to<arma::uvec>::from(ind_mat.row(i)) - 1;
        float norm_val = arma::norm(delta.elem(ind), "fro");
        omega_single(i) = std::pow(norm_val, -kappa);
    }
    // Expand the vector to make later computations easier
    arma::vec omega = arma::vectorise(arma::repmat(omega_single, 1, p).t());
    return omega;
}

// Update the slope parameter vector
arma::vec getBeta(arma::vec invXovY, arma::mat invXov, arma::sp_mat VarLambdat, double varrho, arma::vec v, arma::vec delta)
{
    // See Mehrabani (2023, sec 5.1 step 2a/ 5.2 step 2a)
    arma::vec lagr_vev = VarLambdat * (delta - v / varrho);
    arma::vec beta = invXovY + invXov * lagr_vev;
    return beta;
}

// Apply the soft thresholding operator (see Mehrabani (2023, eq. 5.1))
arma::vec softThreshold(const arma::uvec &ind, const arma::vec &a, const arma::vec &b)
{
    // Extract the elements for the group of parameters in question
    arma::vec a_red = a.elem(ind);
    arma::vec b_red = b.elem(ind);
    // apply the group-wise soft thresholding operator
    arma::vec c = 1 - b_red / arma::norm(a_red, "fro");
    c = arma::clamp(c, 0, arma::datum::inf);
    arma::vec delta = c % a_red;
    return delta;
}

struct DeltaWorker : public Worker
{
    const arma::mat &ind_mat;
    const arma::vec &xi;
    const arma::vec &ada_weights;
    arma::vec &delta;
    unsigned int p;

    DeltaWorker(const arma::mat &ind_mat, const arma::vec &xi, const arma::vec &ada_weights, arma::vec &delta, unsigned int p)
        : ind_mat(ind_mat), xi(xi), ada_weights(ada_weights), delta(delta), p(p) {}

    void operator()(std::size_t begin, std::size_t end)
    {
        for (unsigned int i = begin; i < end; i++)
        {
            arma::uvec ind = arma::conv_to<arma::uvec>::from(ind_mat.row(i)) - 1;
            delta.subvec(i * p, (i + 1) * p - 1) = softThreshold(ind, xi, ada_weights);
        }
    }
};

// Update the parameter difference vector (see Mehrabani (2023, sec 5.1 step 2b))
arma::vec getDelta(const arma::vec &ada_weights, const arma::vec &beta, const arma::vec &v, const arma::sp_mat &Lambda, const double &varrho, const unsigned int &N, const unsigned int &p, const bool &parallel)
{
    // Sum of parameter differences and the Lagrangian parameters
    arma::vec xi = Lambda * beta + v / varrho;
    int n = N * (N - 1) * p / 2;
    // Matrix indicating the elements that belong to the same individual
    arma::mat ind_mat = arma::linspace<arma::mat>(1, n, n);
    ind_mat.reshape(p, n / p);
    ind_mat = ind_mat.t();
    arma::vec delta(n);
    // Vector of complete fusion parameter differences
    if (parallel)
    {
        DeltaWorker deltaWorker(ind_mat, xi, ada_weights, delta, p);
        parallelFor(0, n / p, deltaWorker);
    }
    else
    {
        for (unsigned int i = 0; i < n / p; i++)
        {
            arma::uvec ind = arma::conv_to<arma::uvec>::from(ind_mat.row(i)) - 1;
            delta.subvec(i * p, (i + 1) * p - 1) = softThreshold(ind, xi, ada_weights);
        }
    }
    return delta;
}

// Checks for convergence
bool stoppingCrit(const arma::vec &resid, const double &tol)
{
    // Check if the estimates have sufficiently converged
    bool result = arma::norm(resid, "fro") < tol;
    return result;
}

// Constructs  a block diagonal matrix from a vector of matrices
arma::sp_mat buildBlockDiag(const std::vector<arma::mat> &XList)
{
    // Calculate the total number of rows and columns for the block diagonal matrix
    int total_rows = 0;
    int total_cols = 0;
    for (const auto &X : XList)
    {
        total_rows += X.n_rows;
        total_cols += X.n_cols;
    }
    // Initialize the block diagonal matrix with zeros
    arma::mat block_diag(total_rows, total_cols, arma::fill::zeros);
    // Place each matrix from XList on the diagonal of block_diag
    int row_start = 0;
    int col_start = 0;
    for (const auto &X : XList)
    {
        block_diag.submat(row_start, col_start, row_start + X.n_rows - 1, col_start + X.n_cols - 1) = X;
        row_start += X.n_rows;
        col_start += X.n_cols;
    }
    return sp_mat(block_diag);
}

// Constructs a NT x Kp block predictor matrix
arma::sp_mat buildDiagX_block(const arma::mat &X, const unsigned int &N, arma::uvec &i_index, const arma::uvec &groups)
{
    // Construct a vector with one element per cross-sectional individual
    std::vector<arma::mat> XMat_vec(N);
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        XMat_vec[i] = X.rows(ind_seq);
    }
    // Impose a certain grouping
    unsigned int nGroups = arma::max(groups);
    std::vector<arma::mat> XMat_tilde_vec(nGroups);
    for (unsigned int i = 0; i < nGroups; i++)
    {
        arma::mat groupMatrix;
        for (unsigned int j = 0; j < N; j++)
        {
            if (groups[j] == i + 1)
            {
                groupMatrix = arma::join_cols(groupMatrix, XMat_vec[j]);
            }
        }
        XMat_tilde_vec[i] = groupMatrix;
    }
    // Create the final block matrix
    arma::sp_mat X_tilde = buildBlockDiag(XMat_tilde_vec);
    return X_tilde;
}

// [[Rcpp::export]]
arma::mat buildDiagX_block_dense(const arma::mat &X, const unsigned int &N, arma::uvec &i_index, const arma::uvec &groups)
{
    // Construct a vector with one element per cross-sectional individual
    std::vector<arma::mat> XMat_vec(N);
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        XMat_vec[i] = X.rows(ind_seq);
    }
    // Impose a certain grouping
    unsigned int nGroups = arma::max(groups);
    std::vector<arma::mat> XMat_tilde_vec(nGroups);
    for (unsigned int i = 0; i < nGroups; i++)
    {
        arma::mat groupMatrix;
        for (unsigned int j = 0; j < N; j++)
        {
            if (groups[j] == i + 1)
            {
                groupMatrix = arma::join_cols(groupMatrix, XMat_vec[j]);
            }
        }
        XMat_tilde_vec[i] = groupMatrix;
    }
    // Create the final block matrix
    arma::mat X_tilde = arma::mat(buildBlockDiag(XMat_tilde_vec));
    return X_tilde;
}

struct BetaWorker : public Worker
{
    const arma::uvec &groups;
    const std::vector<arma::mat> &XMat_vec;
    const std::vector<arma::vec> &y_vec;
    const bool &robust;
    arma::mat &beta_mat;
    std::vector<arma::mat> &XMat_tilde_vec;
    std::vector<arma::mat> &XtX_tilde_vec;

    BetaWorker(const arma::uvec &groups, const std::vector<arma::mat> &XMat_vec, const std::vector<arma::vec> &y_vec, const bool &robust, arma::mat &beta_mat, std::vector<arma::mat> &XMat_tilde_vec, std::vector<arma::mat> &XtX_tilde_vec)
        : groups(groups), XMat_vec(XMat_vec), y_vec(y_vec), robust(robust), beta_mat(beta_mat), XMat_tilde_vec(XMat_tilde_vec), XtX_tilde_vec(XtX_tilde_vec) {}

    void operator()(std::size_t begin, std::size_t end)
    {
        for (unsigned int i = begin; i < end; i++)
        {
            arma::mat groupMatrix;
            arma::vec group_y;
            for (unsigned int j = 0; j < groups.n_elem; j++)
            {
                if (groups[j] == i + 1)
                {
                    groupMatrix = arma::join_cols(groupMatrix, XMat_vec[j]);
                    group_y = arma::join_cols(group_y, y_vec[j]);
                }
            }
            XMat_tilde_vec[i] = groupMatrix;
            arma::mat Xt_tilde = groupMatrix.t();
            arma::mat XtX_tilde = Xt_tilde * groupMatrix;
            XtX_tilde_vec[i] = XtX_tilde;
            arma::vec Xty_tilde = Xt_tilde * group_y;
            arma::vec beta;
            if (robust)
            {
                beta = ols_naive(XtX_tilde, Xty_tilde);
            }
            else
            {
                beta = ols_chol(XtX_tilde, Xty_tilde);
            }
            beta_mat.row(i) = beta.t();
        }
    }
};

// Constructs a NT x Kp block predictor matrix, inverse and cross-product
std::tuple<arma::sp_mat, arma::sp_mat, arma::vec> buildDiagX(const arma::mat &X, const arma::vec &y, const unsigned int &N, arma::uvec &i_index, const arma::uvec &groups, const bool &robust, const bool &parallel)
{
    // Construct a vector with one element per cross-sectional individual
    std::vector<arma::mat> XMat_vec(N);
    std::vector<arma::vec> y_vec(N);
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        XMat_vec[i] = X.rows(ind_seq);
        y_vec[i] = y.elem(ind_seq);
    }
    // Iterate over a certain grouping
    unsigned int nGroups = arma::max(groups);
    std::vector<arma::mat> XMat_tilde_vec(nGroups);
    std::vector<arma::mat> XtX_tilde_vec(nGroups);
    arma::mat beta_mat = arma::zeros<arma::mat>(N, X.n_cols);
    if (parallel)
    {
        BetaWorker betaWorker(groups, XMat_vec, y_vec, robust, beta_mat, XMat_tilde_vec, XtX_tilde_vec);
        parallelFor(0, nGroups, betaWorker);
    }
    else
    {
        for (unsigned int i = 0; i < nGroups; i++)
        {
            arma::mat groupMatrix, XtX_tilde, Xt_tilde;
            arma::vec beta, group_y, Xty_tilde;
            for (unsigned int j = 0; j < N; j++)
            {
                if (groups[j] == i + 1)
                {
                    groupMatrix = arma::join_cols(groupMatrix, XMat_vec[j]);
                    group_y = arma::join_cols(group_y, y_vec[j]);
                }
            }
            XMat_tilde_vec[i] = groupMatrix;
            Xt_tilde = groupMatrix.t();
            XtX_tilde = Xt_tilde * groupMatrix;
            XtX_tilde_vec[i] = XtX_tilde;
            Xty_tilde = Xt_tilde * group_y;
            if (robust)
            {
                beta = ols_naive(XtX_tilde, Xty_tilde);
            }
            else
            {
                beta = ols_chol(XtX_tilde, Xty_tilde);
            }
            beta_mat.row(i) = beta.t();
        }
    }

    // Create the final block matrices
    arma::sp_mat X_tilde = buildBlockDiag(XMat_tilde_vec);
    arma::sp_mat XtX = buildBlockDiag(XtX_tilde_vec);
    std::vector<arma::mat> X_tilde_vec(3);
    X_tilde_vec[0] = X_tilde;
    X_tilde_vec[1] = XtX;
    return std::make_tuple(X_tilde, XtX, arma::vectorise(beta_mat.t()));
}

// inverts the error var-cov matrix for each group
arma::sp_mat invertV(const arma::mat &V, const unsigned int &q)
{
    unsigned int N = V.n_rows / q;
    std::vector<arma::mat> W_list(N);
    arma::mat V_red;
    // V is a block matrix of group-wise error var-cov matrices
    for (unsigned int i = 0; i < N; i++)
    {
        // Pick and invert each submatrix individually
        V_red = V.submat(i * q, i * q, (i + 1) * q - 1, (i + 1) * q - 1);
        W_list[i] = arma::pinv(V_red);
    }
    // Put everything back together
    arma::sp_mat W = buildBlockDiag(W_list);
    return W;
}

// Compute the weight matrix for the GMM
arma::sp_mat getW(const arma::mat &X_block, arma::mat &Z_block, const arma::vec &y, const unsigned int &q)
{
    // Initial estimate with an identity weight matrix
    arma::mat Xt = X_block.t() * Z_block * Z_block.t();
    arma::mat XtX = Xt * X_block;
    arma::vec Xty = Xt * y;
    arma::vec beta = arma::inv(XtX + .05 / sqrt(y.n_elem)) * Xty;
    // Obtain the (group-wise) var-cov matrix of the residuals
    arma::vec u = y - X_block * beta;
    arma::mat V = Z_block.t() * u * u.t() * Z_block;
    // Use the inverse as the subsequent weight matrix (tossing all matrices off the diagonal)
    arma::sp_mat W = invertV(V, q);
    return W;
}

struct AlphaWorker : public Worker
{
    const arma::uvec &groups;
    const std::vector<arma::mat> &XMat_vec;
    const std::vector<arma::vec> &y_vec;
    const bool &robust;
    arma::mat &alpha_mat;

    AlphaWorker(const arma::uvec &groups, const std::vector<arma::mat> &XMat_vec, const std::vector<arma::vec> &y_vec, const bool &robust, arma::mat &alpha_mat)
        : groups(groups), XMat_vec(XMat_vec), y_vec(y_vec), robust(robust), alpha_mat(alpha_mat) {}

    void operator()(std::size_t begin, std::size_t end)
    {
        for (unsigned int k = begin; k < end; k++)
        {
            arma::mat groupX;
            arma::vec groupy;
            for (unsigned int j = 0; j < groups.n_elem; j++)
            {
                if (groups[j] == k + 1)
                {
                    groupX = arma::join_cols(groupX, XMat_vec[j]);
                    groupy = arma::join_cols(groupy, y_vec[j]);
                }
            }
            arma::mat groupXt = groupX.t();
            arma::mat groupXtX = groupXt * groupX;
            arma::vec groupXty = groupXt * groupy;
            arma::vec alpha_group;
            if (robust)
            {
                alpha_group = ols_naive(groupXtX, groupXty);
            }
            else
            {
                alpha_group = ols_chol(groupXtX, groupXty);
            }
            alpha_mat.row(k) = alpha_group.t();
        }
    }
};

// Group-wise OLS estimation
arma::mat getGroupwiseOLS(const arma::vec &y, const arma::mat &X, const unsigned int &N, arma::uvec &i_index, const arma::uvec &groups, const unsigned int &p, const bool robust, const bool &parallel)
{
    // Construct a vector with one element per cross-sectional individual
    std::vector<arma::mat> XMat_vec(N);
    std::vector<arma::vec> y_vec(N);
    // Construct a matrix with one row per group
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        XMat_vec[i] = X.rows(ind_seq);
        y_vec[i] = y.elem(ind_seq);
    }

    // Impose a certain grouping
    unsigned int nGroups = arma::max(groups);
    arma::mat alpha_mat = arma::zeros<arma::mat>(nGroups, p);
    if (parallel)
    {
        AlphaWorker alphaWorker(groups, XMat_vec, y_vec, robust, alpha_mat);
        parallelFor(0, nGroups, alphaWorker);
    }
    else
    {
        for (unsigned int k = 0; k < nGroups; k++)
        {
            arma::mat groupX, groupXt, groupXtX;
            arma::vec groupy, alpha_group, groupXty;
            // Obtain which individuals are part of group k
            for (unsigned int j = 0; j < N; j++)
            {
                if (groups[j] == k + 1)
                {
                    groupX = arma::join_cols(groupX, XMat_vec[j]);
                    groupy = arma::join_cols(groupy, y_vec[j]);
                }
            }
            groupXt = groupX.t();
            groupXtX = groupXt * groupX;
            groupXty = groupXt * groupy;
            if (robust)
            {
                alpha_group = ols_naive(groupXtX, groupXty);
            }
            else
            {
                alpha_group = ols_chol(groupXtX, groupXty);
            }
            alpha_mat.row(k) = alpha_group.t();
        }
    }

    return alpha_mat;
}

// Computes the post-lasso estimates based on an estimated grouping
arma::mat getAlpha(const arma::mat &X, const arma::vec &y, arma::mat &Z, const std::string &method, const unsigned int &N, arma::uvec &i_index, const unsigned int &p, const arma::uvec &groups_hat, const bool robust, const bool &parallel)
{
    // Compute the post-lasso estimates
    arma::vec alpha, y_tilde;
    std::vector<arma::vec> y_vec(N);
    std::vector<std::pair<int, int>> group_index_pairs(N);
    arma::uvec ind_seq;
    arma::mat groupZ, W, groupX, groupXt, alpha_mat;
    if (method == "PGMM")
    {
        // Order the individuals in the dependent variable vector according to the groups
        // Create a vector of pairs where the first element of each pair is the group and the second element is the original index
        for (unsigned int i = 0; i < N; i++)
        {
            ind_seq = arma::find(i_index == i + 1);
            y_vec[i] = y.elem(ind_seq);
            group_index_pairs[i] = std::make_pair(groups_hat[i], i);
        }
        // Sort the vector of pairs based on the group
        std::sort(group_index_pairs.begin(), group_index_pairs.end());
        // Create a new y vector where the entries are sorted based on the group structure
        for (unsigned int i = 0; i < group_index_pairs.size(); ++i)
        {
            y_tilde = arma::join_cols(y_tilde, y_vec[group_index_pairs[i].second]);
        }
        // Construct the regressors matrices
        groupX = buildDiagX_block(X, N, i_index, groups_hat);
        groupZ = buildDiagX_block(Z, N, i_index, groups_hat);
        // W = getW(groupX, groupZ, y_tilde, Z.n_cols);
        W = arma::eye<arma::mat>(groupZ.n_cols, groupZ.n_cols);
        groupXt = groupX.t() * groupZ * W * groupZ.t();
        // Compute the post-lasso estimates
        alpha = arma::inv(groupXt * groupX + .05 / sqrt(y.n_elem)) * groupXt * y_tilde;
        alpha_mat = arma::reshape(alpha, p, alpha.n_elem / p).t();
    }
    else
    {
        alpha_mat = getGroupwiseOLS(y, X, N, i_index, groups_hat, p, robust, parallel);
    }
    return alpha_mat;
}

// Splits the observed time periods of a NT x p matrix in half per individual
arma::mat splitMatInHalf(const arma::mat &X, const unsigned int &N, arma::uvec &i_index, const unsigned int &half)
{
    arma::mat X_half;
    arma::uvec ind_seq, ind_seq_half;
    unsigned int half_ind, start, end;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        half_ind = ind_seq.n_elem / 2;
        start = (half == 1) ? 0 : half_ind;
        end = start + half_ind - 1;
        ind_seq_half = ind_seq.subvec(start, end);
        X_half = arma::join_cols(X_half, X.rows(ind_seq_half));
    }
    return X_half;
}

// Splits the observed time periods of a NT vector in half per individual
arma::vec splitVecInHalf(const arma::vec &X, const unsigned int &N, arma::uvec &i_index, const unsigned int &half)
{
    arma::vec X_half;
    arma::uvec ind_seq, ind_seq_half;
    unsigned int start, end;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        unsigned int half_ind = ind_seq.n_elem / 2;
        start = (half == 1) ? 0 : half_ind;
        end = start + half_ind - 1;
        ind_seq_half = ind_seq.subvec(start, end);
        X_half = arma::join_cols(X_half, X.elem(ind_seq_half));
    }
    return X_half;
}

// Obtain the first half of the index vector per individual
arma::uvec splitIndexInHalf(arma::uvec &i_index, const unsigned int &half)
{
    arma::uvec unique_i = arma::unique(i_index);
    arma::uvec i_index_half, ind_seq, half_ind_seq;
    unsigned int start, end, half_ind;
    for (unsigned int i = 0; i < unique_i.n_elem; i++)
    {
        ind_seq = arma::find(i_index == unique_i[i]);
        half_ind = ind_seq.n_elem / 2;
        start = (half == 1) ? 0 : half_ind;
        end = start + half_ind - 1;
        half_ind_seq = ind_seq.subvec(start, end);
        i_index_half = arma::join_cols(i_index_half, i_index.elem(half_ind_seq));
    }
    return i_index_half;
}

// Drops either the first or the last observed time period of a NT x p matrix per individual
arma::mat deleteObsMat(const arma::mat &X, const unsigned int &N, arma::uvec &i_index, const bool first)
{
    arma::mat x_red, X_red;
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        if (first)
        {
            ind_seq.shed_row(0);
        }
        else
        {
            ind_seq.shed_row(ind_seq.n_elem - 1);
        }
        x_red = X.rows(ind_seq);
        X_red = arma::join_cols(X_red, x_red);
    }
    return X_red;
}

// Net out fixed effects
std::vector<arma::mat> netFE(arma::vec &y, arma::mat &X, const std::string &method, const unsigned int &N, arma::uvec &i_index)
{
    arma::mat (*netFEMat)(arma::mat, unsigned int, arma::uvec);
    arma::vec (*netFEVec)(arma::vec, unsigned int, arma::uvec);
    if (method == "PLS")
    {
        netFEMat = &demeanIndMat;
        netFEVec = &demeanIndVec;
    }
    else
    {
        netFEMat = &fdIndMat;
        netFEVec = &fdIndVec;
    }

    arma::mat y_tilde = netFEVec(y, N, i_index);
    arma::mat X_tilde = netFEMat(X, N, i_index);
    std::vector<arma::mat> data(2);
    data[0] = y_tilde;
    data[1] = X_tilde;
    return data;
}

// Deletes one individual observation from the cross-sectional index vector in case the number of observations per individual is uneven
arma::uvec getEvenT_index(arma::uvec &i_index, const unsigned int N)
{
    arma::uvec ind_seq, i_index_trunc;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        if (ind_seq.n_elem % 2 != 0)
        {
            ind_seq.shed_row(ind_seq.n_elem - 1);
            i_index_trunc = arma::join_cols(i_index_trunc, i_index.elem(ind_seq));
        }
        else
        {
            i_index_trunc = arma::join_cols(i_index_trunc, i_index.elem(ind_seq));
        }
    }
    return i_index_trunc;
}

// Deletes one individual observation from the a vector in case the number of observations per individual is uneven
arma::vec getEvenT_vec(const arma::vec &X, const unsigned int &N, arma::uvec i_index)
{
    arma::vec x_red, X_red;
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        x_red = X.elem(ind_seq);
        if (ind_seq.n_elem % 2 != 0)
        {
            x_red.shed_row(x_red.n_elem - 1);
        }
        X_red = arma::join_cols(X_red, x_red);
    }
    return X_red;
}

// Deletes one individual observation from the a matrix in case the number of observations per individual is uneven
arma::mat getEvenT_mat(const arma::mat &X, const unsigned int &N, arma::uvec i_index)
{
    arma::mat x_red, X_red;
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        x_red = X.rows(ind_seq);
        if (ind_seq.n_elem % 2 != 0)
        {
            x_red.shed_row(x_red.n_rows - 1);
        }
        X_red = arma::join_cols(X_red, x_red);
    }
    return X_red;
}

// Compute the Split-panel Jackknife bias corrected estimates
arma::mat spjCorrec(const arma::mat &alpha_mat, arma::mat &X, arma::vec &y, arma::mat &Z, const unsigned int &N, arma::uvec &i_index, const unsigned int &p, const arma::uvec &groups_hat, const std::string &method, const bool &parallel)
{
    // Discard the last time period in case n_periods is not even
    arma::uvec i_index_trunc = getEvenT_index(i_index, N);
    arma::mat X_trunc = getEvenT_mat(X, N, i_index);
    arma::vec y_trunc = getEvenT_vec(y, N, i_index);
    arma::mat Z_trunc;
    if (method == "PGMM")
    {
        Z_trunc = getEvenT_mat(Z, N, i_index);
    }

    // Pull the number of groups
    unsigned int K_hat = arma::max(groups_hat);
    // Initialize
    arma::cube alpha_cube(K_hat, p, 2);
    arma::mat X_half, alpha_half, X_half_tilde, Z_half, Z_half_tilde;
    arma::vec y_half, y_half_tilde;
    arma::uvec i_index_tilde, i_index_half;
    if (method == "PGMM")
    {
        i_index_tilde = deleteOneObsperI(i_index_trunc);
    }
    else
    {
        i_index_tilde = i_index_trunc;
    }
    // Compute the bias correction
    for (unsigned int i = 1; i < 3; i++)
    {
        // Split y and X in half - per individual
        X_half = splitMatInHalf(X_trunc, N, i_index_trunc, i);
        y_half = splitVecInHalf(y_trunc, N, i_index_trunc, i);
        i_index_half = splitIndexInHalf(i_index_trunc, i);
        // Net out fixed effects
        std::vector<arma::mat> data = netFE(y_half, X_half, method, N, i_index_half);
        y_half_tilde = data[0];
        X_half_tilde = data[1];
        if (method == "PGMM")
        {
            Z_half = splitMatInHalf(Z_trunc, N, i_index_trunc, i);
            Z_half_tilde = deleteObsMat(Z_half, N, i_index_half, TRUE);
            i_index_half = deleteOneObsperI(i_index_half);
        }
        // Estimate
        alpha_half = getAlpha(X_half_tilde, y_half_tilde, Z_half_tilde, method, N, i_index_half, p, groups_hat, FALSE, parallel);
        // Store
        alpha_cube.slice(i - 1) = alpha_half;
    }
    // Take the mean across panel halves
    arma::mat alpha_mean = arma::mean(alpha_cube, 2);
    // Apply the correction
    arma::mat alpha_correc = 2 * alpha_mat - alpha_mean;
    return alpha_correc;
}

// Create a vector of all pairs of individuals that are part of the same group
arma::uvec getGroupPairs(int indx, int N)
{
    arma::vec group_indx = arma::cumsum(arma::linspace(N - 1, 1, N - 1));
    // Identify the first individual
    unsigned int first = arma::sum(indx > group_indx) + 1;
    // Identify the second individual
    unsigned int second = indx - std::max(first > 1 ? group_indx(first - 2) : 0, 0.0) + first;
    return arma::uvec({first, second});
}

// Relabel a vector of an estimated grouping to be increasing and consecutive in numeric labels
arma::uvec relabelGroups(const arma::uvec &groups)
{
    std::map<int, int> labels;
    int nextLabel = 1;
    arma::uvec groups_relabel(groups.n_elem);
    for (unsigned int i = 0; i < groups.n_elem; ++i)
    {
        if (groups[i] != 0)
        {
            if (labels.count(groups[i]) == 0)
            {
                labels[groups[i]] = nextLabel++;
            }
            groups_relabel[i] = labels[groups[i]];
        }
    }
    return groups_relabel;
}

// Assigns individuals to groups based on the estimated parameter differences
arma::uvec getGroups(const arma::vec &beta, const arma::vec &y, const arma::mat &X, const arma::sp_mat &Lambda, const unsigned int &p, const unsigned int &N, const double &tol)
{
    // Assign the groups based on which estimates are equal enough
    unsigned int n = N * (N - 1) / 2;
    arma::mat delta = arma::reshape(Lambda * beta, p, n).t();
    arma::vec deltaNorm(n);
    for (unsigned int i = 0; i < n; i++)
    {
        deltaNorm[i] = arma::norm(delta.row(i), "fro");
    }
    arma::uvec identicalUnits = arma::find(deltaNorm < tol);
    arma::uvec groups(N);
    // Construct a vector with all individuals per entry belonging to one group
    arma::umat identicalUnits_mat(2, identicalUnits.n_elem, fill::zeros);
    arma::uvec first_params;
    if (identicalUnits.n_elem == 0)
    {
        groups = arma::linspace<arma::uvec>(1, N, N);
    }
    else
    {
        // Matrix with pairs of individuals that are part of the same group
        for (unsigned int i = 0; i < identicalUnits.n_elem; i++)
        {
            identicalUnits_mat.col(i) = getGroupPairs(identicalUnits[i] + 1, N);
        }
        first_params = unique(identicalUnits_mat.row(0).t());
    }
    std::vector<arma::uvec> groupList(first_params.n_elem);
    // Collapse into a single vector of group memberships
    for (unsigned int i = 0; i < first_params.n_elem; i++)
    {
        arma::uvec second_params_indices = arma::find(identicalUnits_mat.row(0) == first_params[i]);
        arma::uvec second_params(second_params_indices.n_elem);
        for (unsigned int i = 0; i < second_params_indices.n_elem; ++i)
        {
            second_params[i] = identicalUnits_mat(1, second_params_indices[i]);
        }
        arma::uvec group = arma::join_cols(arma::uvec({first_params[i]}), second_params);
        groupList[i] = group;
    }
    // Force transitivity among the groups by pooling groups with at least one shared individual
    unsigned int groupListSize = groupList.size();
    if (groupListSize > 1)
    {
        for (unsigned int i = 0; i < groupListSize; i++)
        {
            for (unsigned int j = fmin(i + 1, groupListSize); j < groupListSize; j++)
            {
                arma::uvec intersect = arma::intersect(groupList[i], groupList[j]);
                if (intersect.n_elem > 0)
                {
                    groupList[i] = arma::unique(arma::join_cols(groupList[i], groupList[j]));
                    groupList[j] = arma::uvec();
                }
            }
        }
    }
    // Remove empty groups
    groupList.erase(std::remove_if(groupList.begin(), groupList.end(), [](const arma::uvec &v)
                                   { return v.is_empty(); }),
                    groupList.end());
    // Compute a vector of group adherences
    arma::uvec groups_hat_raw(N, arma::fill::zeros);
    for (unsigned int i = 0; i < groupList.size(); ++i)
    {
        for (unsigned int j = 0; j < groupList[i].n_elem; ++j)
        {
            groups_hat_raw[groupList[i][j] - 1] = i + 1;
        }
    }
    // Total number of groups
    int K_hat = arma::max(groups_hat_raw);
    // Assign single individuals to their own group
    for (unsigned int i = 0; i < N; i++)
    {
        if (groups_hat_raw[i] == 0)
        {
            groups_hat_raw[i] = ++K_hat;
        }
    }
    // Relabel the groups to be consecutive
    arma::uvec groups_hat = relabelGroups(groups_hat_raw);
    return groups_hat;
}

// Obtain the sum of squared residuals of some post-lasso estimates
float getSSQ(const arma::vec &y, const arma::mat &X, const arma::mat alpha, const arma::uvec groups, const unsigned int &N, arma::uvec &i_index)
{
    arma::uvec individuals = arma::linspace<arma::uvec>(1, N, N);
    arma::mat X_i;
    arma::vec y_i;
    arma::rowvec alpha_vec;
    arma::uvec ind_seq;
    // Compute the ssq
    float ssq_sum = 0;
    for (unsigned int i = 0; i < N; ++i)
    {
        alpha_vec = alpha.row(groups[i] - 1);
        ind_seq = arma::find(i_index == i + 1);
        y_i = y.elem(ind_seq);
        X_i = X.rows(ind_seq);
        ssq_sum += arma::sum(arma::pow(y_i - X_i * alpha_vec.t(), 2));
    }
    float ssq = ssq_sum / y.n_elem;
    return ssq;
}

// Allocate individuals that are part of trivial groups to the best-fitting non-trivial group
arma::uvec mergeTrivialGroups(arma::uvec &groups_hat, const arma::vec &y, const arma::mat &X, arma::mat &Z, const std::string &method, const float &min_group_frac, const unsigned int &N, arma::uvec &i_index, const unsigned int &p, const bool robust, const bool &parallel)
{
    int limit = std::floor(min_group_frac * N);
    // Identify trivial and non-trivial groups
    arma::uvec unique_groups = arma::unique(groups_hat);
    arma::uvec groupCardinality = arma::hist(groups_hat, unique_groups);
    arma::uvec trivialGroups_ind = arma::find(groupCardinality < limit);
    arma::uvec nonTrivialGroups_ind = arma::find(groupCardinality >= limit);
    arma::uvec trivialGroups = unique_groups.elem(trivialGroups_ind);
    arma::uvec nonTrivialGroups = unique_groups.elem(nonTrivialGroups_ind);
    arma::uvec trivialInd;
    if (trivialGroups.n_elem == 0 || nonTrivialGroups.n_elem == 0)
    {
        return groups_hat;
    }
    else
    {
        // Get the indices of individuals part of trivial groups
        for (unsigned int i = 0; i < trivialGroups.n_elem; ++i)
        {
            trivialInd = arma::join_cols(trivialInd, arma::find(groups_hat == trivialGroups(i)));
        }
        // Iterate over all individuals part of trivial groups
        arma::uvec groups_hat_tmp, groups_hat_tmp_label;
        for (unsigned int i = 0; i < trivialInd.n_elem; ++i)
        {
            arma::vec ssq_vec(nonTrivialGroups.n_elem);
            // Find them a new home by iterating over all non-trivial group
            for (unsigned int j = 0; j < nonTrivialGroups.n_elem; ++j)
            {
                groups_hat_tmp = groups_hat;
                groups_hat_tmp(trivialInd(i)) = nonTrivialGroups(j);
                groups_hat_tmp_label = relabelGroups(groups_hat_tmp);
                arma::mat alpha_tmp = getAlpha(X, y, Z, method, N, i_index, p, groups_hat_tmp_label, robust, parallel);
                ssq_vec(j) = getSSQ(y, X, alpha_tmp, groups_hat_tmp_label, N, i_index);
            }
            // Select the non-trivial group with the best fit
            groups_hat(trivialInd(i)) = nonTrivialGroups(arma::index_min(ssq_vec));
        }
        return relabelGroups(groups_hat);
    }
}

// PAGFL routine
Rcpp::List pagfl_algo(arma::vec &y, arma::vec &y_tilde, arma::mat &X, arma::mat &X_tilde, arma::vec &invXcovY, arma::mat &invXcov, const arma::sp_mat &VarLambdat, const arma::sp_mat &Lambda, const std::string &method, arma::mat &Z, arma::mat &Z_tilde, const arma::vec &delta_ini, const arma::vec &omega, arma::vec &v_old, arma::uvec i_index, const arma::uvec &t_index, const unsigned int &N, const unsigned int &n, const unsigned int &p, const unsigned int &q, unsigned int &n_periods, const bool &bias_correc, const double &lambda, const double &min_group_frac, const unsigned int &max_iter, const double &tol_convergence, const double &tol_group, const double &varrho, const bool &parallel, const bool &verbose, const unsigned int &lambda_num, const unsigned int &n_lambda)
{

    //------------------------------//
    // Preliminaries                //
    //------------------------------//

    float lambda_star;
    if (method == "PLS")
    {
        lambda_star = n_periods * lambda / (2 * N);
    }
    else
    {
        lambda_star = std::pow(n_periods, 2) * lambda / (2 * N);
    }
    arma::vec ada_weights = omega * lambda_star / varrho;
    arma::vec resid, v_new, beta;
    arma::vec delta = delta_ini;

    //------------------------------//
    // Run the algorithm            //
    //------------------------------//

    unsigned int iter = 0;
    for (unsigned int i = 0; i < max_iter; i++)
    {
        if (verbose)
        {
          Rcout << "\r" << "Lambda: " << lambda << " (" << lambda_num << "/" << n_lambda << ")" << " - Iteration: " << i + 1 << "/" << max_iter << std::flush;
        }
        // Update the ...
        // parameter estimates (sec. 5.1/ 5.2 step 2a)
        beta = getBeta(invXcovY, invXcov, VarLambdat, varrho, v_old, delta);
        // parameter differences (step 2b)
        delta = getDelta(ada_weights, beta, v_old, Lambda, varrho, N, p, parallel);
        // Lagrangian parameters (step 2c)
        resid = Lambda * beta - delta;
        v_new = v_old + varrho * resid;
        // Check for convergence (step 2d)
        iter++;
        if (stoppingCrit(resid, tol_convergence))
        {
            break;
        }
        v_old = v_new;
    }
    // Create an indicator whether convergence was achieved (also possibly on the final iteration)
    bool convergence = stoppingCrit(resid, tol_convergence);

    //------------------------------//
    // Pull the grouping            //
    //------------------------------//

    // Assign preliminary group adherences
    arma::uvec groups_hat_prelim = getGroups(beta, y_tilde, X_tilde, Lambda, p, N, tol_group);

    // Kick out trivial groups
    arma::uvec groups_hat;
    if (min_group_frac > 0.0)
    {
        groups_hat = mergeTrivialGroups(groups_hat_prelim, y_tilde, X_tilde, Z_tilde, method, min_group_frac, N, i_index, p, FALSE, parallel);
    }
    else
    {
        groups_hat = groups_hat_prelim;
    }
    // Get the total number of groups
    int K_hat = arma::max(groups_hat);

    //------------------------------//
    // Post lasso estimates         //
    //------------------------------//

    arma::mat alpha_mat = getAlpha(X_tilde, y_tilde, Z_tilde, method, N, i_index, p, groups_hat, FALSE, parallel);
    // Apply Split-panel Jackknife bias correction
    if (bias_correc)
    {
        arma::uvec i_index_tilde;
        if (method == "PGMM")
        {
            i_index_tilde = addOneObsperI(i_index);
        }
        else
        {
            i_index_tilde = i_index;
        }
        alpha_mat = spjCorrec(alpha_mat, X, y, Z, N, i_index_tilde, p, groups_hat, method, parallel);
    }

    //------------------------------//
    // Output                       //
    //------------------------------//

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("alpha_hat") = alpha_mat,
        Rcpp::Named("K_hat") = K_hat,
        Rcpp::Named("groups_hat") = groups_hat.t(),
        Rcpp::Named("iter") = iter,
        Rcpp::Named("convergence") = convergence);
    return output;
}

// Compute the BIC IC
Rcpp::List IC(const unsigned int &K, const arma::mat &alpha_hat, const arma::uvec &groups, arma::vec &y_tilde, arma::mat &X_tilde, const double &rho, const unsigned int &N, arma::uvec &i_index, const bool &fit_log)
{
    // Compute the penalty term
    int p = alpha_hat.n_cols;
    double penalty = rho * p * K;

    // Compute the fitness term
    arma::vec fit(y_tilde.n_elem);
    arma::vec resid(y_tilde.n_elem);
    arma::vec y_i, resid_i, fit_i;
    arma::rowvec alpha_vec;
    arma::uvec ind_seq;
    arma::mat X_i;
    arma::vec msr_i;
    float msr_sum = 0;
    for (unsigned int i = 0; i < N; ++i)
    {
        alpha_vec = alpha_hat.row(groups[i] - 1);
        ind_seq = arma::find(i_index == i + 1);
        y_i = y_tilde.elem(ind_seq);
        X_i = X_tilde.rows(ind_seq);
        fit_i = X_i * alpha_vec.t();
        resid_i = y_i - fit_i;
        fit(ind_seq) = fit_i;
        resid(ind_seq) = resid_i;
        msr_i = arma::sum(arma::pow(resid_i, 2));
        msr_sum += msr_i(0, 0);
    }
    float msr = msr_sum / y_tilde.n_elem;

    // Construct the IC
    double IC;
    if (fit_log) {
      IC = log(msr) + penalty;
    } else {
      IC = msr + penalty;
    }

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("IC") = IC,
        Rcpp::Named("fitted") = fit,
        Rcpp::Named("resid") = resid,
        Rcpp::Named("msr") = msr);
    return output;
}

// Combine the PAGFL algo and the IC
// [[Rcpp::export]]
Rcpp::List pagfl_routine(arma::vec &y, arma::mat &X, const std::string &method, arma::mat &Z, arma::uvec &i_index, const arma::uvec &t_index, const unsigned int &N, const bool &bias_correc, const arma::vec &lambda_vec, const double &kappa, const double &min_group_frac, const unsigned int &max_iter, const double &tol_convergence, const double &tol_group, const double &varrho, const double &rho, const bool &parallel, const bool &verbose)
{

    //------------------------------//
    // Preliminaries                //
    //------------------------------//

    unsigned int n_periods = arma::max(t_index);
    // Net out fixed effects
    std::vector<arma::mat> data = netFE(y, X, method, N, i_index);
    arma::vec y_tilde = data[0];
    arma::mat X_tilde = data[1];
    arma::mat Z_tilde;
    unsigned int p = X.n_cols;
    unsigned int n = N * (N - 1) * p / 2;
    unsigned int q = Z.n_cols;
    // Compute some constants
    arma::sp_mat Lambda = buildLambda(p, N);
    arma::sp_mat VarLambdat = varrho * Lambda.t();

    //------------------------------//
    // Initial estimates            //
    //------------------------------//

    arma::sp_mat X_block, Xt, Z_block, W;
    arma::mat XtX;
    arma::vec beta, u, Xty;
    std::tuple<arma::sp_mat, arma::sp_mat, arma::mat> X_block_vec;
    if (method == "PGMM")
    {
        Z_tilde = deleteObsMat(Z, N, i_index, TRUE);
        i_index = deleteOneObsperI(i_index);
        n_periods = n_periods - 1;
        X_block = buildDiagX_block(X_tilde, N, i_index, arma::regspace<arma::uvec>(1, N));
        Z_block = buildDiagX_block(Z_tilde, N, i_index, arma::regspace<arma::uvec>(1, N));
        // W = getW(X_block, Z_block, y_tilde, q);
        W = arma::eye<arma::mat>(Z_block.n_cols, Z_block.n_cols);
        // Recompute the initial non-penalized beta estimate
        Xt = X_block.t() * Z_block * W * Z_block.t();
        XtX = arma::mat(Xt * X_block + .05 / sqrt(y_tilde.n_elem));
        Xty = Xt * y_tilde;
        beta = arma::inv(XtX) * Xty;
    }
    else
    {
        // Build a predictor block matrix
        X_block_vec = buildDiagX(X_tilde, y_tilde, N, i_index, arma::regspace<arma::uvec>(1, N), FALSE, parallel);
        X_block = std::get<0>(X_block_vec);
        XtX = std::get<1>(X_block_vec);
        beta = std::get<2>(X_block_vec);
        Xty = X_block.t() * y_tilde;
    }
    // Pre-invert the predictor var-cov matrix
    arma::mat XtXLambda = arma::mat(XtX + VarLambdat * Lambda);
    arma::mat invXcov = arma::inv(XtXLambda);
    arma::vec invXcovY = invXcov * Xty;

    //------------------------------//
    // Initialize the algorithm     //
    //------------------------------//

    // Parameter differences
    arma::vec delta = Lambda * beta;
    // Lagrangian parameters
    arma::vec v_old(n, fill::zeros);
    // Compute the adaptive weights
    arma::vec omega = getOmega(delta, kappa, N, p);

    //------------------------------//
    // Run the algorithm            //
    //------------------------------//

    Rcpp::List estimOutput, IC_list, output;
    Rcpp::List lambdalist(lambda_vec.n_elem);
    for (unsigned int l = 0; l < lambda_vec.n_elem; l++)
    {
      // Clear the line
      if (verbose)  Rcout << "\r" << std::string(80, ' ') << "\r";
        // Estimate
        estimOutput = pagfl_algo(y, y_tilde, X, X_tilde, invXcovY, invXcov, VarLambdat, Lambda, method, Z, Z_tilde, delta, omega, v_old, i_index, t_index, N, n, p, q, n_periods, bias_correc, lambda_vec[l], min_group_frac, max_iter, tol_convergence, tol_group, varrho, parallel, verbose, l + 1, lambda_vec.n_elem);
        // Compute the Information criterion
        IC_list = IC(Rcpp::as<unsigned int>(estimOutput["K_hat"]), Rcpp::as<arma::mat>(estimOutput["alpha_hat"]), Rcpp::as<arma::uvec>(estimOutput["groups_hat"]), y_tilde, X_tilde, rho, N, i_index, FALSE);
        output = Rcpp::List::create(
            Rcpp::Named("estimOutput") = estimOutput,
            Rcpp::Named("IC") = IC_list);
        lambdalist[l] = output;
    }
    return lambdalist;
}

// Time-varying PAGFL routine
Rcpp::List tv_pagfl_algo(arma::vec &y_tilde, arma::mat &Z_tilde, arma::vec &invZcovY, arma::mat &invZcov, const arma::vec &delta_ini, const arma::vec &omega, arma::vec &v_old, const arma::sp_mat &VarLambdat, const arma::sp_mat &Lambda, const arma::mat &B, const unsigned int d, arma::uvec &i_index, unsigned int n_periods, const unsigned int N, const unsigned int n, const unsigned int p_star, const double lambda, const double &min_group_frac, const unsigned int &max_iter, const double &tol_convergence, const double &tol_group, const double &varrho, const bool &parallel, const bool &verbose, const unsigned int &lambda_num, const unsigned int &n_lambda)
{

    //------------------------------//
    // Initialize                   //
    //------------------------------//

    arma::vec delta = delta_ini;
    float lambda_star = n_periods * lambda / (2 * N);
    arma::vec ada_weights = omega * lambda_star / varrho;

    //------------------------------//
    // Run the algorithm            //
    //------------------------------//

    arma::vec resid, v_new, pi;
    unsigned int iter = 0;
    for (unsigned int i = 0; i < max_iter; i++)
    {
        if (verbose)
        {
            Rcout << "\r" << "Lambda: " << lambda << " (" << lambda_num << "/" << n_lambda << ")" << " - Iteration: " << i + 1 << "/" << max_iter << std::flush;
        }
        // Update the ...
        // parameter estimates (step 2a)
        pi = getBeta(invZcovY, invZcov, VarLambdat, varrho, v_old, delta);
        // parameter differences (step 2b)
        delta = getDelta(ada_weights, pi, v_old, Lambda, varrho, N, p_star, parallel);
        // Lagrangian parameters (step 2c)
        resid = Lambda * pi - delta;
        v_new = v_old + varrho * resid;
        // Check for convergence (step 2d)
        iter++;
        if (stoppingCrit(resid, tol_convergence))
        {
            break;
        }
        v_old = v_new;
    }

    // Create an indicator whether convergence was achieved (also possibly on the final iteration)
    bool convergence = stoppingCrit(resid, tol_convergence);

    //------------------------------//
    // Pull the grouping            //
    //------------------------------//

    // Assign preliminary group adherences
    arma::uvec groups_hat_prelim = getGroups(pi, y_tilde, Z_tilde, Lambda, p_star, N, tol_group);
    // Kick out trivial groups
    arma::uvec groups_hat;
    arma::mat R;
    if (min_group_frac > 0.0)
    {
        groups_hat = mergeTrivialGroups(groups_hat_prelim, y_tilde, Z_tilde, R, "PLS", min_group_frac, N, i_index, p_star, TRUE, parallel);
    }
    else
    {
        groups_hat = groups_hat_prelim;
    }
    // Get the total number of groups
    int K_hat = arma::max(groups_hat);

    //------------------------------//
    // Post lasso estimates         //
    //------------------------------//

    arma::mat xi_mat_vec = getAlpha(Z_tilde, y_tilde, R, "PLS", N, i_index, p_star, groups_hat, TRUE, parallel);

    //------------------------------//
    // Output                       //
    //------------------------------//

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("alpha_hat") = xi_mat_vec,
        Rcpp::Named("K_hat") = K_hat,
        Rcpp::Named("groups_hat") = groups_hat.t(),
        Rcpp::Named("iter") = iter,
        Rcpp::Named("convergence") = convergence);
    return output;
}

// Combine the time-varying PAGFL algo and the IC
// [[Rcpp::export]]
Rcpp::List tv_pagfl_routine(arma::vec &y, arma::mat &X, arma::mat &X_const, const unsigned int &d, const unsigned int &M, arma::uvec &i_index, const arma::uvec &t_index, const unsigned int &N, const unsigned int &p_const, const arma::vec &lambda_vec, const double &kappa, const double &min_group_frac, const unsigned int &max_iter, const double &tol_convergence, const double &tol_group, const double &varrho, const double &rho, const bool &parallel, const bool &verbose)
{

    //------------------------------//
    // Build the B-spline basis     //
    //------------------------------//

    unsigned int n_periods = arma::max(t_index);
    arma::vec knots = arma::linspace(1, n_periods, M + 2);
    arma::vec support = arma::regspace<arma::vec>(1, n_periods);
    arma::mat B = bspline_system(support, d, knots, TRUE);
    arma::mat Z;
    Z = buildZ(X, B, t_index, X.n_cols);
    if (p_const > 0)
    {
        Z = join_rows(Z, X_const);
    }

    //------------------------------//
    // Preliminaries                //
    //------------------------------//

    // Net out fixed effects
    std::vector<arma::mat> data = netFE(y, Z, "PLS", N, i_index);
    arma::vec y_tilde = data[0];
    arma::mat Z_tilde = data[1];

    // Compute some constants
    unsigned int p_star = Z.n_cols;
    unsigned int n = N * (N - 1) * p_star / 2;
    arma::sp_mat Lambda = buildLambda(p_star, N);
    arma::sp_mat VarLambdat = varrho * Lambda.t();

    //------------------------------//
    // Initial estimates            //
    //------------------------------//

    // Build a predictor block matrix (either block-wise or once, depending on the dimensionality)
    std::tuple<arma::sp_mat, arma::sp_mat, arma::mat> Z_block_vec = buildDiagX(Z_tilde, y_tilde, N, i_index, arma::regspace<arma::uvec>(1, N), TRUE, parallel);
    arma::sp_mat Z_block = std::get<0>(Z_block_vec);
    arma::sp_mat ZtZ = std::get<1>(Z_block_vec);
    arma::vec pi = std::get<2>(Z_block_vec);
    arma::mat Zty = Z_block.t() * y_tilde;
    // Pre-invert the predictor var-cov matrix
    arma::mat ZtZLambda = arma::mat(ZtZ + VarLambdat * Lambda);
    arma::mat invZcov = arma::pinv(ZtZLambda);
    arma::vec invZcovY = invZcov * Zty;

    //------------------------------//
    // Initialize the algorithm     //
    //------------------------------//

    // Parameter differences
    arma::vec delta = Lambda * pi;
    // Lagrangian parameters
    arma::vec v_old(n, fill::zeros);
    // Compute the adaptive weights
    arma::vec omega = getOmega(delta, kappa, N, p_star);

    //------------------------------//
    // Run the algorithm            //
    //------------------------------//

    Rcpp::List estimOutput, IC_list, output;
    Rcpp::List lambdalist(lambda_vec.n_elem);
    for (unsigned int l = 0; l < lambda_vec.n_elem; l++)
    {
      // Clear the line
      if (verbose)  Rcout << "\r" << std::string(80, ' ') << "\r";
        // Estimate
        estimOutput = tv_pagfl_algo(y_tilde, Z_tilde, invZcovY, invZcov, delta, omega, v_old, VarLambdat, Lambda, B, d, i_index, n_periods, N, n, p_star, lambda_vec[l], min_group_frac, max_iter, tol_convergence, tol_group, varrho, parallel, verbose, l + 1, lambda_vec.n_elem);
        // Compute the Information Criterion
        IC_list = IC(Rcpp::as<unsigned int>(estimOutput["K_hat"]), Rcpp::as<arma::mat>(estimOutput["alpha_hat"]), Rcpp::as<arma::uvec>(estimOutput["groups_hat"]), y_tilde, Z_tilde, rho, N, i_index, TRUE);
        output = Rcpp::List::create(
            Rcpp::Named("estimOutput") = estimOutput,
            Rcpp::Named("IC") = IC_list);
        lambdalist[l] = output;
    }

    return lambdalist;
}

// Compute time varying coefficients from spline bases
// [[Rcpp::export]]
arma::cube getTVAlpha(const arma::mat &xi, const unsigned int &K_hat, const unsigned int &p, const unsigned int &n_periods, const arma::mat &B)
{
    arma::cube alpha_array(n_periods, p, K_hat, arma::fill::zeros);
    unsigned int J_star = xi.n_cols / p;
    for (unsigned int k = 0; k < K_hat; k++)
    {
        arma::mat xi_mat = xi.row(k);
        xi_mat.reshape(J_star, p);
        alpha_array.slice(k) = B * xi_mat;
    }
    return alpha_array;
}

// [[Rcpp::export]]
arma::cube delete_missing_t(const arma::uvec &i_index, const arma::uvec &t_index, const unsigned int &K_hat, const arma::vec &groups_hat, arma::cube &alpha_hat)
{
    arma::umat df = arma::join_rows(i_index, t_index);
    int n_periods = alpha_hat.n_rows;
    double min_t, max_t;
    arma::uvec g_k, current_indices;
    arma::umat tmp;
    arma::uvec row_indices;
    for (unsigned int k = 0; k < K_hat; k++)
    {
        g_k = arma::find(groups_hat == k + 1);
        arma::uvec row_indices_temp;
        for (unsigned int j = 0; j < g_k.n_elem; j++)
        {
            current_indices = arma::find(df.col(0) == g_k[j] + 1);
            row_indices_temp = arma::join_vert(row_indices_temp, current_indices);
        }
        row_indices = row_indices_temp;
        tmp = df.rows(row_indices);
        min_t = arma::min(tmp.col(1));
        max_t = arma::max(tmp.col(1));
        if (min_t > 1)
            alpha_hat.subcube(0, 0, k, min_t - 2, alpha_hat.n_cols - 1, k).fill(arma::datum::nan);
        if (max_t < n_periods)
            alpha_hat.subcube(max_t, 0, k, n_periods - 1, alpha_hat.n_cols - 1, k).fill(arma::datum::nan);
    }
    return alpha_hat;
}

// [[Rcpp::export]]
arma::vec fitMeasures(unsigned int &N, const unsigned int &k, arma::vec &y, arma::uvec &i_index, const std::string &method, const double &msr)
{
    arma::vec y_tilde;
    if (method == "PLS")
    {
        y_tilde = demeanIndVec(y, N, i_index);
    }
    else
    {
        y_tilde = fdIndVec(y, N, i_index);
    }
    arma::vec ssq = sum(arma::pow(y_tilde, 2));
    unsigned int n = y_tilde.n_elem;
    double r_df = n - N - k;
    double ssr = msr * n;
    float r_se = sqrt(ssr / r_df);
    float r_sq = 1 - ssr / ssq(0, 0);
    float adj_r_sq = 1 - (1 - r_sq) * (n - 1) / r_df;
    arma::vec out = {r_df, r_sq, adj_r_sq, r_se};
    return out;
}

// [[Rcpp::export]]
arma::vec getFE(const arma::vec &y, const arma::uvec &i_index, const unsigned int &N, const std::string &method)
{
    arma::uvec ind_vec;
    arma::vec y_i;
    arma::uvec i_index_tilde = i_index;
    if (method == "PGMM")
    {
        i_index_tilde = deleteOneObsperI(i_index_tilde);
    }
    arma::vec fe_vec(i_index_tilde.n_elem);
    double fe;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_vec = arma::find(i_index_tilde == i + 1);
        y_i = y.elem(ind_vec);
        fe = arma::mean(y_i);
        fe_vec(ind_vec).fill(fe);
    }
    return fe_vec;
}

// [[Rcpp::export]]
Rcpp::List tv_pagfl_oracle_routine(arma::vec &y, arma::mat &X, arma::mat &X_const, const unsigned int &d, const arma::uvec &groups, const unsigned int &M, arma::uvec &i_index, const arma::uvec &t_index, const unsigned int &N, const unsigned int &p_const, const double &rho, const bool &parallel)
{

    //------------------------------//
    // Build the B-spline basis     //
    //------------------------------//

    unsigned int n_periods = arma::max(t_index);
    arma::vec knots = arma::linspace(1, n_periods, M + 2);
    arma::vec support = arma::regspace<arma::vec>(1, n_periods);
    arma::mat B = bspline_system(support, d, knots, TRUE);
    arma::mat Z;
    Z = buildZ(X, B, t_index, X.n_cols);
    if (p_const > 0)
    {
        Z = join_rows(Z, X_const);
    }

    //------------------------------//
    // Preliminaries                //
    //------------------------------//

    unsigned int n_groups = arma::max(groups);
    // Net out fixed effects
    std::vector<arma::mat> data = netFE(y, Z, "PLS", N, i_index);
    arma::vec y_tilde = data[0];
    arma::mat Z_tilde = data[1];
    // Compute some constants
    unsigned int p_star = Z.n_cols;

    //------------------------------//
    // Run the estimation           //
    //------------------------------//

    arma::mat xi_mat_vec = getGroupwiseOLS(y_tilde, Z_tilde, N, i_index, groups, p_star, TRUE, parallel);

    bool convergence = TRUE;
    Rcpp::List estimOutput = Rcpp::List::create(
        Rcpp::Named("alpha_hat") = xi_mat_vec,
        Rcpp::Named("K_hat") = n_groups,
        Rcpp::Named("groups_hat") = groups.t(),
        Rcpp::Named("iter") = 0,
        Rcpp::Named("convergence") = convergence);

    //------------------------------//
    // Compute the IC               //
    //------------------------------//

    Rcpp::List IC_list = IC(n_groups, xi_mat_vec, groups, y_tilde, Z_tilde, rho, N, i_index, TRUE);

    //------------------------------//
    // Output                       //
    //------------------------------//

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("estimOutput") = estimOutput,
        Rcpp::Named("IC") = IC_list);
    return output;
}

// [[Rcpp::export]]
Rcpp::List pagfl_oracle_routine(arma::vec &y, arma::mat &X, const arma::uvec &groups, const std::string &method, arma::mat &Z, arma::uvec i_index, const arma::uvec &t_index, const unsigned int &N, const bool &bias_correc, const double &rho, const bool &parallel)
{

    //------------------------------//
    // Preliminaries                //
    //------------------------------//

    unsigned int n_groups = arma::max(groups);
    unsigned int n_periods = arma::max(t_index);
    // Net out fixed effects
    std::vector<arma::mat> data = netFE(y, X, method, N, i_index);
    arma::vec y_tilde = data[0];
    arma::mat X_tilde = data[1];
    arma::mat Z_tilde;
    unsigned int p = X.n_cols;

    if (method == "PGMM")
    {
        Z_tilde = deleteObsMat(Z, N, i_index, TRUE);
        i_index = deleteOneObsperI(i_index);
        n_periods = n_periods - 1;
    }

    //------------------------------//
    // Run the estimation           //
    //------------------------------//

    arma::mat alpha_mat = getAlpha(X_tilde, y_tilde, Z_tilde, method, N, i_index, p, groups, FALSE, parallel);
    // Apply Split-panel Jackknife bias correction
    if (bias_correc)
    {
        arma::uvec i_index_tilde;
        if (method == "PGMM")
        {
            i_index_tilde = addOneObsperI(i_index);
        }
        else
        {
            i_index_tilde = i_index;
        }
        alpha_mat = spjCorrec(alpha_mat, X, y, Z, N, i_index_tilde, p, groups, method, parallel);
    }

    //------------------------------//
    // Compute the IC               //
    //------------------------------//

    Rcpp::List IC_list = IC(n_groups, alpha_mat, groups, y_tilde, X_tilde, rho, N, i_index, FALSE);

    //------------------------------//
    // Output                       //
    //------------------------------//

    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("alpha_hat") = alpha_mat,
        Rcpp::Named("IC") = IC_list);
    return output;
}
