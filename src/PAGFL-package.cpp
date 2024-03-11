#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// Delete one observation per individual unit to the index vector
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
arma::uvec addOneObsperI(const arma::uvec &vec) {
    arma::uvec unique_values = arma::unique(vec);
    arma::uvec result = vec;
    arma::uvec indices;
    for (unsigned int i = 0; i < unique_values.n_elem; ++i) {
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
arma::vec demeanIndVec(arma::vec x, unsigned int N, arma::uvec i_index)
{
    arma::vec x_tilde = x;
    arma::uvec ind_seq;
    arma::vec x_red;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        x_red = x.elem(ind_seq);
        x_tilde.elem(ind_seq) = x_red - arma::mean(x_red);
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
arma::mat buildLambda(const unsigned int &p, const unsigned int &N)
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
    arma::mat Lambda = arma::kron(D, arma::eye<arma::mat>(p, p));
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
arma::vec getBeta(arma::mat invXcov, arma::mat XtY, double varrho, arma::mat VarLambdat, arma::vec v, arma::vec delta)
{
    // See Mehrabani (2023, sec 5.1 step 2a/ 5.2 step 2a)
    arma::mat XYcov = XtY + VarLambdat * (delta - v / varrho);
    arma::vec beta = invXcov * XYcov;
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

// Update the parameter difference vector (see Mehrabani (2023, sec 5.1 step 2b))
arma::vec getDelta(const arma::vec &ada_weights, const arma::vec &beta, const arma::vec &v, const arma::mat &Lambda, const double &varrho, const unsigned int &N, const unsigned int &p)
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
    for (unsigned int i = 0; i < n / p; i++)
    {
        arma::uvec ind = arma::conv_to<arma::uvec>::from(ind_mat.row(i)) - 1;
        delta.subvec(i * p, (i + 1) * p - 1) = softThreshold(ind, xi, ada_weights);
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
arma::mat buildBlockDiag(const std::vector<arma::mat> &XList)
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
    return block_diag;
}

// Constructs a NT x Kp block predictor matrix
// [[Rcpp::export]]
arma::mat buildDiagX(const arma::mat &X, const unsigned int &N, arma::uvec &i_index, const arma::uvec &groups)
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
    arma::mat X_tilde = buildBlockDiag(XMat_tilde_vec);
    return X_tilde;
}

// inverts the error var-cov matrix for each group
arma::mat invertV(const arma::mat &V, const unsigned int &q)
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
    arma::mat W = buildBlockDiag(W_list);
    return W;
}

// Compute the weight matrix for the GMM
arma::mat getW(const arma::mat &X_block, arma::mat &Z_block, const arma::vec &y, const unsigned int &q)
{
    // Initial estimate with an identity weight matrix
    arma::mat Xt = X_block.t() * Z_block * Z_block.t();
    arma::mat XtX = Xt * X_block;
    arma::mat Xty = Xt * y;
    arma::vec beta = arma::inv(XtX + .05 / sqrt(y.n_elem)) * Xty;
    // Obtain the (group-wise) var-cov matrix of the residuals
    arma::vec u = y - X_block * beta;
    arma::mat V = Z_block.t() * u * u.t() * Z_block;
    // Use the inverse as the subsequent weight matrix (tossing all matrices off the diagonal)
    arma::mat W = invertV(V, q);
    return W;
}

// Computes the post-lasso estimates based on an estimated grouping
// [[Rcpp::export]]
arma::mat getAlpha(const arma::mat &X, const arma::vec &y, arma::mat &Z, const std::string &method, const unsigned int &N, arma::uvec &i_index, const unsigned int &p, const arma::uvec &groups_hat)
{
    arma::mat groupX = buildDiagX(X, N, i_index, groups_hat);
    arma::mat groupXt = groupX.t();

    // Order the individuals in the dependent variable vector according to the groups
    std::vector<arma::vec> y_vec(N);
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        y_vec[i] = y.elem(ind_seq);
    }
    // Create a vector of pairs where the first element of each pair is the group and the second element is the original index
    std::vector<std::pair<int, int>> group_index_pairs(N);
    for (unsigned int i = 0; i < N; ++i)
    {
        group_index_pairs[i] = std::make_pair(groups_hat[i], i);
    }
    // Sort the vector of pairs based on the group
    std::sort(group_index_pairs.begin(), group_index_pairs.end());
    // Create a new y vector where the entries are sorted based on the group structure
    arma::vec y_tilde;
    for (unsigned int i = 0; i < group_index_pairs.size(); ++i)
    {
        y_tilde = arma::join_cols(y_tilde, y_vec[group_index_pairs[i].second]);
    }

    // Compute the post-lasso estimates
    arma::vec alpha;
    arma::mat groupZ, W;
    if (method == "PGMM")
    {
        groupZ = buildDiagX(Z, N, i_index, groups_hat);
        // W = getW(groupX, groupZ, y_tilde, Z.n_cols);
        W = arma::eye<arma::mat>(groupZ.n_cols, groupZ.n_cols);
        groupXt = groupXt * groupZ * W * groupZ.t();
        // Compute the post-lasso estimates
        alpha = arma::inv(groupXt * groupX + .05 / sqrt(y.n_elem)) * groupXt * y_tilde;
    }
    else
    {
        alpha = arma::pinv(groupXt * groupX) * groupXt * y_tilde;
    }
    arma::mat alpha_mat = arma::reshape(alpha, p, alpha.n_elem / p).t();
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
// [[Rcpp::export]]
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

// Drops the last observed time period of a NT vector per individual
arma::vec deleteLastObsVec(const arma::vec &X, const unsigned int &N, arma::uvec i_index)
{
    arma::vec x_red, X_red;
    arma::uvec ind_seq;
    for (unsigned int i = 0; i < N; i++)
    {
        ind_seq = arma::find(i_index == i + 1);
        ind_seq.shed_row(ind_seq.n_elem - 1);
        x_red = X.elem(ind_seq);
        X_red = arma::join_cols(X_red, x_red);
    }
    return X_red;
}

// Net out fixed effects
// [[Rcpp::export]]
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

    arma::mat X_tilde = netFEMat(X, N, i_index);
    arma::mat y_tilde = netFEVec(y, N, i_index);
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
// [[Rcpp::export]]
arma::mat spjCorrec(const arma::mat &alpha_mat, arma::mat &X, arma::vec &y, arma::mat &Z, const unsigned int &N, arma::uvec &i_index, const unsigned int &p, const arma::uvec &groups_hat, const std::string &method)
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
        alpha_half = getAlpha(X_half_tilde, y_half_tilde, Z_half_tilde, method, N, i_index_half, p, groups_hat);
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
arma::uvec getGroups(const arma::vec &beta, const arma::vec &y, const arma::mat &X, const arma::mat &Lambda, const unsigned int &p, const unsigned int &N, const double &tol)
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
    arma::mat X_block = buildDiagX(X, N, i_index, individuals);
    // Initialize the alpha matrix
    arma::mat alpha_mat(groups.n_elem, alpha.n_cols);
    // Fill with one group per row
    for (unsigned int i = 0; i < groups.n_elem; ++i)
    {
        alpha_mat.row(i) = alpha.row(groups[i] - 1);
    }
    // Compute the ssq
    arma::vec alpha_vec = arma::vectorise(alpha_mat.t());
    float ssq = arma::mean(arma::pow(y - X_block * alpha_vec, 2));
    return ssq;
}

// Allocate individuals that are part of trivial groups to the best-fitting non-trivial group
arma::uvec mergeTrivialGroups(arma::uvec &groups_hat, const arma::vec &y, const arma::mat &X, arma::mat &Z, const std::string &method, const float &min_group_frac, const unsigned int &N, arma::uvec &i_index, const unsigned int &p)
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
                arma::mat alpha_tmp = getAlpha(X, y, Z, method, N, i_index, p, groups_hat_tmp_label);
                ssq_vec(j) = getSSQ(y, X, alpha_tmp, groups_hat_tmp_label, N, i_index);
            }
            // Select the non-trivial group with the best fit
            groups_hat(trivialInd(i)) = nonTrivialGroups(arma::index_min(ssq_vec));
        }
        return relabelGroups(groups_hat);
    }
}

// PAGFL routine
// [[Rcpp::export]]
Rcpp::List pagfl_algo(arma::vec &y, arma::mat &X, const std::string &method, arma::mat &Z, arma::uvec &i_index, const arma::uvec &t_index, const unsigned int &N, const bool &bias_correc, const double &lambda, const double &kappa, const double &min_group_frac, const unsigned int &max_iter, const double &tol_convergence, const double &tol_group, const double &varrho)
{

    //------------------------------//
    // Preliminaries                //
    //------------------------------//

    unsigned int n_periods = arma::max(t_index);
    unsigned int p = X.n_cols;
    unsigned int n = N * (N - 1) * p / 2;
    unsigned int q = Z.n_cols;
    // Compute some constants
    arma::mat Lambda = buildLambda(p, N);
    arma::mat VarLambdat = varrho * Lambda.t();
    float lambda_star;
    if (method == "PLS")
    {
        lambda_star = n_periods * lambda / (2 * N);
    }
    else
    {
        lambda_star = std::pow(n_periods, 2) * lambda / (2 * N);
    }

    // Net out fixed effects
    std::vector<arma::mat> data = netFE(y, X, method, N, i_index);
    arma::vec y_tilde = data[0];
    arma::mat X_tilde = data[1];
    arma::mat Z_tilde;
    if (method == "PGMM")
    {
        Z_tilde = deleteObsMat(Z, N, i_index, TRUE);
        i_index = deleteOneObsperI(i_index);
        n_periods = n_periods - 1;
    }
    arma::mat X_block = buildDiagX(X_tilde, N, i_index, arma::regspace<arma::uvec>(1, N));

    //------------------------------//
    // Initial estimates            //
    //------------------------------//

    arma::mat Xt = X_block.t();
    arma::mat XtX = Xt * X_block;
    arma::vec u, beta;
    arma::mat W, Z_block;
    if (method == "PGMM")
    {
        Z_block = buildDiagX(Z_tilde, N, i_index, arma::regspace<arma::uvec>(1, N));
        W = getW(X_block, Z_block, y_tilde, q);
        W = arma::eye<arma::mat>(Z_block.n_cols, Z_block.n_cols);
        // Recompute the initial non-penalized beta estimate
        Xt = Xt * Z_block * W * Z_block.t();
        XtX = Xt * X_block + .05 / sqrt(y_tilde.n_elem);
    }
    arma::mat Xty = Xt * y_tilde;
    beta = arma::inv(XtX) * Xty;
    // Pre-invert the predictor var-cov matrix
    arma::mat invXcov = arma::inv(XtX + VarLambdat * Lambda);
    // Parameter differences
    arma::vec delta = Lambda * beta;
    // Lagrangian parameters
    arma::vec v_old(n, fill::zeros);
    // Compute the adaptive weights
    arma::vec omega = getOmega(delta, kappa, N, p);
    arma::vec ada_weights = omega * lambda_star / varrho;

    //------------------------------//
    // Run the algorithm            //
    //------------------------------//

    arma::vec resid, v_new;
    unsigned int iter = 0;
    for (unsigned int i = 0; i < max_iter; i++)
    {
        // Update the ...
        // parameter estimates (sec. 5.1/ 5.2 step 2a)
        beta = getBeta(invXcov, Xty, varrho, VarLambdat, v_old, delta);
        // parameter differences (step 2b)
        delta = getDelta(ada_weights, beta, v_old, Lambda, varrho, N, p);
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
    // Assign preliminary group adherences
    arma::uvec groups_hat_prelim = getGroups(beta, y_tilde, X_tilde, Lambda, p, N, tol_group);
    // Kick out trivial groups
    arma::uvec groups_hat;
    if (min_group_frac > 0.0)
    {
        groups_hat = mergeTrivialGroups(groups_hat_prelim, y_tilde, X_tilde, Z_tilde, method, min_group_frac, N, i_index, p);
    }
    else
    {
        groups_hat = groups_hat_prelim;
    }
    // Get the total number of groups
    int K_hat = arma::max(groups_hat);
    // Post lasso estimates
    arma::mat alpha_mat = getAlpha(X_tilde, y_tilde, Z_tilde, method, N, i_index, p, groups_hat);
    // Apply Split-panel Jackknife bias correction
    if (bias_correc)
    {
        if (method == "PGMM")
        {
            n_periods = n_periods + 1;
            i_index = addOneObsperI(i_index);
        }
        alpha_mat = spjCorrec(alpha_mat, X, y, Z, N, i_index, p, groups_hat, method);
    }
    // Return the estimates
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("alpha_hat") = alpha_mat,
        Rcpp::Named("K_hat") = K_hat,
        Rcpp::Named("groups_hat") = groups_hat.t(),
        Rcpp::Named("iter") = iter,
        Rcpp::Named("convergence") = convergence);
    return output;
}

// Compute the BIC IC
// [[Rcpp::export]]
double IC(const Rcpp::List &estimOutput, arma::vec &y, arma::mat &X, const double &rho, const std::string &method, const unsigned int &N, arma::uvec &i_index)
{
    // Net out fixed effects
    std::vector<arma::mat> data = netFE(y, X, method, N, i_index);
    arma::vec y_tilde = data[0];
    arma::mat X_tilde = data[1];
    if (method == "PGMM")
    {
        i_index = deleteOneObsperI(i_index);
    }
    arma::mat X_block = buildDiagX(X_tilde, N, i_index, arma::regspace<arma::uvec>(1, N));
    // Compute the penalty term
    int K = Rcpp::as<int>(estimOutput["K_hat"]);
    arma::mat alpha_hat = Rcpp::as<arma::mat>(estimOutput["alpha_hat"]);
    int p = alpha_hat.n_cols;
    double penalty = rho * p * K;
    // Pull the group composition
    arma::uvec groups = Rcpp::as<arma::uvec>(estimOutput["groups_hat"]);
    // Initialize the alpha matrix
    arma::mat alpha_mat(groups.n_elem, alpha_hat.n_cols);
    // Fill with one group per row
    for (unsigned int i = 0; i < groups.n_elem; ++i)
    {
        alpha_mat.row(i) = alpha_hat.row(groups[i] - 1);
    }
    // Compute the fitness term
    arma::vec alpha = arma::vectorise(alpha_mat.t());
    double fitness = arma::mean(arma::pow(y_tilde - X_block * alpha, 2));
    // Construct the IC
    double IC = fitness + penalty;
    return IC;
}

// Time-varying PAGFL routine
// [[Rcpp::export]]
Rcpp::List dyn_pagfl_algo(arma::vec &y, arma::mat &Z, const arma::mat &B, const unsigned int &d, const unsigned int &J, arma::uvec &i_index, const arma::uvec &t_index, const unsigned int N, const double &lambda, const double &kappa, const double &min_group_frac, const unsigned int &max_iter, const double &tol_convergence, const double &tol_group, const double &varrho)
{

    //------------------------------//
    // Preliminaries                //
    //------------------------------//

    unsigned int p_star = Z.n_cols;
    unsigned int n = N * (N - 1) * p_star / 2;
    unsigned int n_periods = arma::max(t_index);
    // Compute some constants
    arma::mat Lambda = buildLambda(p_star, N);
    arma::mat VarLambdat = varrho * Lambda.t();
    float lambda_star = n_periods * lambda / (2 * N);

    // Net out fixed effects
    std::vector<arma::mat> data = netFE(y, Z, "PLS", N, i_index);
    arma::vec y_tilde = data[0];
    arma::mat Z_tilde = data[1];
    arma::mat Z_block = buildDiagX(Z_tilde, N, i_index, arma::regspace<arma::uvec>(1, N));

    //------------------------------//
    // Initial estimates            //
    //------------------------------//

    arma::mat Zt = Z_block.t();
    arma::mat ZtZ = Zt * Z_block;
    arma::mat Zty = Zt * y_tilde;
    arma::vec pi = arma::pinv(ZtZ) * Zty;
    // Pre-invert the predictor var-cov matrix
    arma::mat invZcov = arma::pinv(ZtZ + VarLambdat * Lambda);
    // Parameter differences
    arma::vec delta = Lambda * pi;
    // Lagrangian parameters
    arma::vec v_old(n, fill::zeros);
    // Compute the adaptive weights
    arma::vec omega = getOmega(delta, kappa, N, p_star);
    arma::vec ada_weights = omega * lambda_star / varrho;

    //------------------------------//
    // Run the algorithm            //
    //------------------------------//

    arma::vec resid, v_new;
    unsigned int iter = 0;
    for (unsigned int i = 0; i < max_iter; i++)
    {
        // Update the ...
        // parameter estimates (step 2a)
        pi = getBeta(invZcov, Zty, varrho, VarLambdat, v_old, delta);
        // parameter differences (step 2b)
        delta = getDelta(ada_weights, pi, v_old, Lambda, varrho, N, p_star);
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
    // Assign preliminary group adherences
    arma::uvec groups_hat_prelim = getGroups(pi, y_tilde, Z_tilde, Lambda, p_star, N, tol_group);
    // Kick out trivial groups
    arma::uvec groups_hat;
    arma::mat R;
    if (min_group_frac > 0.0)
    {
        groups_hat = mergeTrivialGroups(groups_hat_prelim, y_tilde, Z_tilde, R, "PLS", min_group_frac, N, i_index, p_star);
    }
    else
    {
        groups_hat = groups_hat_prelim;
    }
    // Get the total number of groups
    int K_hat = arma::max(groups_hat);
    // Post lasso estimates
    arma::mat xi_mat_vec = getAlpha(Z_tilde, y_tilde, R, "PLS", N, i_index, p_star, groups_hat);
    // Return the estimates
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("alpha_hat") = xi_mat_vec,
        Rcpp::Named("K_hat") = K_hat,
        Rcpp::Named("groups_hat") = groups_hat.t(),
        Rcpp::Named("iter") = iter,
        Rcpp::Named("convergence") = convergence);
    return output;
}

// Compute time varying coefficients from spline bases
// [[Rcpp::export]]
arma::cube getDynAlpha(const arma::mat &xi, const unsigned int &K_hat, const unsigned int &p, const unsigned int &n_periods, const arma::mat &B)
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

// Project the regressors on the spline bases
// [[Rcpp::export]]
arma::mat buildZ(const arma::mat &X, const arma::mat &B, const arma::uvec &t_index, const unsigned int &J, const unsigned int &d, const unsigned int &p)
{
    arma::mat Z(X.n_rows, (J + d) * p, arma::fill::zeros);
    unsigned int t;
    for (unsigned int i = 0; i < t_index.n_elem; ++i)
    {
        t = t_index[i] - 1;
        Z.row(i) = arma::kron(X.row(i), B.row(t));
    }
    return Z;
}
