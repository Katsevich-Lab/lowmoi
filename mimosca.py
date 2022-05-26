import numpy as np
from scipy import sparse
from sklearn.linear_model import ElasticNet
from statsmodels.distributions.empirical_distribution import ECDF

##############
# MY FUNCTIONS
##############
def get_dense_array(pieces):
    data = pieces[0]
    indices = pieces[1]
    indptr = pieces[2]
    n_row = pieces[3]
    n_col = pieces[4]
    x = sparse.csc_array((data, indices, indptr), (n_row, n_col))
    ret = x.toarray()
    return(ret)


def run_mimosca(gene_mat_pieces, gRNA_mat_pieces, cov_ind, num_iters=1000):
    # convert the gene and gRNA matrices into dense numpy matrices
    gene_mat_dense = get_dense_array(gene_mat_pieces)
    gRNA_mat_dense = get_dense_array(gRNA_mat_pieces)
    # estimate the coefficients
    B_mat = fit_lm(gRNA_mat_dense, gene_mat_dense)
    # compute the p-values for gRNA at cov_ind
    null_distribution = shuffle_and_fit(gRNA_mat_dense, gene_mat_dense, cov_ind, num_iters)
    p_vals = calc_p_vals(B_mat, null_distribution, cov_ind)
    return(p_vals)


###################
# CHRIS'S FUNCTIONS
###################
def fit_lm(X, y, l1_ratio=0.5, alpha=0.0005, max_iter=10000, z_score=False):

	lmfit = ElasticNet(precompute=True, l1_ratio=l1_ratio, alpha=alpha, max_iter=max_iter)
	lmfit.fit(X, y)

	return lmfit.coef_


def shuffle_and_fit(X, y, cov_ind, num_iters):
	# Can pre-allocate "all_nulls" for speed
	for iter in np.arange(num_iters):
		X_shuffled = X.copy()
		X_shuffled[:, cov_ind] = np.random.permutation(X_shuffled[:, cov_ind])
		lm_coefs = fit_lm(X_shuffled, y)

		if iter == 0:
			all_nulls = lm_coefs[:, cov_ind].flatten()
		else:
			all_nulls = np.append(all_nulls, lm_coefs[:, cov_ind].flatten())
		del X_shuffled

	return all_nulls


def calc_p_vals(beta_mat, null_distrib, cov_ind):

	curr_coeffs = beta_mat[:, cov_ind].flatten()
	curr_ECDF = ECDF(null_distrib)
	p_vals = np.ones(curr_coeffs.size)

	neg_inds = np.where(curr_coeffs < 0)[0]
	pos_inds = np.where(curr_coeffs >= 0)[0]

	p_vals[neg_inds] = curr_ECDF(curr_coeffs[neg_inds])
	p_vals[pos_inds] = 1-curr_ECDF(curr_coeffs[pos_inds])

	return p_vals
