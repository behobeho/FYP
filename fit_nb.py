import numpy as np
#from fast_negative_binomial import optimise_all_genes, optimise_all_genes_zi


def _init_nb_moments(counts_dense, eps=1e-12, r_max=1e6):
    X = np.asarray(counts_dense, dtype=np.float64)
    mu = X.mean(axis=0)
    var = X.var(axis=0, ddof=1) if X.shape[0] > 1 else np.zeros_like(mu)
    mu = np.maximum(mu, eps)
    denom = var - mu

    r = np.empty_like(mu)
    poisson_like = denom <= 0
    r[poisson_like] = float(r_max)
    r[~poisson_like] = (mu[~poisson_like] * mu[~poisson_like]) / denom[~poisson_like]
    r = np.clip(r, 1e-8, float(r_max))
    return mu, r


def _nb_p0_from_mu_r(mu, r, eps=1e-12):
    mu = np.asarray(mu, dtype=np.float64)
    r = np.asarray(r, dtype=np.float64)
    p = r / (r + mu + eps)
    return p**r


def _init_alpha_from_zeros_nb_weight(counts_dense, mu0, r0, eps=1e-12):
    X = np.asarray(counts_dense, dtype=np.int32)
    n_cells = X.shape[0]
    z = (X == 0).sum(axis=0).astype(np.float64) / float(n_cells)

    mu0 = np.asarray(mu0, dtype=np.float64)
    r0 = np.asarray(r0, dtype=np.float64)

    p0 = _nb_p0_from_mu_r(mu0, r0, eps=eps)
    denom = 1.0 - p0
    denom = np.maximum(denom, eps)

    alpha = (1.0 - z) / denom
    alpha = np.clip(alpha, 1e-8, 1.0 - 1e-8)
    return alpha


def fit_nb_marginals(
    counts_dense,
    zi=False,
    r_init=10.0,
    zi_lr=1e-2,
    max_iter=1000,
    zi_alpha_init=0.05,
    zi_alpha_init_mode="moments",
):

    n_genes = counts_dense.shape[1]

    if zi_alpha_init_mode == "moments":
        mu0, r0 = _init_nb_moments(counts_dense)
        r0 = np.where(np.isfinite(r0), r0, float(r_init))
        alpha0 = _init_alpha_from_zeros_nb_weight(counts_dense, mu0, r0)
    elif zi_alpha_init_mode == "zeros":
        mu0 = counts_dense.mean(axis=0, dtype=np.float64)
        mu0[mu0 <= 0] = 1e-12
        r0 = np.full(n_genes, float(r_init), dtype=np.float64)
        if zi_alpha_init is None:
            alpha0 = _init_alpha_from_zeros_nb_weight(counts_dense, mu0, r0)
        else:
            alpha0 = np.full(n_genes, float(zi_alpha_init), dtype=np.float64)
    else:
        raise ValueError("zi_alpha_init_mode must be 'moments' or 'zeros'.")

    if not zi:
        m0_all, r_all = optimise_all_genes(
            counts_dense.T.astype(np.int32, copy=False),
            list(np.asarray(mu0, dtype=np.float64)),
            list(np.asarray(r0, dtype=np.float64)),
            int(max_iter),
        )

        return (
            np.asarray(m0_all, dtype=np.float64),
            np.asarray(r_all, dtype=np.float64),
        )

    else:
        m0_all, r_all, alpha_all = optimise_all_genes_zi(
            counts_dense.T.astype(np.int32, copy=False),
            list(np.asarray(mu0, dtype=np.float64)),
            list(np.asarray(r0, dtype=np.float64)),
            list(np.asarray(alpha0, dtype=np.float64)),
            float(zi_lr),
            int(max_iter),
        )

        return (
            np.asarray(m0_all, dtype=np.float64),
            np.asarray(r_all, dtype=np.float64),
            np.clip(np.asarray(alpha_all, dtype=np.float64), 1e-8, 1.0 - 1e-8),
        )


if __name__ == "__main__":
    d = np.random.randint(0, 100, size=(1000, 100))

    ms, rs = fit_nb_marginals(d)
