"""
Module: esym_data.py
        Data preparation for the Esym(rho) DNNs, shared by
        train_esym_dnn.py, MR_EOS_v3.ipynb and figures_v2.ipynb

        Follows MR_EOS_v2.ipynb / figures.ipynb:
        - NS sequences from tov_ml_cli.x (ns_paper_v*.h5): M, R, lambda
        - the same quality filter on the sequences
        - DNN input: 50 random masses in [1, 2] Msun (sorted) followed by
          R(M) or Lambda(M) at those masses (100 values)
        - DNN output: Esym(rho) at 100 densities in [0.07, 0.8] fm^-3

        Changes with respect to the notebooks:
        - the random masses come from a seeded generator (reproducible)
        - the curves are interpolated on the stable branch (M up to Mmax),
          where M(rho_c) is increasing as np.interp requires
        - the split is a seeded random permutation: n_train EOSs for
          training, the rest divided equally into validation and test
        - the input scaler is fitted on the training set only
"""
import h5py
import numpy as np

# --- Constants (cgs), as in the notebooks ---
c     = 2.99792458e10          # Speed of light
g     = 6.67259e-8             # Gravitational constant
m_sun = 1.9892e33              # Mass of the Sun
rho0  = 0.16                   # density of normal nuclear matter [fm^-3]

# --- Conversion factors, as in the notebooks ---
fm3cm3      = 1.0e+39                  # fm-3 to cm-3
mevfm3gcm3  = (1.0/1.7827)*1.0e-12     # MeV fm-3 to g cm-3
dycm2mevfm3 = (197.33/3.1616)*1.0e-35  # dyn cm-2 to MeV fm-3

# --- DNN output grid ---
RHOX = np.linspace(0.07, 0.8, 100)     # [fm^-3]; RHOX[-1] = 0.8 fm^-3 = 5 rho0

# --- Data files from this project ---
EOS_FILES = ['eos_paper_v1.h5', 'eos_paper_v2.h5', 'eos_paper_v3.h5', 'eos_paper_v4.h5']
NS_FILES  = ['ns_paper_v1.h5', 'ns_paper_v2.h5', 'ns_paper_v3.h5', 'ns_paper_v4.h5']

N_POINTS = 50      # points per input curve
M_LOW    = 1.0     # input mass range [Msun]
M_HIGH   = 2.0


def e_0(rho, J_0=0.0):
    """Energy of symmetric nuclear matter"""
    E_sat = -15.8              # Energy per nucleon [MeV]
    K_0   = 240.0              # K_0 [MeV]
    rho_0 = 0.16               # saturation density [fm^-1]
    x     = (rho - rho_0) / (3.0*rho_0)
    E     = E_sat + (1.0/2.0)*K_0*(x**2.0) + (1.0/6.0)*J_0*(x**3.0)
    return E


def e_sym(rho, L=58.7, K_sym=0.0, J_sym=0.0):
    """Symmetry energy"""
    E_sym0 = 31.7            # Symmetry energy at saturation density [MeV]
    rho_0 = 0.16             # saturation density [fm^-1]
    x     = (rho - rho_0) / (3.0*rho_0)
    E     = E_sym0 + L*x + (1.0/2.0)*K_sym*(x**2.0) + (1.0/6.0)*J_sym*(x**3.0)
    return E


def rc_params():
    """Set rcParams for the plots (as in figures.ipynb)"""
    import matplotlib as mpl
    params = {
        'axes.linewidth':1.5,
        'xtick.major.size':6.5,
        'xtick.major.width':1.5,
        'xtick.minor.size':3.5,
        'xtick.minor.width':1.5,
        'ytick.major.size':6.5,
        'ytick.major.width':1.5,
        'ytick.minor.size':3.5,
        'ytick.minor.width':1.5,
        'xtick.labelsize':14,
        'ytick.labelsize':14,
        'xtick.direction':'in',
        'ytick.direction':'in',
        'xtick.top':True,
        'ytick.right':True
        }
    mpl.rcParams.update(params)


def lambda_dimensionless(m, l):
    """Lambda = lambda / (G M / c^2)^5 * c^2 / G, for M in Msun and
    lambda in 1e36 g cm^2 s^2 (as in the notebooks)"""
    ms = m * m_sun
    cl = ((g**4.0) / (c**10.0)) * (ms**5.0) * 1e-36
    return l / cl


def load_ns(files=NS_FILES):
    """NS sequences and EOS parameters, concatenated in file order.
    Returns a dict with m, r, l (lambda), eospar and the column names."""
    m, r, l, par = [], [], [], []
    cols = None
    for fn in files:
        with h5py.File(fn, 'r') as f:
            ns = f['ns'][:]
            m.append(ns[:, 0, :])
            r.append(ns[:, 1, :])
            l.append(ns[:, 2, :])
            par.append(f['eospar'][:])
            c_ = f['eospar'].attrs.get('columns', b'')
            c_ = c_.decode() if isinstance(c_, bytes) else str(c_)
            cols = [s.strip() for s in c_.split(',')]
    return {'m': np.concatenate(m), 'r': np.concatenate(r), 'l': np.concatenate(l),
            'eospar': np.concatenate(par), 'columns': cols}


def quality_mask(m, r, l):
    """The filter of MR_EOS_v2.ipynb / figures.ipynb (True = keep):
    reject a sequence if R jumps up by more than 0.03 km, lambda > 5.5 or
    lambda jumps up by more than 1.0 between neighbouring points, if it has
    NaN masses or negative lambda, or if its maximum mass is below 2 Msun."""
    bad = (np.diff(r, axis=1) > 0.03).any(axis=1)
    bad |= (l[:, :-1] > 5.5).any(axis=1)
    bad |= (np.diff(l, axis=1) > 1.0).any(axis=1)
    bad |= np.isnan(m).any(axis=1)
    bad |= (l < 0.0).any(axis=1)
    bad |= ~(np.nanmax(np.where(np.isnan(m), -np.inf, m), axis=1) >= 2.0)
    return ~bad


def esym_params(eospar, columns):
    """(L, Ksym, Jsym) arrays from eospar (2 columns [L, Ksym] or 4 columns
    [J0, L, Ksym, Jsym])"""
    if eospar.shape[1] == 2:
        return eospar[:, 0], eospar[:, 1], np.zeros(len(eospar))
    if eospar.shape[1] == 4:
        return eospar[:, 1], eospar[:, 2], eospar[:, 3]
    raise ValueError('eospar must have 2 or 4 columns, got %d (%s)' % (eospar.shape[1], columns))


def load_dataset(files=NS_FILES):
    """Filtered sequences: m, r, l, Lambda (N, 50) and L, Ksym, Jsym (N,)"""
    d = load_ns(files)
    keep = quality_mask(d['m'], d['r'], d['l'])
    L, Ksym, Jsym = esym_params(d['eospar'], d['columns'])
    out = {k: d[k][keep] for k in ('m', 'r', 'l')}
    out['Lambda'] = lambda_dimensionless(out['m'], out['l'])
    out['L'], out['Ksym'], out['Jsym'] = L[keep], Ksym[keep], Jsym[keep]
    out['n_total'] = len(keep)
    out['n_kept'] = int(keep.sum())
    return out


def split_indices(n, n_train=40000, seed=42):
    """Seeded random split: n_train for training, the rest divided equally
    into validation and test"""
    if n <= n_train + 2:
        raise ValueError('only %d sequences, need more than n_train = %d' % (n, n_train))
    perm = np.random.default_rng(seed).permutation(n)
    n_val = (n - n_train) // 2
    return perm[:n_train], perm[n_train:n_train+n_val], perm[n_train+n_val:]


def stable_branch(m, y):
    """Part of the sequence up to the maximum mass"""
    k = int(np.argmax(m))
    return m[:k+1], y[:k+1]


def make_inputs(data, idx, kind='MR', seed=0):
    """DNN inputs and outputs for the sequences idx.
    kind = 'MR': x = [M_1..M_50, R(M_1)..R(M_50)]
    kind = 'MLambda': x = [M_1..M_50, Lambda(M_1)..Lambda(M_50)]
    y = Esym(RHOX) for the EOS parameters of each sequence."""
    ykey = {'MR': 'r', 'MLambda': 'Lambda'}[kind]
    rng = np.random.default_rng(seed)
    x = np.zeros((len(idx), 2*N_POINTS))
    y = np.zeros((len(idx), len(RHOX)))
    for n, i in enumerate(idx):
        mx = np.sort(rng.uniform(low=M_LOW, high=M_HIGH, size=(N_POINTS,)))
        ms, ys = stable_branch(data['m'][i], data[ykey][i])
        x[n, :N_POINTS] = mx
        x[n, N_POINTS:] = np.interp(mx, ms, ys)
        y[n, :] = e_sym(RHOX, data['L'][i], data['Ksym'][i], data['Jsym'][i])
    return x, y


class Scaler:
    """Min-max scaling of each input column, fitted on the training set
    (same transform as sklearn's MinMaxScaler: x*scale + min)"""

    def __init__(self, scale=None, min_=None):
        self.scale = scale
        self.min_ = min_

    def fit(self, x):
        from sklearn.preprocessing import MinMaxScaler
        s = MinMaxScaler().fit(x)
        self.scale, self.min_ = s.scale_.copy(), s.min_.copy()
        return self

    def transform(self, x):
        return x*self.scale + self.min_

    def save(self, fn):
        np.savez(fn, scale=self.scale, min_=self.min_)

    @classmethod
    def load(cls, fn):
        d = np.load(fn)
        return cls(d['scale'], d['min_'])


def build_sets(kind='MR', files=NS_FILES, n_train=40000, seed=42):
    """Train / validation / test sets for kind = 'MR' or 'MLambda'.
    Uses the seeds seed (split), seed+1, seed+2, seed+3 (random masses)."""
    data = load_dataset(files)
    idx = split_indices(data['n_kept'], n_train, seed)
    sets = {'data': data, 'kind': kind}
    for name, ix, s in zip(('train', 'val', 'test'), idx, (seed+1, seed+2, seed+3)):
        x, y = make_inputs(data, ix, kind, s)
        sets['x_' + name + '_raw'] = x
        sets['y_' + name] = y
        sets['idx_' + name] = ix
    scaler = Scaler().fit(sets['x_train_raw'])
    for name in ('train', 'val', 'test'):
        sets['x_' + name] = scaler.transform(sets['x_' + name + '_raw'])
    sets['scaler'] = scaler
    return sets
