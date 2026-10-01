# esymML

**Nuclear symmetry energy from neutron-star observations with deep neural networks.**

Code for P. G. Krastev, *Translating Neutron Star Observations to Nuclear Symmetry Energy via
Deep Neural Networks*, Galaxies **10**, 16 (2022),
[doi:10.3390/galaxies10010016](https://doi.org/10.3390/galaxies10010016).

A deep neural network (DNN) learns the density dependence of the nuclear symmetry energy,
E<sub>sym</sub>(ρ) up to 5ρ<sub>0</sub>, from neutron-star observables: the mass–radius relation
M(R), or the mass–tidal-deformability relation M(Λ). The training data are built from a
parametric equation of state (EOS):

```
 parametric EOS          neutron stars             DNN training               figures
 eos_v3_cli.x / .py  ->  tov_ml_cli.x          ->  train_esym_dnn.py      ->  figures_v2.ipynb
 eos_*.h5               ns_*.h5                   models/                     figures/
 (L, Ksym sampled;      (M, R, lambda, k2, I,     (M-R and M-Lambda           (MR_EOS_v3.ipynb:
  beta-stable npe-mu)    beta for 50 stars/EOS)    input)                      training analysis)
```

## Contents

| File | Purpose |
|---|---|
| `eos_v3_cli.f90`, `minpack_hybrd.f`, `Makefile.eos` | EOS generator, Fortran (`eos_v3_cli.x`) |
| `eos_v3_cli.py`, `beta_numba.py` | EOS generator, Python; same output as the Fortran version |
| `EOS_P_A18dvUIX.txt` | APR crust table, used below 0.075 fm<sup>−3</sup> (read at run time) |
| `tov_ml_cli.f90`, `libtov.f90`, `libtov_h.f90`, `Makefile.tov` | Neutron-star sequences from an EOS file (`tov_ml_cli.x`); solver from [tovSolve](https://github.com/pkrastev/tovSolve) |
| `esym_data.py` | Data preparation for the DNNs (filter, inputs, split, scaler) |
| `train_esym_dnn.py` | DNN training (M–R or M–Λ input) |
| `MR_EOS_v3.ipynb` | Training notebook: data, loss curves, test errors, predictions |
| `figures_v2.ipynb` | The figures of the paper, from the new data and models |
| `models/` | Trained DNNs, input scalers, test sets, loss histories, metrics |
| `figures/` | Figures made by `figures_v2.ipynb` |
| `make_data.sh` | Regenerates the data files (EOS tables and NS sequences) |
| `environment.yml` | Python environment |

## Requirements

- Linux (x86-64, glibc).
- **Fortran:** Intel `ifx` (oneAPI) and an HDF5 library with Fortran support built with `ifx`.
  The Makefiles use the environment variables `HDF5_INCLUDE` and `HDF5_LIB`
  (set, e.g., by `module load hdf5` on HPC systems).
- **Python:** the conda environment in `environment.yml` (Python 3.10, NumPy 1.26, SciPy 1.12,
  h5py, numba, TensorFlow 2.17 / Keras 3, scikit-learn, pandas, matplotlib, Jupyter).
  `eos_v3_cli.x` links the same LAPACK (OpenBLAS) as SciPy in this environment.
- A GPU is optional. The models in `models/` were trained on an NVIDIA A100 (about 40 minutes per model).

## Build

```bash
conda env create -f environment.yml
conda activate esymml
module load intel hdf5                 # or: make ifx and an ifx-built HDF5 available,
                                       # with HDF5_INCLUDE and HDF5_LIB set
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

make -f Makefile.eos LAPACK="-L$CONDA_PREFIX/lib -llapack"     # eos_v3_cli.x
make -f Makefile.tov                                           # tov_ml_cli.x
```

On the FASRC Cannon cluster: `module load python`, activate the environment, then
`module load intel mpich hdf5`. Both programs build in the same directory, and
`make -f Makefile.<eos|tov> clean` removes only the files of that build.

## Usage

### 1. EOS tables

```bash
./eos_v3_cli.x -n 100000 --paper -s 101 -o eos_paper_v1.h5
./eos_v3_cli.x -n 100000 --paper -s 102 -o eos_paper_v2.h5
./eos_v3_cli.x -n 100000 --paper -s 103 -o eos_paper_v3.h5
./eos_v3_cli.x -n 100000 --paper -s 104 -o eos_paper_v4.h5
```

| Option | Meaning |
|---|---|
| `-n NTRY` | number of Monte Carlo samples to try (default 100000) |
| `-o OUTPUT` | output HDF5 file (default `eos_set_test.h5`) |
| `-s SEED` | random seed (default: none) |
| `--paper` | the setup of the paper: J<sub>0</sub> = J<sub>sym</sub> = 0, only L and K<sub>sym</sub> are sampled; tables up to 10<sup>0.015</sup> fm<sup>−3</sup> |

`python eos_v3_cli.py` takes the same options and writes the same file (see Reproducibility).
Both read `EOS_P_A18dvUIX.txt` from the current directory. About 15% of the samples are
accepted with `--paper`; one run of 100K samples takes about 2 minutes with `eos_v3_cli.x`.

Output: `eos` (N, 3, 500) with the energy density (g cm<sup>−3</sup>), pressure (dyn cm<sup>−2</sup>)
and baryon number density (cm<sup>−3</sup>) at 500 densities, logarithmically spaced from
10<sup>−14.2</sup> to 10<sup>1.121</sup> fm<sup>−3</sup> (10<sup>0.015</sup> with `--paper`), and `eospar`:
(N, 2) = [L, K<sub>sym</sub>] with `--paper`, (N, 4) = [J<sub>0</sub>, L, K<sub>sym</sub>, J<sub>sym</sub>] otherwise (MeV).

### 2. Neutron-star sequences

```bash
./tov_ml_cli.x -i eos_paper_v1.h5 -o ns_paper_v1.h5
./tov_ml_cli.x -i eos_paper_v2.h5 -o ns_paper_v2.h5
./tov_ml_cli.x -i eos_paper_v3.h5 -o ns_paper_v3.h5
./tov_ml_cli.x -i eos_paper_v4.h5 -o ns_paper_v4.h5
```

| Option | Meaning |
|---|---|
| `-i INPUT`, `-o OUTPUT` | EOS file, output file |
| `--rho-min`, `--rho-max`, `--nrho` | central baryon densities (default 0.13–1.02 fm<sup>−3</sup>, 50 stars) |
| `--tol` | relative tolerance of the ODE solver (default 10<sup>−8</sup>) |

The number of EOSs is read from the input file. Output: `ns` (N, 3, 50) = M (M<sub>☉</sub>),
R (km), λ (10<sup>36</sup> g cm<sup>2</sup> s<sup>2</sup>); `k2`, `I` (10<sup>45</sup> g cm<sup>2</sup>),
`beta` (N, 50); `rhoc` (fm<sup>−3</sup>); and `eospar`, copied from the input. A star whose
integration fails is stored as NaN. About 0.1 s per EOS (25–35 minutes per file above).

### 3. DNN training

```bash
python train_esym_dnn.py --input MR          # M-R input
python train_esym_dnn.py --input MLambda     # M-Lambda input
```

The script reads `ns_paper_v1–v4.h5` from the current directory (`--ns-files` to change them)
and writes to `models/`:
`esym_<input>.keras` and `.weights.h5` (model and weights), `_scaler.npz` (input scaler),
`_test.npz` (test set and predictions), `_history.csv`, `_metrics.json`.
Other options: `--n-train`, `--epochs`, `--patience`, `--batch-size`, `--lr`, `--seed`, `--out-dir`.

### 4. Notebooks

- `MR_EOS_v3.ipynb`: data, training curves, test errors and predictions. It loads the models
  in `models/`; set `RUN_TRAINING = True` to train in the notebook instead.
- `figures_v2.ipynb`: the figures of the paper, written to `figures/`.

Both need the `ns_paper_*.h5` files (and `figures_v2.ipynb` also `eos_paper_v1.h5`) in this directory.

## Physics and DNN setup

**EOS.** Energy per nucleon of asymmetric matter, E(ρ, δ) = E<sub>0</sub>(ρ) + E<sub>sym</sub>(ρ) δ², with
x = (ρ − ρ<sub>0</sub>)/3ρ<sub>0</sub>:

- E<sub>0</sub>(ρ) = E<sub>sat</sub> + K<sub>0</sub>x²/2 + J<sub>0</sub>x³/6, with E<sub>sat</sub> = −15.8 MeV, K<sub>0</sub> = 240 MeV, ρ<sub>0</sub> = 0.16 fm<sup>−3</sup>
- E<sub>sym</sub>(ρ) = E<sub>sym</sub>(ρ<sub>0</sub>) + Lx + K<sub>sym</sub>x²/2 + J<sub>sym</sub>x³/6, with E<sub>sym</sub>(ρ<sub>0</sub>) = 31.7 MeV
- sampled uniformly: L ∈ [30.6, 86.8] MeV, K<sub>sym</sub> ∈ [−400, 100] MeV
  (and J<sub>0</sub> ∈ [−800, 400], J<sub>sym</sub> ∈ [−200, 800] MeV without `--paper`)

Neutron-star matter is n, p, e, μ in β-equilibrium, with the APR crust below 0.075 fm<sup>−3</sup>.
An EOS is rejected if E<sub>0</sub> or E<sub>sym</sub> is negative at 1.35 fm<sup>−3</sup>, if the β-equilibrium
solve fails, for pure neutron matter at the highest densities, or for negative pressure or negative speed of sound squared
(without `--paper` also if the pressure at the top of the table is below 1200 MeV fm<sup>−3</sup>).

**Neutron stars.** Structure, tidal Love number and moment of inertia from tovSolve in the
pseudo-enthalpy formalism (Lindblom 1992), with an adaptive Dormand–Prince 5(4) integrator.

**DNNs.** Input: 50 random masses in 1–2 M<sub>☉</sub> (sorted) and R or Λ at those masses
(100 values, min–max scaled with the training set). Output: E<sub>sym</sub> at 100 densities in
0.07–0.8 fm<sup>−3</sup>. Network: 12 dense layers of 100 (ReLU; linear input and output layers),
Adam with AMSgrad, learning rate 0.003, MSE loss, batch size 500, at most 2000 epochs, early
stopping on the validation loss (patience 100, best weights restored). Before training, sequences
with jumps in R or λ, NaN, negative λ or a maximum mass below 2 M<sub>☉</sub> are removed.

## Results

The models in `models/` were trained on the four `ns_paper` files: 60,144 sequences, of which
49,357 pass the filter. They were split into 40,000 for training, 4,678 for validation and
4,679 for testing (seed 42). Absolute error of E<sub>sym</sub> on the test set:

| Input | at 5ρ<sub>0</sub> (mean ± std) | at 5ρ<sub>0</sub> (max) | all densities (mean) |
|---|---|---|---|
| M–R | 0.18 ± 0.16 MeV | 1.6 MeV | 0.07 MeV |
| M–Λ | 0.47 ± 0.44 MeV | 5.2 MeV | 0.16 MeV |

Both models ran the full 2000 epochs (best epochs 1991 and 1995): the validation loss was
still decreasing slowly, so early stopping did not end the training.

## Data

The data files are not in this repository (`eos_paper_v1–v4.h5`, about 180 MB each, and
`ns_paper_v1–v4.h5`, about 36 MB each). They are regenerated, after building the two programs, with

```bash
./make_data.sh
```

which runs these commands:

```bash
./eos_v3_cli.x -n 100000 --paper -s 101 -o eos_paper_v1.h5
./eos_v3_cli.x -n 100000 --paper -s 102 -o eos_paper_v2.h5
./eos_v3_cli.x -n 100000 --paper -s 103 -o eos_paper_v3.h5
./eos_v3_cli.x -n 100000 --paper -s 104 -o eos_paper_v4.h5

./tov_ml_cli.x -i eos_paper_v1.h5 -o ns_paper_v1.h5
./tov_ml_cli.x -i eos_paper_v2.h5 -o ns_paper_v2.h5
./tov_ml_cli.x -i eos_paper_v3.h5 -o ns_paper_v3.h5
./tov_ml_cli.x -i eos_paper_v4.h5 -o ns_paper_v4.h5
```

This takes about 2–2.5 hours on one core (about 2 minutes per EOS file, 25–35 minutes per NS file).
The files used to train the models in `models/` contain 15,070, 14,925, 15,181 and 14,968 EOSs;
`make_data.sh` checks these numbers and warns if they differ.

With the toolchain described above (Intel `ifx` with the flags in the Makefiles, glibc, and the
environment in `environment.yml`, which pins OpenBLAS 0.3.26), the regenerated files are
bit-for-bit identical to the ones the models were trained on. With a different compiler or
LAPACK library the tables can differ in the last bits; the results are then the same within
numerical precision.

## Reproducibility

- For the same seed, `eos_v3_cli.x` and `eos_v3_cli.py` give bit-for-bit identical files when
  NumPy's AVX-512 math is disabled:
  `NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR"`.
  With AVX-512, NumPy computes `pow`, `log` and `log10` with different rounding, which changes the
  last bits of the tables (up to ~3×10<sup>−9</sup> relative) but not which EOSs are accepted.
  The Fortran version needs its strict floating-point flags (`Makefile.eos`), takes `pow`, `log`
  and `log10` from glibc at run time, and must be linked with the same LAPACK as SciPy.
- `libtov.f90` is compiled with `-no-vec` (`Makefile.tov`), so that an EOS read from HDF5 gives
  the same stars as the same EOS read from a text file by tovSolve.
- Relation to the code used for the paper: the paper's data were made with an earlier version of
  the EOS generator (`eos_v3.py`), a different TOV solver (`tov_ml.f90`) and the notebooks
  `MR_EOS_v2.ipynb` and `figures.ipynb`, which are not included here; the comments in the code
  refer to them. `--paper` reproduces the physical setup of the paper (model, constants,
  parameter ranges, density grid), but not its exact set of EOSs, because the script that
  generated them was not preserved. The TOV solver here is
  [tovSolve](https://github.com/pkrastev/tovSolve) (pseudo-enthalpy formalism).

## Citation

If you use this code, please cite the paper (see also `CITATION.cff`):

```bibtex
@article{Krastev2022,
  author  = {Krastev, Plamen G.},
  title   = {Translating Neutron Star Observations to Nuclear Symmetry Energy via Deep Neural Networks},
  journal = {Galaxies},
  volume  = {10},
  number  = {1},
  pages   = {16},
  year    = {2022},
  doi     = {10.3390/galaxies10010016}
}
```

## License

MIT License, see `LICENSE`. The repository contains code derived from MINPACK, NumPy, SciPy and
the Mersenne Twister, under their own licenses; see `THIRD_PARTY_NOTICES.md`. The crust table
`EOS_P_A18dvUIX.txt` is based on the APR equation of state (A. Akmal, V. R. Pandharipande and
D. G. Ravenhall, Phys. Rev. C **58**, 1804 (1998)).
