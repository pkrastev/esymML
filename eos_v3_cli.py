#!/usr/bin/env python
"""
Program: eos_v3_cli.py
         Generate parametric EOSs via Monte-Carlo sampling
         We vary J0, L, Ksym, and Jsym

         Command-line (and optimized) version of eos_v3.py:
         python eos_v3_cli.py -n <number of instances to try> -o <output file>
                              [-s <random seed>] [--paper]

         --paper: setup of the Galaxies 2022 paper: J0 = Jsym = 0 (only L and
         Ksym are sampled), output tables up to 10**0.015 fm^-3, no check on
         the pressure at the top of the table (skip 4), eospar = [L, Ksym].

         The beta-equilibrium solve runs in Numba (beta_numba.py, a port of
         the MINPACK solver behind scipy.optimize.root). For a given random
         seed the output is identical to eos_v3.py (with --paper: to eos_v3.py
         changed in the same way).
"""
import h5py
import time
import random
import argparse
import numpy as np

from scipy.interpolate import interp1d

from beta_numba import solve_beta

# --- Parameters ---
# --- (i) physical / mathematical constants ---
mn   = 939.5656328            # neutron mass
mp   = 938.272328             # proton mass
mm   = 105.658389             # muon mass
pi   = 3.1415926535897932384  # pi
rho0 = 0.16                   # density of normal nuclear matter

# --- (ii) conversion factors ---
ev           = 1.60217733e-12          # erg  (eV to erg)
yr           = 3.155693e7              # s    (years to seconds)
fm3cm3       = 1.0e+39                 # fm-3 to cm-3
mevfm3gcm3   = (1.0/1.7827)*1.0e-12    # MeV fm-3 to g cm-2
dycm2mevfm3  = (197.33/3.1616)*1.0e-35 #
ergcm3mevfm3 = (ev*1.0e6)*fm3cm3       #
hbar_c       = 197.33                  # h_bar * c

# --- EOS parameters (sampling ranges) ---
# --- Symmetric nuclear matter ---
J0_low    = -800.0
J0_high   = 400.0
# --- Symmetry Energy ---
L_low     = 30.6
L_high    = 86.8
Ksym_low  = -400.0
Ksym_high = 100.0
Jsym_low  = -200.0
Jsym_high = 800.0

# --- Constants used in the beta-equilibrium equations ---
dm    = mn - mp
mm2   = mm**2
hbc2  = hbar_c**2
third = 1./3.
two3  = 2./3.

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

def derivative(func, x0, dx):
    """First derivative with a 5-point central difference.
    Replaces scipy.misc.derivative(func, x0, dx, n=1, order=5), which was
    removed from SciPy; same weights and summation order."""
    weights = np.array([1, -8, 0, 8, -1]) / 12.0
    val = 0.0
    for k in range(5):
        if weights[k] != 0.0:
            val += weights[k] * func(x0 + (k - 2) * dx)
    return val / dx

def log_grid(dlo, dhi, nstep):
    """Points equally spaced in log10 between 10**dlo and 10**dhi"""
    zzx = np.zeros(nstep)
    dstp = 0.0
    if ( nstep != 1 ):
        dstp = ( dhi - dlo ) / ( nstep - 1)
    for j in np.arange(nstep):
        zz = dlo + float(j)*dstp
        zz = 10.0**zz
        zzx[j] = zz
    return zzx

# --- Density grids (the same for every EOS) ---
rhox = np.linspace(0.02,1.35, 500)       # grid for beta-equilibrium
xx0  = np.linspace(0.04, 1.34, 500)      # grid for pressure / energy density
kk   = np.array(np.where(xx0 <= 0.09)).size
zzx  = log_grid(-14.2, 1.121, 500)       # grid for the output tables (--paper: up to 10**0.015)
rlog = np.log10(log_grid(-14.0, 0.09, 500))  # grid for the speed of sound
lzzx = np.log10(zzx)

# --- Per-density factors for the beta-equilibrium equations ---
tx   = np.array([3.0*pi*pi*rhox[i3] for i3 in range(len(rhox))])
kfex = np.array([hbar_c*( ( 3.0*pi*pi*rhox[i3] )**(1./3.)) for i3 in range(len(rhox))])
cst  = np.array([dm, mm2, hbc2, third, two3, 0.5])

# --- Crust (APR EOS), set by load_crust() ---
crust = None

# --- Check the pressure at the top of the table (skip 4); off with --paper ---
check_ptop = True

def load_crust(fname="EOS_P_A18dvUIX.txt"):
    """Load the APR EOS (Rodrigo pipeline) used below 0.075 fm^-3"""
    darr = np.loadtxt(fname, skiprows=0)
    eden_apr = 10.0**darr[:,0]*mevfm3gcm3
    p_apr    = 10.0**darr[:,1]/ergcm3mevfm3
    rho_apr  = darr[:,2]/fm3cm3
    jj = np.array(np.where(rho_apr < 0.075)).size
    return rho_apr[:jj], p_apr[:jj], eden_apr[:jj]

def compute_eos(par):
    """Build one EOS table from (J0, L, Ksym, Jsym).
    Returns (status, eos): status is 'ok' or the reason the EOS was rejected."""
    J0, L, Ksym, Jsym = par
    rho_apr, p_apr, eden_apr = crust

    e0x   = e_0(rhox, J0)
    esymx = e_sym(rhox, L, Ksym, Jsym)

    if (e0x[-1] > e0x[-1] + esymx[-1]) or e0x[-1] < 0.0 or esymx[-1] < 0.0:
        return 'skip 1', None

    # --- Particle concentrations ---
    esym  = esymx
    rhozz = rhox
    ok, ypr, yn, ye, ymu = solve_beta(esym, rhozz, tx, kfex, cst)
    if not ok:
        return 'solver failure', None

    if yn[-1] == yn[-2] == 1.0:
        return 'skip 2', None

    # --- beta-stable matter (total energy) ---
    e_beta   = e0x + esym * ( 1.0 - 2.0*ypr )**2 + mn*yn + mp*ypr
    eel_beta = (hbar_c*( ( 3.0*pi*pi*ye)**(4./3.))/(4.*pi*pi))*rhozz**(1./3.)
    kfmu     = hbar_c*( 3.0*pi*pi*rhozz*ymu )**(1./3.)
    mumu     = np.sqrt(mm*mm + hbar_c*hbar_c*(3.0*pi*pi*rhozz*ymu)**(2./3.))
    emu_beta = (mumu*kfmu*(mumu*mumu - 0.5*mm*mm) - 0.5*mm**4.0*np.log((mumu + kfmu)/mm))/rhozz / (4.0*pi*pi*hbar_c**3.)
    etot     = e_beta + eel_beta + emu_beta          # Total energy

    # --- Total pressure ---
    f  = interp1d(rhozz, etot, kind = 'quadratic')
    xx = xx0
    df = derivative(f, xx, dx=1e-6)
    pz = xx * xx * df

    # --- Energy density ---
    eden_tmp = rhozz*etot
    f2       = interp1d(rhozz, eden_tmp, kind = 'quadratic')
    edenz    = f2(xx)

    pz = pz[kk:]
    edenz = edenz[kk:]
    xx = xx[kk:]

    rhoxx  = np.concatenate((rho_apr, xx))
    pxx    = np.concatenate((p_apr, pz))
    edenxx = np.concatenate((eden_apr, edenz))

    # --- Interpolated values ---
    eden = np.interp(zzx,rhoxx,edenxx)
    p    = np.interp(zzx,rhoxx,pxx)

    if np.any(p < 0.0):
        return 'skip 3', None

    if check_ptop and p[-1] < 1200.0:
        return 'skip 4', None

    # --- Speed of sound ---
    r  = rlog
    rh = 1e-6
    rp = r + rh
    rm = r - rh
    pp = np.interp(rp,lzzx,np.log10(p))
    pm = np.interp(rm,lzzx,np.log10(p))
    edp = np.interp(rp,lzzx,np.log10(eden))
    edm = np.interp(rm,lzzx,np.log10(eden))
    pp = 10.0**pp
    pm = 10.0**pm
    edp = 10.0**edp
    edm = 10.0**edm

    # Simple 2-point derivative approxiamtion...........................
    cs2 = ( pp - pm ) / (edp - edm )

    if np.any(cs2 < 0.0):
        return 'skip 5', None

    # --- Write final tables ---
    eos = np.zeros((3,500))
    eos[0,:] = eden/mevfm3gcm3  # Energy density
    eos[1,:] = p/dycm2mevfm3    # Pressure
    eos[2,:] = zzx*fm3cm3       # Number density

    return 'ok', eos

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Generate parametric EOSs via Monte-Carlo sampling")
    parser.add_argument("-n", "--ntry", type=int, default=100000,
                        help="number of instances (Monte-Carlo samples) to try (default: 100000)")
    parser.add_argument("-o", "--output", type=str, default="eos_set_test.h5",
                        help="name of the output HDF5 file (default: eos_set_test.h5)")
    parser.add_argument("--paper", action="store_true",
                        help="paper setup: J0 = Jsym = 0, tables up to 10**0.015 fm^-3, "
                             "no skip 4, eospar = [L, Ksym]")
    parser.add_argument("-s", "--seed", type=int, default=None,
                        help="random seed, for reproducible runs (default: none)")
    args = parser.parse_args()
    if args.ntry < 1:
        parser.error("--ntry must be a positive integer")
    return args

def main():
    args = parse_args()
    t_start = time.time()

    global crust, zzx, lzzx, check_ptop
    crust = load_crust()
    if args.paper:
        zzx  = log_grid(-14.2, 0.015, 500)
        lzzx = np.log10(zzx)
        check_ptop = False

    # --- Monte Carlo sampling of the EOS parameters ---
    # (same draw order as eos_v3.py: J0, Jsym, L, Ksym)
    random.seed(args.seed)
    pars = []
    for i in range(args.ntry):
        if args.paper:
            # --- only L and Ksym are sampled (in this order) ---
            J0 = 0.0
            Jsym = 0.0
        else:
            J0 = random.uniform(J0_low, J0_high)
            Jsym = random.uniform(Jsym_low, Jsym_high)
        L = random.uniform(L_low, L_high)
        Ksym = random.uniform(Ksym_low, Ksym_high)
        pars.append((J0, L, Ksym, Jsym))

    print('Trying {0} instances{1}'.format(args.ntry, ' (paper setup)' if args.paper else ''), flush=True)

    eos_arr = []          # --> EOS
    eos_par_arr = []      # --> EOS parameters (J0, L, Ksym, Jsym)
    counts = {}
    nprint = max(1, args.ntry // 20)

    for i, par in enumerate(pars):
        status, eos = compute_eos(par)
        counts[status] = counts.get(status, 0) + 1
        if status == 'ok':
            eos_arr.append(eos)
            if args.paper:
                eos_par_arr.append(np.array(pars[i][1:3]))    # [L, Ksym]
            else:
                eos_par_arr.append(np.array(pars[i]))         # [J0, L, Ksym, Jsym]
        if (i + 1) % nprint == 0:
            print('iter: {0}/{1}  accepted: {2}  ({3:.1f} s)'.format(
                i + 1, args.ntry, len(eos_arr), time.time() - t_start), flush=True)

    # --- END OF LOOP OVER NUMBER OF SAMPLES ---

    eos_par_arr = np.array(eos_par_arr)
    eos_arr = np.array(eos_arr)

    # --- Write data ---
    h5f   = h5py.File(args.output,'w')
    data1 = h5f.create_dataset("eos", data=eos_arr, dtype='f8')
    data2 = h5f.create_dataset("eospar", data=eos_par_arr, dtype='f8')
    h5f.close()

    print('Rejected:', ', '.join('{0}: {1}'.format(k, v) for k, v in sorted(counts.items()) if k != 'ok'))
    print('Accepted {0} of {1} instances; written to {2}  ({3:.1f} s)'.format(
        len(eos_arr), args.ntry, args.output, time.time() - t_start))

if __name__ == "__main__":
    main()
