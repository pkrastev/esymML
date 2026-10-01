"""
Module: beta_numba.py
        Numba version of the beta-equilibrium solve used in eos_v3.py.

        hybrd() and its helpers are a line-by-line port of MINPACK as shipped
        with SciPy 1.12 (scipy/optimize/minpack/*.f), called with the same
        options as scipy.optimize.root(method='hybr') uses by default, so the
        particle fractions are bit-for-bit identical to the SciPy results.
        Indices in the ported routines are kept 1-based (as in the Fortran)
        and shifted by one on array access.
"""
import math
import numpy as np
from numba import njit

# --- MINPACK machine constants (dpmpar) ---
EPSMCH = 2.22044604926e-16    # dpmpar(1)
GIANT  = 1.79769313485e+308   # dpmpar(3)
RDWARF = 3.834e-20
RGIANT = 1.304e19

# --- scipy.optimize.root(method='hybr') defaults ---
XTOL   = 1.49012e-08
FACTOR = 100.0
EPSFCN = 2.220446049250313e-16   # np.finfo(float).eps

# ---------------------------------------------------------------------
# Residuals of the beta-equilibrium equations (same expressions as
# frac_e / frac_emu in eos_v3.py). Return False where the original code
# gets a numpy warning (negative base of a fractional power, or a
# non-finite result), which rejects the EOS.
# cst = [mn - mp, mm**2, hbar_c**2, 1/3, 2/3, 0.5]
# ---------------------------------------------------------------------
@njit(cache=True)
def fcn(n, x, fvec, c4, t, kfe, cst):
    dm   = cst[0]
    mm2  = cst[1]
    hbc2 = cst[2]
    if n == 3:
        # --- electrons only ---
        if x[2] < 0.0:
            return False
        f1 = c4*( 1.0 - 2.0*x[0] ) + dm - kfe*math.pow(x[2], cst[3])
        if not math.isfinite(f1):
            return False
        fvec[0] = f1
        fvec[1] = x[0] - x[2]
        fvec[2] = 1.0 - x[0] - x[1]
    else:
        # --- electrons plus muons ---
        tx3 = t*x[3]
        if x[2] < 0.0 or tx3 < 0.0:
            return False
        mue = kfe*math.pow(x[2], cst[3])
        f1 = c4*(1.0 - 2.0*x[0]) + dm - mue
        f2 = math.pow(mm2 + hbc2*math.pow(tx3, cst[4]), cst[5]) - mue
        if not math.isfinite(f1 + f2):
            return False
        fvec[0] = f1
        fvec[1] = f2
        fvec[2] = x[0] - x[2] - x[3]
        fvec[3] = 1.0 - x[0] - x[1]
    return True

# ---------------------------------------------------------------------
# MINPACK routines
# ---------------------------------------------------------------------
@njit(cache=True)
def enorm(n, x):
    """Euclidean norm of x(1:n)"""
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    x1max = 0.0
    x3max = 0.0
    floatn = float(n)
    agiant = RGIANT/floatn
    for i in range(n):
        xabs = abs(x[i])
        if xabs > RDWARF and xabs < agiant:
            s2 = s2 + xabs*xabs
        elif xabs > RDWARF:
            if xabs > x1max:
                q = x1max/xabs
                s1 = 1.0 + s1*(q*q)
                x1max = xabs
            else:
                q = xabs/x1max
                s1 = s1 + q*q
        else:
            if xabs > x3max:
                q = x3max/xabs
                s3 = 1.0 + s3*(q*q)
                x3max = xabs
            elif xabs != 0.0:
                q = xabs/x3max
                s3 = s3 + q*q
    if s1 != 0.0:
        return x1max*math.sqrt(s1+(s2/x1max)/x1max)
    if s2 != 0.0:
        if s2 >= x3max:
            return math.sqrt(s2*(1.0+(x3max/s2)*(x3max*s3)))
        return math.sqrt(x3max*((s2/x3max)+(x3max*s3)))
    return x3max*math.sqrt(s3)

@njit(cache=True)
def fdjac1(n, x, fvec, fjac, wa1, c4, t, kfe, cst):
    """Forward-difference Jacobian (dense case, ml + mu + 1 >= n)"""
    eps = math.sqrt(max(EPSFCN, EPSMCH))
    for j in range(1, n+1):
        temp = x[j-1]
        h = eps*abs(temp)
        if h == 0.0:
            h = eps
        x[j-1] = temp + h
        if not fcn(n, x, wa1, c4, t, kfe, cst):
            return False
        x[j-1] = temp
        for i in range(1, n+1):
            fjac[i-1, j-1] = (wa1[i-1] - fvec[i-1])/h
    return True

@njit(cache=True)
def qrfac(m, n, a, rdiag, acnorm, wa):
    """QR factorization without pivoting"""
    for j in range(1, n+1):
        acnorm[j-1] = enorm(m, a[:, j-1])
        rdiag[j-1] = acnorm[j-1]
        wa[j-1] = rdiag[j-1]
    minmn = min(m, n)
    for j in range(1, minmn+1):
        ajnorm = enorm(m-j+1, a[j-1:, j-1])
        if ajnorm != 0.0:
            if a[j-1, j-1] < 0.0:
                ajnorm = -ajnorm
            for i in range(j, m+1):
                a[i-1, j-1] = a[i-1, j-1]/ajnorm
            a[j-1, j-1] = a[j-1, j-1] + 1.0
            for k in range(j+1, n+1):
                s = 0.0
                for i in range(j, m+1):
                    s = s + a[i-1, j-1]*a[i-1, k-1]
                temp = s/a[j-1, j-1]
                for i in range(j, m+1):
                    a[i-1, k-1] = a[i-1, k-1] - temp*a[i-1, j-1]
        rdiag[j-1] = -ajnorm

@njit(cache=True)
def qform(m, n, q, wa):
    """Accumulate the orthogonal matrix Q from its factored form"""
    minmn = min(m, n)
    for j in range(2, minmn+1):
        for i in range(1, j):
            q[i-1, j-1] = 0.0
    for j in range(n+1, m+1):
        for i in range(1, m+1):
            q[i-1, j-1] = 0.0
        q[j-1, j-1] = 1.0
    for l in range(1, minmn+1):
        k = minmn - l + 1
        for i in range(k, m+1):
            wa[i-1] = q[i-1, k-1]
            q[i-1, k-1] = 0.0
        q[k-1, k-1] = 1.0
        if wa[k-1] != 0.0:
            for j in range(k, m+1):
                s = 0.0
                for i in range(k, m+1):
                    s = s + q[i-1, j-1]*wa[i-1]
                temp = s/wa[k-1]
                for i in range(k, m+1):
                    q[i-1, j-1] = q[i-1, j-1] - temp*wa[i-1]

@njit(cache=True)
def dogleg(n, r, diag, qtb, delta, x, wa1, wa2):
    """Dogleg step (r is the packed upper triangle, stored by rows)"""
    jj = (n*(n + 1))//2 + 1
    for k in range(1, n+1):
        j = n - k + 1
        jp1 = j + 1
        jj = jj - k
        l = jj + 1
        s = 0.0
        for i in range(jp1, n+1):
            s = s + r[l-1]*x[i-1]
            l = l + 1
        temp = r[jj-1]
        if temp == 0.0:
            l = j
            for i in range(1, j+1):
                temp = max(temp, abs(r[l-1]))
                l = l + n - i
            temp = EPSMCH*temp
            if temp == 0.0:
                temp = EPSMCH
        x[j-1] = (qtb[j-1] - s)/temp
    for j in range(1, n+1):
        wa1[j-1] = 0.0
        wa2[j-1] = diag[j-1]*x[j-1]
    qnorm = enorm(n, wa2)
    if qnorm <= delta:
        return
    l = 1
    for j in range(1, n+1):
        temp = qtb[j-1]
        for i in range(j, n+1):
            wa1[i-1] = wa1[i-1] + r[l-1]*temp
            l = l + 1
        wa1[j-1] = wa1[j-1]/diag[j-1]
    gnorm = enorm(n, wa1)
    sgnorm = 0.0
    alpha = delta/qnorm
    if gnorm != 0.0:
        for j in range(1, n+1):
            wa1[j-1] = (wa1[j-1]/gnorm)/diag[j-1]
        l = 1
        for j in range(1, n+1):
            s = 0.0
            for i in range(j, n+1):
                s = s + r[l-1]*wa1[i-1]
                l = l + 1
            wa2[j-1] = s
        temp = enorm(n, wa2)
        sgnorm = (gnorm/temp)/temp
        alpha = 0.0
        if sgnorm < delta:
            bnorm = enorm(n, qtb)
            dq = delta/qnorm
            sd = sgnorm/delta
            temp = (bnorm/gnorm)*(bnorm/qnorm)*sd
            temp = temp - dq*(sd*sd) + math.sqrt((temp-dq)*(temp-dq) + (1.0-dq*dq)*(1.0-sd*sd))
            alpha = (dq*(1.0 - sd*sd))/temp
    temp = (1.0 - alpha)*min(sgnorm, delta)
    for j in range(1, n+1):
        x[j-1] = temp*wa1[j-1] + alpha*x[j-1]

@njit(cache=True)
def r1updt(m, n, s, u, v, w):
    """Rank-1 update of the packed lower-trapezoidal matrix s"""
    jj = (n*(2*m - n + 1))//2 - (m - n)
    l = jj
    for i in range(n, m+1):
        w[i-1] = s[l-1]
        l = l + 1
    nm1 = n - 1
    for nmj in range(1, nm1+1):
        j = n - nmj
        jj = jj - (m - j + 1)
        w[j-1] = 0.0
        if v[j-1] != 0.0:
            if abs(v[n-1]) < abs(v[j-1]):
                cotan = v[n-1]/v[j-1]
                sin = 0.5/math.sqrt(0.25+0.25*(cotan*cotan))
                cos = sin*cotan
                tau = 1.0
                if abs(cos)*GIANT > 1.0:
                    tau = 1.0/cos
            else:
                tan = v[j-1]/v[n-1]
                cos = 0.5/math.sqrt(0.25+0.25*(tan*tan))
                sin = cos*tan
                tau = sin
            v[n-1] = sin*v[j-1] + cos*v[n-1]
            v[j-1] = tau
            l = jj
            for i in range(j, m+1):
                temp = cos*s[l-1] - sin*w[i-1]
                w[i-1] = sin*s[l-1] + cos*w[i-1]
                s[l-1] = temp
                l = l + 1
    for i in range(1, m+1):
        w[i-1] = w[i-1] + v[n-1]*u[i-1]
    sing = False
    for j in range(1, nm1+1):
        if w[j-1] != 0.0:
            if abs(s[jj-1]) < abs(w[j-1]):
                cotan = s[jj-1]/w[j-1]
                sin = 0.5/math.sqrt(0.25+0.25*(cotan*cotan))
                cos = sin*cotan
                tau = 1.0
                if abs(cos)*GIANT > 1.0:
                    tau = 1.0/cos
            else:
                tan = w[j-1]/s[jj-1]
                cos = 0.5/math.sqrt(0.25+0.25*(tan*tan))
                sin = cos*tan
                tau = sin
            l = jj
            for i in range(j, m+1):
                temp = cos*s[l-1] + sin*w[i-1]
                w[i-1] = -sin*s[l-1] + cos*w[i-1]
                s[l-1] = temp
                l = l + 1
            w[j-1] = tau
        if s[jj-1] == 0.0:
            sing = True
        jj = jj + (m - j + 1)
    l = jj
    for i in range(n, m+1):
        s[l-1] = w[i-1]
        l = l + 1
    if s[jj-1] == 0.0:
        sing = True
    return sing

@njit(cache=True)
def r1mpyq(m, n, a, v, w):
    """Multiply a (m x n) by the Givens rotations stored in v and w"""
    nm1 = n - 1
    if nm1 < 1:
        return
    cos = 0.0
    sin = 0.0
    for nmj in range(1, nm1+1):
        j = n - nmj
        if abs(v[j-1]) > 1.0:
            cos = 1.0/v[j-1]
            sin = math.sqrt(1.0-cos*cos)
        if abs(v[j-1]) <= 1.0:
            sin = v[j-1]
            cos = math.sqrt(1.0-sin*sin)
        for i in range(1, m+1):
            temp = cos*a[i-1, j-1] - sin*a[i-1, n-1]
            a[i-1, n-1] = sin*a[i-1, j-1] + cos*a[i-1, n-1]
            a[i-1, j-1] = temp
    for j in range(1, nm1+1):
        if abs(w[j-1]) > 1.0:
            cos = 1.0/w[j-1]
            sin = math.sqrt(1.0-cos*cos)
        if abs(w[j-1]) <= 1.0:
            sin = w[j-1]
            cos = math.sqrt(1.0-sin*sin)
        for i in range(1, m+1):
            temp = cos*a[i-1, j-1] + sin*a[i-1, n-1]
            a[i-1, n-1] = -sin*a[i-1, j-1] + cos*a[i-1, n-1]
            a[i-1, j-1] = temp

@njit(cache=True)
def hybrd(n, x, c4, t, kfe, cst):
    """MINPACK hybrd with mode=1, nprint=0, ml=mu=n-1 (SciPy defaults).
    Solves in place for x. Returns False if a function evaluation failed."""
    maxfev = 200*(n + 1)
    fvec = np.zeros(n)
    diag = np.zeros(n)
    fjac = np.zeros((n, n))
    r    = np.zeros((n*(n + 1))//2)
    qtf  = np.zeros(n)
    wa1  = np.zeros(n)
    wa2  = np.zeros(n)
    wa3  = np.zeros(n)
    wa4  = np.zeros(n)
    qtf2 = qtf.reshape((1, n))

    if not fcn(n, x, fvec, c4, t, kfe, cst):
        return False
    nfev = 1
    fnorm = enorm(n, fvec)
    msum = n
    it = 1
    ncsuc = 0
    ncfail = 0
    nslow1 = 0
    nslow2 = 0
    delta = 0.0
    xnorm = 0.0

    while True:
        jeval = True
        if not fdjac1(n, x, fvec, fjac, wa1, c4, t, kfe, cst):
            return False
        nfev = nfev + msum
        qrfac(n, n, fjac, wa1, wa2, wa3)
        if it == 1:
            for j in range(1, n+1):
                diag[j-1] = wa2[j-1]
                if wa2[j-1] == 0.0:
                    diag[j-1] = 1.0
            for j in range(1, n+1):
                wa3[j-1] = diag[j-1]*x[j-1]
            xnorm = enorm(n, wa3)
            delta = FACTOR*xnorm
            if delta == 0.0:
                delta = FACTOR
        for i in range(1, n+1):
            qtf[i-1] = fvec[i-1]
        for j in range(1, n+1):
            if fjac[j-1, j-1] != 0.0:
                s = 0.0
                for i in range(j, n+1):
                    s = s + fjac[i-1, j-1]*qtf[i-1]
                temp = -s/fjac[j-1, j-1]
                for i in range(j, n+1):
                    qtf[i-1] = qtf[i-1] + fjac[i-1, j-1]*temp
        for j in range(1, n+1):
            l = j
            for i in range(1, j):
                r[l-1] = fjac[i-1, j-1]
                l = l + n - i
            r[l-1] = wa1[j-1]
        qform(n, n, fjac, wa1)
        for j in range(1, n+1):
            diag[j-1] = max(diag[j-1], wa2[j-1])

        while True:
            dogleg(n, r, diag, qtf, delta, wa1, wa2, wa3)
            for j in range(1, n+1):
                wa1[j-1] = -wa1[j-1]
                wa2[j-1] = x[j-1] + wa1[j-1]
                wa3[j-1] = diag[j-1]*wa1[j-1]
            pnorm = enorm(n, wa3)
            if it == 1:
                delta = min(delta, pnorm)
            if not fcn(n, wa2, wa4, c4, t, kfe, cst):
                return False
            nfev = nfev + 1
            fnorm1 = enorm(n, wa4)
            actred = -1.0
            if fnorm1 < fnorm:
                q = fnorm1/fnorm
                actred = 1.0 - q*q
            l = 1
            for i in range(1, n+1):
                s = 0.0
                for j in range(i, n+1):
                    s = s + r[l-1]*wa1[j-1]
                    l = l + 1
                wa3[i-1] = qtf[i-1] + s
            temp = enorm(n, wa3)
            prered = 0.0
            if temp < fnorm:
                q = temp/fnorm
                prered = 1.0 - q*q
            ratio = 0.0
            if prered > 0.0:
                ratio = actred/prered
            if ratio < 0.1:
                ncsuc = 0
                ncfail = ncfail + 1
                delta = 0.5*delta
            else:
                ncfail = 0
                ncsuc = ncsuc + 1
                if ratio >= 0.5 or ncsuc > 1:
                    delta = max(delta, pnorm/0.5)
                if abs(ratio-1.0) <= 0.1:
                    delta = pnorm/0.5
            if ratio >= 1.0e-4:
                for j in range(1, n+1):
                    x[j-1] = wa2[j-1]
                    wa2[j-1] = diag[j-1]*x[j-1]
                    fvec[j-1] = wa4[j-1]
                xnorm = enorm(n, wa2)
                fnorm = fnorm1
                it = it + 1
            nslow1 = nslow1 + 1
            if actred >= 1.0e-3:
                nslow1 = 0
            if jeval:
                nslow2 = nslow2 + 1
            if actred >= 0.1:
                nslow2 = 0
            # --- tests for convergence / termination ---
            if delta <= XTOL*xnorm or fnorm == 0.0:
                return True
            if nfev >= maxfev:
                return True
            if 0.1*max(0.1*delta, pnorm) <= EPSMCH*xnorm:
                return True
            if nslow2 == 5 or nslow1 == 10:
                return True
            if ncfail == 2:
                break
            for j in range(1, n+1):
                s = 0.0
                for i in range(1, n+1):
                    s = s + fjac[i-1, j-1]*wa4[i-1]
                wa2[j-1] = (s - wa3[j-1])/pnorm
                wa1[j-1] = diag[j-1]*((diag[j-1]*wa1[j-1])/pnorm)
                if ratio >= 1.0e-4:
                    qtf[j-1] = s
            r1updt(n, n, r, wa1, wa2, wa3)
            r1mpyq(n, n, fjac, wa2, wa3)
            r1mpyq(1, n, qtf2, wa2, wa3)
            jeval = False

# ---------------------------------------------------------------------
# Beta equilibrium on the whole density grid
# ---------------------------------------------------------------------
@njit(cache=True)
def solve_beta(esym, rho, tx, kfex, cst):
    """Particle fractions at every density (same initial guesses as eos_v3.py).
    Returns (ok, ypr, yn, ye, ymu); ok is False if the EOS is rejected."""
    N = rho.shape[0]
    ypr = np.zeros(N)
    yn  = np.zeros(N)
    ye  = np.zeros(N)
    ymu = np.zeros(N)
    for i3 in range(N):
        c4 = 4.0*esym[i3]
        if rho[i3] <= 0.16:
            x = np.array([1.0e-2, 1.0, 1.0e-2])
            if not hybrd(3, x, c4, tx[i3], kfex[i3], cst):
                return False, ypr, yn, ye, ymu
            ypr[i3] = x[0]
            yn[i3]  = x[1]
            ye[i3]  = x[2]
            ymu[i3] = 0.0
        else:
            x = np.array([1.0e-2, 1.0, 1.0e-2, 1.0e-2])
            if not hybrd(4, x, c4, tx[i3], kfex[i3], cst):
                return False, ypr, yn, ye, ymu
            ypr[i3] = x[0]
            yn[i3]  = x[1]
            ye[i3]  = x[2]
            ymu[i3] = x[3]
    return True, ypr, yn, ye, ymu
