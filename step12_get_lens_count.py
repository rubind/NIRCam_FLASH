"""
microlensing_lmc_fast.py

Clean, fast pipeline:
  1) Fast analytic stage: tau, rate [1/s] using u_th(rho) lookup (no per-s bisection)
  2) Fast MC stage: efficiency epsilon (dimensionless) + (A_peak, x) pairs for detected hits
  3) Summary: expected detected events = Σ (rate_i * T_i * epsilon_i), plus uncertainties

Notes:
- Uses finite-source parameter rho = (R_* * x) / R_E  (projected source to lens plane).
- NFW normalized to local density rho_local_GeV_cm3 (GeV/c^2 per cm^3).
- Efficiency MC uses interval overlap with cadence, no lightcurve synthesis.
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any
import sys

# NumPy 2.4 removed the deprecated np.trapz alias.  Keep compatibility with
# both newer NumPy releases and older environments that lack np.trapezoid.
_trapezoid = getattr(np, "trapezoid", None)
if _trapezoid is None:
    _trapezoid = np.trapz

# ----------------------------
# Constants
# ----------------------------
G = 6.67430e-11
c = 299792458.0
Msun = 1.98847e30
kpc = 3.085677581e19

# LMC line-of-sight (approx)
LMC_l = np.deg2rad(280.0)
LMC_b = np.deg2rad(-33.0)

R0 = 8.2 * kpc
DS_DEFAULT = 50.0 * kpc

# ----------------------------
# Data model
# ----------------------------
@dataclass
class Star:
    m_short: float
    m_long: float
    R_star_m: float
    times_s: np.ndarray
    name: str = "star"
    DS_m: float = DS_DEFAULT

@dataclass
class MCConfig:
    n_events_try: int = 200_000
    sigma_v: float = 160e3
    Ath: float = 1.02

# ----------------------------
# Halo model (NFW)
# ----------------------------
def gev_cm3_to_kg_m3(rho_GeV_cm3: float) -> float:
    """
    Convert GeV/c^2 per cm^3 -> kg/m^3.
    1 GeV/c^2 = 1.78266192e-27 kg
    1 cm^-3 = 1e6 m^-3
    """
    return rho_GeV_cm3 * 1.78266192e-27 * 1e6

@dataclass
class NFWModel:
    r_s_m: float = 16.0 * kpc
    rho_local_GeV_cm3: float = 0.4
    R0_m: float = R0
    l_rad: float = LMC_l
    b_rad: float = LMC_b

    def __post_init__(self):
        self.rho_local = gev_cm3_to_kg_m3(self.rho_local_GeV_cm3)
        x = self.R0_m / self.r_s_m
        # Choose rho_s so that rho_NFW(R0)=rho_local
        self.rho_s = self.rho_local * x * (1 + x) ** 2

    def rho(self, r_m: np.ndarray) -> np.ndarray:
        x = r_m / self.r_s_m
        return self.rho_s / (x * (1 + x) ** 2)

    def r_galactocentric(self, s_m: np.ndarray) -> np.ndarray:
        # r^2 = R0^2 + s^2 - 2 R0 s cos(b) cos(l)
        return np.sqrt(self.R0_m**2 + s_m**2 - 2 * self.R0_m * s_m * np.cos(self.b_rad) * np.cos(self.l_rad))

    def density_along_los(self, s_m: np.ndarray) -> np.ndarray:
        return self.rho(self.r_galactocentric(s_m))

@dataclass
class M31NFWModel:
    """M31-centered NFW halo along a line of sight through projected offset b_proj_m."""
    r_s_m: float = 25.0 * kpc
    rho_s_Msun_kpc3: float = 4.96e6
    D_center_m: float = 780.0 * kpc
    b_proj_m: float = 0.0

    def __post_init__(self):
        self.rho_s = self.rho_s_Msun_kpc3 * Msun / kpc**3

    def rho(self, r_m: np.ndarray) -> np.ndarray:
        # The NFW cusp is integrable in the rate and optical-depth integrals.
        r_safe = np.maximum(np.asarray(r_m, dtype=float), 1e-12 * self.r_s_m)
        x = r_safe / self.r_s_m
        return self.rho_s / (x * (1 + x) ** 2)

    def density_along_los(self, s_m: np.ndarray) -> np.ndarray:
        r = np.sqrt((self.D_center_m - np.asarray(s_m))**2 + self.b_proj_m**2)
        return self.rho(r)

@dataclass
class CompositeHalo:
    """Sum independent lens populations while retaining a common f_pbh."""
    components: Tuple[Any, ...]

    def density_along_los(self, s_m: np.ndarray) -> np.ndarray:
        return np.sum([component.density_along_los(s_m) for component in self.components], axis=0)

# ----------------------------
# Microlensing geometry
# ----------------------------
def R_E(M_kg: float, s_m: np.ndarray, DS_m: float) -> np.ndarray:
    x = s_m / DS_m
    return np.sqrt((4 * G * M_kg / c**2) * DS_m * x * (1 - x))

def A_point(u: np.ndarray) -> np.ndarray:
    u = np.maximum(u, 1e-12)
    return (u*u + 2) / (u * np.sqrt(u*u + 4))

_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(96)

def A_finite(u: float, rho: float) -> float:
    """
    Uniform-disk finite-source magnification using a nonsingular 1-D overlap integral.

    Polar coordinates are centered on the lens.  The factor r*A_point(r) is
    finite at r=0, avoiding the point-lens singularity that biases a direct
    two-dimensional source grid.
    """
    u = abs(float(u))
    rho = float(rho)
    if rho <= 0:
        return float(A_point(np.array([u]))[0])

    # Exact central value for a uniform disk.
    if u == 0.0:
        return float(np.sqrt(rho * rho + 4.0) / rho)

    # Finite-source corrections are negligible here and the overlap integral
    # loses relative precision as its width approaches machine precision.
    if rho < 1e-4:
        return float(A_point(np.array([u]))[0])

    # Full annuli lie inside the source for r < rho-u when u<rho.
    total = 0.0
    if rho > u:
        r_full = rho - u
        total += np.pi * r_full * np.sqrt(r_full * r_full + 4.0)

    # Partial annuli span |rho-u| < r < rho+u.  Gauss-Legendre nodes avoid
    # both endpoints, including the integrable contact singularities.
    r_lo = abs(rho - u)
    r_hi = rho + u
    half = 0.5 * (r_hi - r_lo)
    mid = 0.5 * (r_hi + r_lo)
    r = mid + half * _GL_NODES
    z = (r*r + u*u - rho*rho) / (2.0 * r * u)
    angle = 2.0 * np.arccos(np.clip(z, -1.0, 1.0))
    rA = (r*r + 2.0) / np.sqrt(r*r + 4.0)
    total += half * float(np.sum(_GL_WEIGHTS * angle * rA))
    return total / (np.pi * rho * rho)

def u_thresh_finite(Ath: float, rho: float) -> float:
    """
    Find u where A_finite(u,rho)=Ath by bisection; returns 0 if Amax<Ath.
    """
    Amax = A_finite(0.0, rho)
    if Amax <= Ath:
        return 0.0

    lo, hi = 0.0, max(1.0, rho + 1.0)
    # Ensure bracket
    while A_finite(hi, rho) > Ath and hi < 1e6:
        hi *= 2.0

    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if A_finite(mid, rho) > Ath:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

def build_uth_table(Ath: float, rho_min: float = 1e-4, rho_max: float = 1e4, n_rho: int = 512) -> Tuple[np.ndarray, np.ndarray]:
    if Ath <= 1.0:
        raise ValueError("Ath must be greater than 1.")
    rho_cut = 2.0 / np.sqrt(Ath * Ath - 1.0)
    rho_grid = np.logspace(np.log10(rho_min), np.log10(rho_max), n_rho)
    if rho_min < rho_cut < rho_max:
        rho_grid = np.unique(np.sort(np.append(rho_grid, rho_cut)))
    u_grid = np.zeros_like(rho_grid)
    for i, rho in enumerate(rho_grid):
        if rho < rho_cut:
            u_grid[i] = u_thresh_finite(Ath, float(rho))
    return rho_grid, u_grid

def interp_uth(rho_vals: np.ndarray, rho_grid: np.ndarray, u_grid: np.ndarray) -> np.ndarray:
    # Log interpolation in rho
    rho_clipped = np.clip(rho_vals, rho_grid[0], rho_grid[-1])
    return np.interp(np.log10(rho_clipped), np.log10(rho_grid), u_grid)

# ----------------------------
# Fast analytic stage
# ----------------------------
def mean_v_perp(sigma_v: float) -> float:
    """Mean speed of a Rayleigh distribution made from two N(0,sigma_v) components."""
    return float(np.sqrt(np.pi / 2.0) * sigma_v)

def analytic_tau_and_rate_fast(
    star: Star,
    M_kg: float,
    halo: Any,
    Ath: float,
    f_pbh: float = 1.0,
    sigma_v: float = 160e3,
    nsamp_s: int = 600,
    uth_table: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> Tuple[float, float]:
    """
    Returns (tau, rate) where:
      tau  = optical depth (dimensionless)
      rate = event rate [1/s] for amplification >= Ath (thresholded tube)
    """
    if uth_table is None:
        uth_table = build_uth_table(Ath=Ath)
    rho_grid, u_grid = uth_table

    DS_m = star.DS_m
    # Cosine spacing resolves both the nearby MW halo and an integrable M31
    # cusp near the source much better than a uniform distance grid.
    z = np.linspace(1e-5, 1.0 - 1e-5, nsamp_s)
    x = 0.5 * (1.0 - np.cos(np.pi * z))
    s = DS_m * x
    RE = R_E(M_kg, s, DS_m)

    # Finite-source parameter with projection to lens plane
    rho_star = (star.R_star_m * x) / RE
    u_th = interp_uth(rho_star, rho_grid, u_grid)

    sigma = np.pi * (u_th * RE) ** 2
    n_lens = f_pbh * halo.density_along_los(s) / M_kg
    vbar = mean_v_perp(sigma_v)

    tau = float(_trapezoid(n_lens * sigma, x=s))
    rate = float(_trapezoid(n_lens * (2.0 * u_th * RE) * vbar, x=s))
    return tau, rate

def run_analytic_fast(
    stars: List[Star],
    M_kg: float,
    halo: Any,
    Ath: float = 1.02,
    f_pbh: float = 1.0,
    sigma_v: float = 160e3,
    nsamp_s: int = 600,
    uth_table: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Returns dict name -> {'tau':..., 'rate':...}
    """
    if uth_table is None:
        uth_table = build_uth_table(Ath=Ath)
    out: Dict[str, Dict[str, float]] = {}
    for s in stars:
        tau, rate = analytic_tau_and_rate_fast(
            s, M_kg, halo, Ath, f_pbh=f_pbh, sigma_v=sigma_v, nsamp_s=nsamp_s, uth_table=uth_table
        )
        out[s.name] = {"tau": tau, "rate": rate}
    return out

# ----------------------------
# Fast MC stage (efficiency only)
# ----------------------------
def make_s_sampler(
    halo: Any,
    star: Star,
    M_kg: float,
    uth_table: Tuple[np.ndarray, np.ndarray],
    n_grid: int = 4096,
) -> Tuple[np.ndarray, np.ndarray]:
    # q(s) is proportional to the marginal spatial event rate
    # rho(s)*u_th(s)*sqrt(x(1-x)).  It is normalizable at an M31 NFW cusp
    # and avoids wasting trials in finite-source regions with u_th=0.
    DS_m = star.DS_m
    z = np.linspace(1e-5, 1.0 - 1e-5, n_grid)
    x = 0.5 * (1.0 - np.cos(np.pi * z))
    s = DS_m * x
    geom = np.sqrt(x * (1.0 - x))
    RE = R_E(M_kg, s, DS_m)
    rho_star = star.R_star_m * x / RE
    u_th = interp_uth(rho_star, uth_table[0], uth_table[1])
    w = halo.density_along_los(s) * geom * u_th
    cdf = np.zeros_like(s)
    cdf[1:] = np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(s))
    if cdf[-1] > 0.0:
        cdf /= cdf[-1]
    else:
        cdf = np.linspace(0.0, 1.0, n_grid)
    return s, cdf

def sample_s_from_cdf(n: int, s_grid: np.ndarray, cdf_grid: np.ndarray) -> np.ndarray:
    u = np.random.random(n)
    return np.interp(u, cdf_grid, s_grid)

def draw_v_perp(n: int, sigma_v: float) -> np.ndarray:
    return np.sqrt(np.random.normal(0, sigma_v, n) ** 2 + np.random.normal(0, sigma_v, n) ** 2)

def any_hit_interval(t_sorted: np.ndarray, t0s: np.ndarray, dt_half: np.ndarray) -> np.ndarray:
    left = np.searchsorted(t_sorted, t0s - dt_half, side="left")
    right = np.searchsorted(t_sorted, t0s + dt_half, side="right")
    return (right - left) > 0

def simulate_star_mc_fast_eff(
    star: Star,
    M_kg: float,
    halo: Any,
    cfg: MCConfig,
    s_sampler: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    uth_table: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    record_pairs: bool = True,
) -> Dict[str, Any]:
    """
    Returns efficiency only (dimensionless), plus optional (A_peak, x) for detected hits.
    """
    t_sorted = np.asarray(star.times_s, dtype=float)
    if t_sorted.ndim != 1 or len(t_sorted) < 2:
        raise ValueError(f"{star.name}: times_s must be a 1D array with >=2 entries.")
    t_sorted = np.sort(t_sorted)
    T_obs = float(t_sorted[-1] - t_sorted[0])

    if uth_table is None:
        uth_table = build_uth_table(Ath=cfg.Ath)
    rho_grid, u_grid = uth_table

    if s_sampler is None:
        s_sampler = make_s_sampler(halo, star, M_kg, uth_table)
    s_grid, cdf = s_sampler

    N = int(cfg.n_events_try)
    s = sample_s_from_cdf(N, s_grid, cdf)
    x = s / star.DS_m

    RE = R_E(M_kg, s, star.DS_m)
    rho = (star.R_star_m * x) / RE
    u_th = interp_uth(rho, rho_grid, u_grid)

    feasible = u_th > 0.0
    n_feasible = int(np.count_nonzero(feasible))
    if n_feasible == 0:
        return {
            "name": star.name,
            "epsilon": 0.0,
            "sigma_epsilon": 0.0,
            "n_trials": N,
            "n_feasible": 0,
            "n_hits": 0,
            "pairs_Apeak_x": [],
        }

    vperp = draw_v_perp(N, cfg.sigma_v)
    tE = RE / np.maximum(vperp, 1.0)

    # Event trajectories crossing the lensing tube have u0 uniform on [0, u_th].
    # The sqrt draw (p(u0) proportional to u0) would instead describe optical depth.
    u0 = np.zeros(N)
    u0_vals = u_th[feasible] * np.random.random(n_feasible)
    u0[feasible] = u0_vals

    # Duration above threshold
    dt = np.zeros(N)
    dt[feasible] = tE[feasible] * np.sqrt(np.maximum(u_th[feasible] ** 2 - u0[feasible] ** 2, 0.0))

    # The analytic rate counts event peaks per unit time.  Draw peaks over the
    # same observing baseline used later in lambda = rate * T_obs * epsilon.
    t0 = np.random.uniform(t_sorted[0], t_sorted[-1], size=N)

    hits = any_hit_interval(t_sorted, t0, dt) & feasible

    # NEW: enforce true peak magnification threshold
    if np.any(hits):
        A_peak = np.array([A_finite(float(u), float(r)) for u, r in zip(u0[hits], rho[hits])], dtype=float)
        keep = A_peak >= cfg.Ath
        # shrink hits to those truly above threshold
        hit_idx = np.where(hits)[0]
        hits[:] = False
        hits[hit_idx[keep]] = True

    n_hits = int(np.count_nonzero(hits))

    # q(s) already contains rho(s)*u_th(s)*sqrt(x(1-x)); since R_E is
    # proportional to sqrt(x(1-x)), only the event-rate speed factor remains.
    event_weight = np.where(feasible, vperp, 0.0)
    weight_sum = float(np.sum(event_weight))
    hit_weight = float(np.sum(event_weight[hits]))
    epsilon = hit_weight / weight_sum if weight_sum > 0.0 else 0.0

    # Self-normalized importance-sampling variance.  At epsilon=0 or 1 the
    # empirical variance collapses, so use a Jeffreys effective-count fallback
    # rather than incorrectly reporting zero Monte Carlo uncertainty.
    weight_sq_sum = float(np.sum(event_weight**2))
    n_effective = weight_sum**2 / weight_sq_sum if weight_sq_sum > 0.0 else 0.0
    if 0.0 < epsilon < 1.0:
        centered = hits.astype(float) - epsilon
        var_epsilon = float(np.sum((event_weight * centered) ** 2) / weight_sum**2)
        if n_effective > 1.0:
            var_epsilon *= n_effective / (n_effective - 1.0)
        sigma_epsilon = float(np.sqrt(max(var_epsilon, 0.0)))
    elif n_effective > 0.0:
        alpha = n_effective * epsilon + 0.5
        beta = n_effective * (1.0 - epsilon) + 0.5
        sigma_epsilon = float(np.sqrt(alpha * beta / ((alpha + beta)**2 * (alpha + beta + 1.0))))
    else:
        sigma_epsilon = 0.0

    pairs: List[Tuple[float, float]] = []
    if record_pairs and n_hits > 0:
        # A_peak at closest approach is A_finite(u0, rho)
        # This is only evaluated for hits, so usually cheap.
        # Resample detected pairs by their event-rate weights so downstream
        # histograms represent the same event population as epsilon.
        hit_idx = np.where(hits)[0]
        hit_prob = event_weight[hit_idx] / np.sum(event_weight[hit_idx])
        pair_idx = np.random.choice(hit_idx, size=n_hits, replace=True, p=hit_prob)
        A_peak = np.array([A_finite(float(u), float(r)) for u, r in zip(u0[pair_idx], rho[pair_idx])], dtype=float)
        pairs = list(zip(A_peak.tolist(), x[pair_idx].tolist()))

    return {
        "name": star.name,
        "epsilon": float(epsilon),
        "sigma_epsilon": float(sigma_epsilon),
        "n_trials": N,
        "n_feasible": n_feasible,
        "n_hits": n_hits,
        "n_effective": float(n_effective),
        "pairs_Apeak_x": pairs,
    }

def run_mc_fast_eff(
    stars: List[Star],
    M_kg: float,
    halo: Any,
    Ath: float = 1.02,
    n_events_try: int = 200_000,
    sigma_v: float = 160e3,
    uth_table: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    record_pairs: bool = True,
) -> Tuple[List[Tuple], List[Tuple[str, float, float]]]:
    """
    Returns:
      rows_eff: [(name, epsilon, sigma_epsilon, n_trials, n_feasible, n_hits), ...]
      pairs:    [(name, A_peak, x), ...]  for detected hits only
    """
    cfg = MCConfig(n_events_try=n_events_try, sigma_v=sigma_v, Ath=Ath)
    if uth_table is None:
        uth_table = build_uth_table(Ath=Ath)

    rows: List[Tuple] = []
    all_pairs: List[Tuple[str, float, float]] = []
    for s in stars:
        sampler = make_s_sampler(halo, s, M_kg, uth_table)
        res = simulate_star_mc_fast_eff(
            s, M_kg, halo, cfg,
            s_sampler=sampler,
            uth_table=uth_table,
            record_pairs=record_pairs,
        )
        rows.append((res["name"], res["epsilon"], res["sigma_epsilon"], res["n_trials"], res["n_feasible"], res["n_hits"]))
        if record_pairs and res["pairs_Apeak_x"]:
            all_pairs.extend([(s.name, A, x) for (A, x) in res["pairs_Apeak_x"]])
    return rows, all_pairs

# ----------------------------
# Summary stage: combine analytic rate + MC efficiency + observing time
# ----------------------------
def summarize_run_with_rates(
    stars: List[Star],
    rows_eff: List[Tuple],                 # (name, epsilon, sigma_epsilon, n_trials, n_feasible, n_hits)
    rates_by_name: Dict[str, Dict[str, float]],  # from run_analytic_fast or your own cache
) -> Dict[str, Any]:
    """
    Computes:
      lambda_i = rate_i * T_i * epsilon_i
      N_exp = Σ lambda_i
      P(≥1) = 1 - exp(-Σ lambda_i)

    Uncertainty includes MC efficiency uncertainty only (and whatever sigma_rate you include in rates_by_name).
    If you want, include sigma_rate in rates_by_name[name]['sigma_rate'].
    """
    eff = {r[0]: {"epsilon": float(r[1]), "sigma_epsilon": float(r[2]), "n_trials": int(r[3]), "n_feasible": int(r[4]), "n_hits": int(r[5])}
           for r in rows_eff}

    per_star = []
    lambdas = []
    sigma_lambdas = []
    T_total = 0.0

    for s in stars:
        if s.name not in eff:
            raise KeyError(f"No MC efficiency for star '{s.name}'.")
        if s.name not in rates_by_name:
            raise KeyError(f"No analytic rate for star '{s.name}' in rates_by_name.")

        T_obs = float(np.max(s.times_s) - np.min(s.times_s))
        T_total += T_obs

        epsilon = eff[s.name]["epsilon"]
        sigma_epsilon = eff[s.name]["sigma_epsilon"]

        rate = float(rates_by_name[s.name].get("rate", 0.0))
        tau = float(rates_by_name[s.name].get("tau", 0.0))
        sigma_rate = float(rates_by_name[s.name].get("sigma_rate", 0.0))

        lam = rate * T_obs * epsilon
        sigma_lam = float(np.sqrt((T_obs * epsilon * sigma_rate) ** 2 + (T_obs * rate * sigma_epsilon) ** 2))

        p = 1.0 - np.exp(-lam)
        sigma_p = (1.0 - p) * sigma_lam

        per_star.append({
            "name": s.name,
            "T_obs": T_obs,
            "tau": tau,
            "rate": rate,
            "sigma_rate": sigma_rate,
            "epsilon": epsilon,
            "sigma_epsilon": sigma_epsilon,
            "lambda_val": lam,
            "sigma_lambda": sigma_lam,
            "p": p,
            "sigma_p": sigma_p,
            "n_trials": eff[s.name]["n_trials"],
            "n_feasible": eff[s.name]["n_feasible"],
            "n_hits": eff[s.name]["n_hits"],
        })

        lambdas.append(lam)
        sigma_lambdas.append(sigma_lam)

    Lambda = float(np.sum(lambdas))
    sigma_Lambda = float(np.sqrt(np.sum(np.array(sigma_lambdas) ** 2)))
    N_exp = Lambda
    sigma_N_exp = sigma_Lambda

    P_ge1 = 1.0 - np.exp(-Lambda)
    sigma_P_ge1 = (1.0 - P_ge1) * sigma_Lambda

    return {
        "total_observing_time_s": float(T_total),
        "Lambda": Lambda,
        "sigma_Lambda": sigma_Lambda,
        "expected_events": N_exp,
        "sigma_expected_events": sigma_N_exp,
        "P_at_least_one": P_ge1,
        "sigma_P_at_least_one": sigma_P_ge1,
        "per_star": per_star,
    }

def run_population_breakdown(
    stars: List[Star],
    M_kg: float,
    populations: Dict[str, Tuple[Any, float]],
    Ath: float = 1.02,
    f_pbh: float = 1.0,
    n_events_try: int = 200_000,
    nsamp_s: int = 1200,
    uth_table: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    record_pairs: bool = True,
) -> Dict[str, Any]:
    """Run and retain separate analytic/MC results for independent lens populations.

    populations maps a label to (halo_model, Rayleigh component sigma_v).
    Poisson means add across populations; efficiencies must not be averaged.
    """
    if uth_table is None:
        uth_table = build_uth_table(Ath=Ath)

    by_population: Dict[str, Dict[str, Any]] = {}
    tagged_pairs: List[Tuple[str, str, float, float]] = []
    for label, (halo, sigma_v) in populations.items():
        rates = run_analytic_fast(stars, M_kg, halo, Ath=Ath, f_pbh=f_pbh,
                                  sigma_v=sigma_v, nsamp_s=nsamp_s, uth_table=uth_table)
        rows_eff, pairs = run_mc_fast_eff(stars, M_kg, halo, Ath=Ath,
                                          n_events_try=n_events_try, sigma_v=sigma_v,
                                          uth_table=uth_table, record_pairs=record_pairs)
        summary = summarize_run_with_rates(stars, rows_eff, rates)
        by_population[label] = {
            "halo": halo,
            "sigma_v": sigma_v,
            "rates": rates,
            "rows_eff": rows_eff,
            "pairs_Apeak_x": pairs,
            "summary": summary,
        }
        tagged_pairs.extend((label, star_name, A_peak, x) for star_name, A_peak, x in pairs)

    Lambda = float(sum(item["summary"]["Lambda"] for item in by_population.values()))
    sigma_Lambda = float(np.sqrt(sum(item["summary"]["sigma_Lambda"]**2 for item in by_population.values())))
    P_ge1 = 1.0 - np.exp(-Lambda)
    return {
        "populations": by_population,
        "Lambda": Lambda,
        "sigma_Lambda": sigma_Lambda,
        "expected_events": Lambda,
        "sigma_expected_events": sigma_Lambda,
        "P_at_least_one": P_ge1,
        "sigma_P_at_least_one": (1.0 - P_ge1) * sigma_Lambda,
        "pairs_population_star_Apeak_x": tagged_pairs,
    }

# ----------------------------
# Histogram helper (optional)
# ----------------------------
def histogram_amplifications(
    pairs: List[Tuple[str, float, float]],
    bins: int = 40,
    range_deltaA: Optional[Tuple[float, float]] = None,
    density: bool = False,
):
    """
    pairs: [(name, A_peak, x), ...]
    Returns (hist, edges, deltaA)
    """
    if len(pairs) == 0:
        return np.array([]), np.array([]), np.array([])

    A_pk = np.array([p[1] for p in pairs], dtype=float)
    dA = A_pk - 1.0

    if range_deltaA is None:
        lo = float(np.percentile(dA, 0.5))
        hi = float(np.percentile(dA, 99.5))
        if lo == hi:
            lo, hi = float(np.min(dA)), float(np.max(dA) + 1e-6)
        range_deltaA = (lo, hi)

    hist, edges = np.histogram(dA, bins=bins, range=range_deltaA, density=density)
    return hist, edges, dA



# ----------------------------
# Example usage (copy into your notebook, not required to run here)
# ----------------------------
M31_l = np.deg2rad(121.17)
M31_b = np.deg2rad(-21.57)
DS_M31 = 780.0 * kpc

if __name__ == "__main__":
    # Example: 3 stars, 400 hours, 20s cadence (as you used)

    
    star_radius_rsol = float(sys.argv[1])
    star_hours = float(sys.argv[2])
    one_sigma_phot = float(sys.argv[3])
    BH_mass = float(sys.argv[4])
  
    times_s = np.arange(0, star_hours * 3600, 10.737*2)

    #stars = [
    #    Star(m_short=22.1, m_long=21.9, R_star_m=star_radius_rsol * 6.957e8, times_s=times_s, name="LMC-1"),
    #]
    stars_m31 = [Star(m_short=22.1, m_long=21.9, R_star_m=star_radius_rsol  * 6.957e8, times_s=times_s, name="M31-1", DS_m=DS_M31)]

#Andromeda


#t = np.arange(0, 0.63*7*6.4*1e6 * 3600, 270.0)
#t = np.arange(0, 3600, 270.0)

halo_mw_to_m31 = NFWModel(r_s_m=16*kpc, rho_local_GeV_cm3=0.4, l_rad=M31_l, b_rad=M31_b)
halo_m31 = M31NFWModel(D_center_m=DS_M31)
Ath = 1.0 + one_sigma_phot*4
M = BH_mass * Msun
f_pbh = 1.0
sigma_v_mw = 160e3
# LensCalcPy uses exp(-v^2/v_disp^2) with v_disp=250 km/s for M31,
# equivalent to Rayleigh component sigma=v_disp/sqrt(2).
sigma_v_m31 = 250e3 / np.sqrt(2.0)


uth_table = build_uth_table(Ath=Ath)

# Separate Monte Carlo runs are required because the two populations have
# different distance and velocity distributions.
population_result = run_population_breakdown(
    stars_m31, M,
    populations={
        "Milky Way": (halo_mw_to_m31, sigma_v_mw),
        "M31": (halo_m31, sigma_v_m31),
    },
    Ath=Ath, f_pbh=f_pbh, n_events_try=200_000, uth_table=uth_table,
)

for label, result in population_result["populations"].items():
    row = result["rows_eff"][0]
    print(f"{label:10s}: rate={result['rates']['M31-1']['rate']*3600:.6g}/hr, "
          f"epsilon={row[1]:.5g} ± {row[2]:.2g}, "
          f"N_exp={result['summary']['expected_events']:.6g}")
print("combined:", population_result["Lambda"], "±", population_result["sigma_Lambda"],
      " P(>=1)=", population_result["P_at_least_one"])
