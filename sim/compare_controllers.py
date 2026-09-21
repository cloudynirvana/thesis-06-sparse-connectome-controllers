#!/usr/bin/env python3
"""
Closed-loop comparison of lumped adaptive-therapy (AT) vs sparse
Kenyon-cell-style policies on a toy two-clone cancer ODE.

Computational research only. Not a medical device, not a protocol,
not a dose, and not a claim about fly neurons as treatment.

Plant is a competitive logistic pair (sensitive S, resistant R) with a
scalar input u in [0, 1]. Controllers are architecture classes, not
therapies.

Usage:
    python3 sim/compare_controllers.py
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

RNG_SEED = 20260921
DT = 0.05
T_FINAL = 200.0
K_CAPACITY = 1.0


@dataclass(frozen=True)
class PlantParams:
    r_s: float = 0.28
    r_r: float = 0.16
    a_sr: float = 1.0
    a_rs: float = 1.6
    d_s: float = 0.55
    d_r: float = 0.06
    k: float = K_CAPACITY


@dataclass
class ATParams:
    t_on: float = 0.50
    t_off: float = 0.25
    u_on: float = 1.0


@dataclass
class AffineParams:
    w_s: float = 2.4
    w_r: float = -1.8
    b: float = -0.55


@dataclass
class KCParams:
    n_pn: int = 8
    n_kc: int = 96
    n_claw: int = 3
    top_k: int = 8
    ridge: float = 1e-2


def rhs(state: np.ndarray, u: float, p: PlantParams) -> np.ndarray:
    s, r = state
    s = max(s, 0.0)
    r = max(r, 0.0)
    load = (s + p.a_sr * r) / p.k
    load_r = (r + p.a_rs * s) / p.k
    ds = p.r_s * s * (1.0 - load) - p.d_s * u * s
    dr = p.r_r * r * (1.0 - load_r) - p.d_r * u * r
    return np.array([ds, dr], dtype=float)


def rk4_step(state: np.ndarray, u: float, p: PlantParams, dt: float) -> np.ndarray:
    k1 = rhs(state, u, p)
    k2 = rhs(state + 0.5 * dt * k1, u, p)
    k3 = rhs(state + 0.5 * dt * k2, u, p)
    k4 = rhs(state + dt * k3, u, p)
    nxt = state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return np.maximum(nxt, 0.0)


def total_burden(state: np.ndarray) -> float:
    return float(state[0] + state[1])


def pn_features(state: np.ndarray, u_prev: float) -> np.ndarray:
    s, r = state
    t = s + r
    frac_s = s / t if t > 1e-8 else 0.0
    frac_r = r / t if t > 1e-8 else 0.0
    return np.array(
        [s, r, t, 1.0 - min(t, 1.0), frac_s, frac_r, u_prev, 1.0],
        dtype=float,
    )


class LumpedAT:
    name = "lumped_at"

    def __init__(self, params: ATParams):
        self.p = params
        self.u = 0.0

    def reset(self) -> None:
        self.u = 0.0

    def act(self, state: np.ndarray) -> float:
        t = total_burden(state)
        if t >= self.p.t_on:
            self.u = self.p.u_on
        elif t <= self.p.t_off:
            self.u = 0.0
        return self.u


class AffineComposition:
    name = "affine_sr"

    def __init__(self, params: AffineParams):
        self.p = params

    def reset(self) -> None:
        return

    def act(self, state: np.ndarray) -> float:
        s, r = state
        raw = self.p.w_s * s + self.p.w_r * r + self.p.b
        return float(np.clip(raw, 0.0, 1.0))


class SparseKC:
    name = "sparse_kc"

    def __init__(self, params: KCParams, rng: np.random.Generator):
        self.p = params
        self.u_prev = 0.0
        claws = np.array(
            [rng.choice(params.n_pn, size=params.n_claw, replace=False) for _ in range(params.n_kc)]
        )
        self.claws = claws
        self.weights = rng.normal(0.0, 1.0, size=(params.n_kc, params.n_claw))
        self.bias = rng.normal(0.0, 0.25, size=params.n_kc)
        self.readout = np.zeros(params.n_kc, dtype=float)
        self.readout_bias = 0.0

    def reset(self) -> None:
        self.u_prev = 0.0

    def code(self, state: np.ndarray) -> np.ndarray:
        pn = pn_features(state, self.u_prev)
        drive = np.einsum("ij,ij->i", self.weights, pn[self.claws]) + self.bias
        k = self.p.top_k
        idx = np.argpartition(drive, -k)[-k:]
        h = np.zeros_like(drive)
        h[idx] = np.maximum(drive[idx], 0.0)
        nrm = np.linalg.norm(h)
        if nrm > 1e-12:
            h = h / nrm
        return h

    def act(self, state: np.ndarray) -> float:
        h = self.code(state)
        u = float(np.clip(h @ self.readout + self.readout_bias, 0.0, 1.0))
        self.u_prev = u
        return u

    def fit_readout(self, states: np.ndarray, targets: np.ndarray) -> None:
        codes = np.stack([self.code(s) for s in states], axis=0)
        x = np.concatenate([codes, np.ones((len(states), 1))], axis=1)
        xtx = x.T @ x + self.p.ridge * np.eye(x.shape[1])
        w = np.linalg.solve(xtx, x.T @ targets)
        self.readout = w[:-1]
        self.readout_bias = float(w[-1])
        self.reset()


def gatenby_like_target(state: np.ndarray) -> float:
    """Composition-aware treat-for-stability target used only to train the KC readout.

    Treat when burden is high *and* a competitive sensitive majority remains.
    This is an in-silico training target, not a clinical rule.
    """
    s, r = state
    t = s + r
    frac_s = s / t if t > 1e-8 else 0.0
    if t >= 0.50 and frac_s >= 0.35:
        return 1.0
    if t <= 0.22:
        return 0.0
    if frac_s < 0.20:
        return 0.0
    return 0.45 if t > 0.35 else 0.0


def simulate(controller, state0: np.ndarray, plant: PlantParams, t_final: float = T_FINAL) -> dict:
    controller.reset()
    n = int(round(t_final / DT))
    s = np.array(state0, dtype=float)
    t_hist = np.zeros(n + 1)
    s_hist = np.zeros(n + 1)
    r_hist = np.zeros(n + 1)
    u_hist = np.zeros(n + 1)
    t_hist[0] = total_burden(s)
    s_hist[0], r_hist[0] = s
    u_hist[0] = controller.act(s)
    for i in range(n):
        u = u_hist[i]
        s = rk4_step(s, u, plant, DT)
        t_hist[i + 1] = total_burden(s)
        s_hist[i + 1], r_hist[i + 1] = s
        u_hist[i + 1] = controller.act(s)
    switches = int(np.sum(np.abs(np.diff(u_hist > 0.5)) > 0))
    duty = float(np.mean(u_hist))
    t_mean = float(np.mean(t_hist))
    t_max = float(np.max(t_hist))
    frac_r_end = float(r_hist[-1] / max(t_hist[-1], 1e-12))
    time_r_dom = None
    for i, (ss, rr) in enumerate(zip(s_hist, r_hist)):
        if rr > ss and (ss + rr) > 0.05:
            time_r_dom = i * DT
            break
    return {
        "t": t_hist,
        "s": s_hist,
        "r": r_hist,
        "u": u_hist,
        "switches": switches,
        "duty": duty,
        "t_mean": t_mean,
        "t_max": t_max,
        "frac_r_end": frac_r_end,
        "time_r_dom": time_r_dom,
        "u0": float(u_hist[0]),
    }


def isoline_u_variance(controller, t_value: float, n: int = 81) -> dict:
    fracs = np.linspace(0.05, 0.95, n)
    us = []
    for f in fracs:
        s = t_value * f
        r = t_value * (1.0 - f)
        controller.reset()
        us.append(controller.act(np.array([s, r])))
    us = np.array(us, dtype=float)
    return {
        "mean": float(np.mean(us)),
        "std": float(np.std(us)),
        "min": float(np.min(us)),
        "max": float(np.max(us)),
        "range": float(np.max(us) - np.min(us)),
    }


def retune_at(traj_t: np.ndarray, grid_on=None, grid_off=None) -> ATParams:
    """Fit hysteresis thresholds to match a binary ON/OFF trace of another policy."""
    if grid_on is None:
        grid_on = np.linspace(0.30, 0.80, 26)
    if grid_off is None:
        grid_off = np.linspace(0.05, 0.45, 21)
    target = (traj_t > np.median(traj_t)).astype(float)
    # Use the *burden* series as the only AT observation; match duty cycle
    # of the reference u via closed-loop later. Here we just scan thresholds
    # against the reference burden path as a 1-D reduction.
    best = None
    best_err = np.inf
    for t_on in grid_on:
        for t_off in grid_off:
            if t_off >= t_on:
                continue
            u = 0.0
            pred = np.zeros_like(traj_t)
            for i, t in enumerate(traj_t):
                if t >= t_on:
                    u = 1.0
                elif t <= t_off:
                    u = 0.0
                pred[i] = u
            err = float(np.mean((pred - target) ** 2))
            if err < best_err:
                best_err = err
                best = ATParams(t_on=float(t_on), t_off=float(t_off))
    assert best is not None
    return best


def decision_region_complexity(controller, n: int = 41) -> dict:
    grid = np.linspace(0.02, 0.90, n)
    field = np.zeros((n, n))
    for i, s in enumerate(grid):
        for j, r in enumerate(grid):
            controller.reset()
            field[i, j] = 1.0 if controller.act(np.array([s, r])) > 0.5 else 0.0
    # 4-connected component count of the ON set
    seen = np.zeros_like(field, dtype=bool)
    comps = 0
    for i in range(n):
        for j in range(n):
            if field[i, j] < 0.5 or seen[i, j]:
                continue
            comps += 1
            stack = [(i, j)]
            seen[i, j] = True
            while stack:
                a, b = stack.pop()
                for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    na, nb = a + da, b + db
                    if 0 <= na < n and 0 <= nb < n and not seen[na, nb] and field[na, nb] > 0.5:
                        seen[na, nb] = True
                        stack.append((na, nb))
    return {"on_fraction": float(np.mean(field)), "on_components": int(comps)}


def code_hamming(kc: SparseKC, a: np.ndarray, b: np.ndarray) -> int:
    kc.reset()
    ha = kc.code(a) > 0
    kc.reset()
    hb = kc.code(b) > 0
    return int(np.sum(ha != hb))


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    plant = PlantParams()
    at = LumpedAT(ATParams())
    affine = AffineComposition(AffineParams())
    kc = SparseKC(KCParams(), rng)

    # Train KC readout on a composition-aware in-silico target (not AT, not a clinic).
    ss = rng.uniform(0.02, 0.85, size=800)
    rr = rng.uniform(0.02, 0.85, size=800)
    train_states = np.stack([ss, rr], axis=1)
    targets = np.array([gatenby_like_target(st) for st in train_states])
    kc.fit_readout(train_states, targets)

    ics = {
        "matched_T_sensitive": np.array([0.48, 0.12]),  # T=0.60, S-rich
        "matched_T_resistant": np.array([0.12, 0.48]),  # T=0.60, R-rich
        "mid_mix": np.array([0.30, 0.20]),
    }

    controllers = {"lumped_at": at, "affine_sr": affine, "sparse_kc": kc}

    isoline = {}
    regions = {}
    closed = {}
    for name, ctl in controllers.items():
        isoline[name] = {
            "T0.60": isoline_u_variance(ctl, 0.60),
            "T0.40": isoline_u_variance(ctl, 0.40),
        }
        regions[name] = decision_region_complexity(ctl)
        closed[name] = {}
        for ic_name, ic in ics.items():
            out = simulate(ctl, ic, plant)
            closed[name][ic_name] = {
                k: (None if v is None else (v if not isinstance(v, np.ndarray) else None))
                for k, v in out.items()
                if k not in {"t", "s", "r", "u"}
            }
            closed[name][ic_name]["t_end"] = float(out["t"][-1])
            closed[name][ic_name]["s_end"] = float(out["s"][-1])
            closed[name][ic_name]["r_end"] = float(out["r"][-1])

    # Immediate action at matched burden
    immediate = {}
    for name, ctl in controllers.items():
        immediate[name] = {}
        for ic_name, ic in ics.items():
            ctl.reset()
            immediate[name][ic_name] = float(ctl.act(ic))

    # Retune AT on the S-rich closed-loop burden of KC, then transfer to R-rich IC.
    kc.reset()
    kc_srich = simulate(kc, ics["matched_T_sensitive"], plant)
    tuned = retune_at(kc_srich["t"])
    at_tuned = LumpedAT(tuned)
    transfer = {
        "tuned_thresholds": asdict(tuned),
        "srich": {},
        "rrich": {},
    }
    for label, ic_key in (("srich", "matched_T_sensitive"), ("rrich", "matched_T_resistant")):
        kc.reset()
        ref = simulate(kc, ics[ic_key], plant)
        tun = simulate(at_tuned, ics[ic_key], plant)
        n = min(len(ref["t"]), len(tun["t"]))
        transfer[label] = {
            "rmse_T": float(np.sqrt(np.mean((ref["t"][:n] - tun["t"][:n]) ** 2))),
            "rmse_u": float(np.sqrt(np.mean((ref["u"][:n] - tun["u"][:n]) ** 2))),
            "kc_duty": ref["duty"],
            "at_duty": tun["duty"],
            "kc_u0": ref["u0"],
            "at_u0": tun["u0"],
            "kc_time_r_dom": ref["time_r_dom"],
            "at_time_r_dom": tun["time_r_dom"],
        }

    # Isoline disagreement with lumped AT (reset maps; hysteresis band uses u=0)
    def isoline_series(controller, t_value: float, n: int = 81) -> np.ndarray:
        fracs = np.linspace(0.05, 0.95, n)
        us = []
        for f in fracs:
            controller.reset()
            us.append(controller.act(np.array([t_value * f, t_value * (1.0 - f)])))
        return np.linspace(0.05, 0.95, n), np.array(us, dtype=float)

    fracs60, u_at_60 = isoline_series(at, 0.60)
    _, u_aff_60 = isoline_series(affine, 0.60)
    _, u_kc_60 = isoline_series(kc, 0.60)

    def pearson(x, y) -> float:
        if np.std(x) < 1e-12 or np.std(y) < 1e-12:
            return 0.0
        return float(np.corrcoef(x, y)[0, 1])

    isoline_id = {
        "T0.60_pearson_u_vs_sensitive_fraction": {
            "lumped_at": pearson(fracs60, u_at_60),
            "affine_sr": pearson(fracs60, u_aff_60),
            "sparse_kc": pearson(fracs60, u_kc_60),
        },
        "T0.60_fraction_disagree_AT_threshold_0.5": {
            "affine_sr": float(np.mean((u_aff_60 > 0.5) != (u_at_60 > 0.5))),
            "sparse_kc": float(np.mean((u_kc_60 > 0.5) != (u_at_60 > 0.5))),
        },
        "T0.60_mean_abs_u_minus_AT": {
            "affine_sr": float(np.mean(np.abs(u_aff_60 - u_at_60))),
            "sparse_kc": float(np.mean(np.abs(u_kc_60 - u_at_60))),
        },
    }

    # Nearby-state code distance (pattern separation signature)
    a = np.array([0.50, 0.10])
    b = np.array([0.10, 0.50])  # same T, opposite composition
    kc.reset()
    ua = kc.act(a.copy())
    kc.reset()
    ub = kc.act(b.copy())
    affine.reset()
    aa = affine.act(a)
    affine.reset()
    ab = affine.act(b)
    at.reset()
    ata = at.act(a)
    at.reset()
    atb = at.act(b)
    separation = {
        "delta_T": float(abs(total_burden(a) - total_burden(b))),
        "kc_hamming_active": code_hamming(kc, a, b),
        "kc_du": float(abs(ua - ub)),
        "affine_du": float(abs(aa - ab)),
        "at_du": float(abs(ata - atb)),
        "kc_ua": ua,
        "kc_ub": ub,
        "affine_ua": aa,
        "affine_ub": ab,
        "at_ua": ata,
        "at_ub": atb,
    }

    summary = {
        "seed": RNG_SEED,
        "plant": asdict(plant),
        "disclaimer": "In-silico controller-architecture comparison on a toy ODE. Not treatment.",
        "immediate_u_at_matched_T": immediate,
        "isoline_u": isoline,
        "isoline_identifiability": isoline_id,
        "decision_regions": regions,
        "closed_loop": closed,
        "at_retune_transfer": transfer,
        "nearby_state_separation": separation,
    }

    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "results.json"
    out_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
