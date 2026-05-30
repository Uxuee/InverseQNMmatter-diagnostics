from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


"""Residual diagnostics relative to an equatorial Kerr photon-ring baseline.

The rotating eikonal baseline uses the azimuthal phase number m:

    omega_QNM = m*Omega_phi - i*(n + 1/2)*lambda.

The geodesic photon-ring frequency Omega_phi is signed by the branch
convention. Residuals are computed from positive QNM frequency magnitudes,
using Re(omega_QNM)/abs(m) and abs(Omega_phi).
"""


# Demo/default inputs. Units use G = c = 1.
M = 1.0
a = 0.7
ell = 2
m = 2
n = 0
branch = 1  # +1 prograde, -1 retrograde


def validate_kerr_spin(M: float, a: float) -> None:
    if M <= 0.0:
        raise ValueError("M must be positive.")
    if abs(a) > M:
        raise ValueError("Kerr spin must satisfy abs(a) <= M.")


def validate_branch(branch: int) -> None:
    if branch not in (-1, 1):
        raise ValueError("branch must be +1 for prograde or -1 for retrograde.")


def validate_azimuthal_number(m: int, branch: int) -> None:
    if m == 0:
        raise ValueError(
            "The equatorial Kerr eikonal baseline requires nonzero azimuthal number m."
        )
    if np.sign(m) != branch:
        raise ValueError(
            "Use sign(m) consistent with the branch: m > 0 for prograde, "
            "m < 0 for retrograde."
        )


def kerr_photon_radius(M: float, a: float, branch: int) -> float:
    validate_kerr_spin(M, a)
    validate_branch(branch)
    argument = np.clip(-branch * a / M, -1.0, 1.0)
    return 2.0 * M * (1.0 + np.cos((2.0 / 3.0) * np.arccos(argument)))


def kerr_photon_frequency_signed(M: float, a: float, branch: int) -> float:
    r_ph = kerr_photon_radius(M, a, branch)
    sqrt_M = np.sqrt(M)
    return branch * sqrt_M / (r_ph**1.5 + branch * a * sqrt_M)


def kerr_delta(r: float, M: float, a: float) -> float:
    return r**2 - 2.0 * M * r + a**2


def kerr_impact_parameter(M: float, a: float, branch: int) -> float:
    r_ph = kerr_photon_radius(M, a, branch)
    if abs(a) < 1.0e-8 * M:
        return branch * np.sqrt(27.0) * M

    numerator = r_ph**3 - 3.0 * M * r_ph**2 + a**2 * r_ph + a**2 * M
    denominator = a * (r_ph - M)
    return -numerator / denominator


def kerr_radial_potential(r: float, M: float, a: float, impact_parameter: float) -> float:
    E = 1.0
    L = impact_parameter
    delta = kerr_delta(r, M, a)
    return ((r**2 + a**2) * E - a * L) ** 2 - delta * (L - a * E) ** 2


def second_derivative_r(fun, r: float, step: float) -> float:
    return (fun(r + step) - 2.0 * fun(r) + fun(r - step)) / step**2


def kerr_dt_dtau(r: float, M: float, a: float, impact_parameter: float) -> float:
    """Coordinate-time geodesic factor dt/dtau for equatorial null Kerr orbits."""

    E = 1.0
    L = impact_parameter
    delta = kerr_delta(r, M, a)
    P = E * (r**2 + a**2) - a * L
    return ((r**2 + a**2) * P / delta + a * (L - a * E)) / r**2


def kerr_lyapunov_exponent(M: float, a: float, branch: int) -> float:
    """Numerically compute coordinate-time Lyapunov exponent for the photon ring.

    The expression follows lambda_K = sqrt(R''/(2*r_ph^4*(dt/dtau)^2)).
    The Schwarzschild limit is handled smoothly by using b = +/-sqrt(27) M.
    """

    r_ph = kerr_photon_radius(M, a, branch)
    b = kerr_impact_parameter(M, a, branch)
    step = max(1.0e-5 * M, 1.0e-6 * r_ph)
    radial_potential = lambda radius: kerr_radial_potential(radius, M, a, b)
    R_second = second_derivative_r(radial_potential, r_ph, step)
    dt_dtau = kerr_dt_dtau(r_ph, M, a, b)
    lambda_sq = R_second / (2.0 * r_ph**4 * dt_dtau**2)

    if lambda_sq < -1.0e-10:
        raise RuntimeError(
            "Computed negative lambda_K^2. Check radial-potential sign convention."
        )
    return float(np.sqrt(max(lambda_sq, 0.0)))


def extract_equatorial_kerr_observables(
    omega_qnm: complex,
    m: int,
    n: int,
) -> tuple[float, float]:
    """Extract positive-frequency eikonal observables for a Kerr baseline.

    The equatorial Kerr eikonal relation uses the azimuthal number m:
    omega_QNM = m*Omega_phi - i*(n+1/2)*lambda.
    For residual diagnostics we compare positive observed QNM frequency
    magnitudes, so Omega_obs = Re(omega_qnm)/abs(m).
    """

    if m == 0:
        raise ValueError(
            "The equatorial Kerr eikonal baseline requires nonzero azimuthal number m."
        )
    if n <= -0.5:
        raise ValueError("n must satisfy n + 1/2 > 0.")
    Omega_obs = np.real(omega_qnm) / abs(m)
    if Omega_obs < 0.0:
        raise ValueError("Expected positive Re(omega_qnm) for residual diagnostics.")
    lambda_obs = -np.imag(omega_qnm) / (n + 0.5)
    return Omega_obs, lambda_obs


def compute_kerr_residuals(
    omega_qnm: complex,
    M: float,
    a: float,
    m: int,
    n: int,
    branch: int,
) -> dict[str, float]:
    validate_branch(branch)
    validate_azimuthal_number(m, branch)
    Omega_obs, lambda_obs = extract_equatorial_kerr_observables(omega_qnm, m, n)
    r_ph = kerr_photon_radius(M, a, branch)
    Omega_Kerr_signed = kerr_photon_frequency_signed(M, a, branch)
    Omega_Kerr_abs = abs(Omega_Kerr_signed)
    lambda_Kerr = kerr_lyapunov_exponent(M, a, branch)
    return {
        "Omega_obs": Omega_obs,
        "lambda_obs": lambda_obs,
        "r_ph": r_ph,
        "Omega_Kerr_signed": Omega_Kerr_signed,
        "Omega_Kerr_abs": Omega_Kerr_abs,
        "lambda_Kerr": lambda_Kerr,
        "A_Kerr": Omega_obs / Omega_Kerr_abs - 1.0,
        "B_Kerr": lambda_obs / lambda_Kerr - 1.0,
    }


def synthetic_residual_demo(
    M: float,
    a: float,
    ell: int,
    m: int,
    n: int,
    branch: int,
    A_Kerr_injected: float = 0.01,
    B_Kerr_injected: float = -0.02,
) -> pd.DataFrame:
    validate_azimuthal_number(m, branch)
    Omega_Kerr_signed = kerr_photon_frequency_signed(M, a, branch)
    Omega_Kerr_abs = abs(Omega_Kerr_signed)
    lambda_Kerr = kerr_lyapunov_exponent(M, a, branch)
    Omega_obs = Omega_Kerr_abs * (1.0 + A_Kerr_injected)
    lambda_obs = lambda_Kerr * (1.0 + B_Kerr_injected)
    omega_qnm = abs(m) * Omega_obs - 1j * (n + 0.5) * lambda_obs
    recovered = compute_kerr_residuals(omega_qnm, M, a, m, n, branch)
    schwarzschild_expected = 1.0 / (3.0 * np.sqrt(3.0) * M)
    schwarzschild_lambda = kerr_lyapunov_exponent(M, 0.0, branch)

    return pd.DataFrame(
        [
            {
                "M": M,
                "a": a,
                "a_over_M": a / M,
                "ell": ell,
                "m": m,
                "n": n,
                "branch": branch,
                "omega_re": omega_qnm.real,
                "omega_im": omega_qnm.imag,
                "Omega_Kerr_signed": Omega_Kerr_signed,
                "Omega_Kerr_abs": Omega_Kerr_abs,
                "Omega_obs": recovered["Omega_obs"],
                "A_Kerr_injected": A_Kerr_injected,
                "B_Kerr_injected": B_Kerr_injected,
                "A_Kerr_recovered": recovered["A_Kerr"],
                "B_Kerr_recovered": recovered["B_Kerr"],
                "A_abs_error": abs(recovered["A_Kerr"] - A_Kerr_injected),
                "B_abs_error": abs(recovered["B_Kerr"] - B_Kerr_injected),
                "r_ph": recovered["r_ph"],
                "lambda_Kerr": recovered["lambda_Kerr"],
                "schwarzschild_lambda_expected": schwarzschild_expected,
                "schwarzschild_lambda_numeric": schwarzschild_lambda,
                "schwarzschild_lambda_abs_error": abs(
                    schwarzschild_lambda - schwarzschild_expected
                ),
            }
        ]
    )


def scan_kerr_baselines(M: float) -> pd.DataFrame:
    rows = []
    for branch_value in [1, -1]:
        for spin in np.linspace(0.0, 0.98, 80):
            a_value = spin * M
            rows.append(
                {
                    "a_over_M": spin,
                    "branch": branch_value,
                    "branch_label": "prograde" if branch_value == 1 else "retrograde",
                    "r_ph_over_M": kerr_photon_radius(M, a_value, branch_value) / M,
                    "M_Omega_phi_signed": M
                    * kerr_photon_frequency_signed(M, a_value, branch_value),
                    "M_Omega_K_abs": M
                    * abs(kerr_photon_frequency_signed(M, a_value, branch_value)),
                    "M_lambda_K": M * kerr_lyapunov_exponent(M, a_value, branch_value),
                }
            )
    return pd.DataFrame(rows)


def schwarzschild_limit_check(M: float) -> pd.DataFrame:
    expected = 1.0 / (3.0 * np.sqrt(3.0))
    lambda_numeric = M * kerr_lyapunov_exponent(M, 0.0, branch=1)
    return pd.DataFrame(
        [
            {
                "M": M,
                "a": 0.0,
                "a_over_M": 0.0,
                "M_lambda_Kerr_numeric": lambda_numeric,
                "M_lambda_Kerr_expected": expected,
                "abs_error": abs(lambda_numeric - expected),
            }
        ]
    )


def plot_branch_scan(
    scan: pd.DataFrame,
    y_column: str,
    ylabel: str,
    title: str,
    filename: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5), constrained_layout=True)
    for label, group in scan.groupby("branch_label"):
        ordered = group.sort_values("a_over_M")
        ax.plot(ordered["a_over_M"], ordered[y_column], label=label)
    if y_column == "M_lambda_K":
        ax.axhline(
            1.0 / (3.0 * np.sqrt(3.0)),
            color="black",
            linestyle="--",
            linewidth=1.0,
            label="Schwarzschild limit",
        )
    ax.set_xlabel("a/M")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.savefig(filename, dpi=180)
    plt.close(fig)


def plot_synthetic_recovery(demo: pd.DataFrame, filename: Path) -> None:
    labels = ["A_Kerr", "B_Kerr"]
    injected = [demo["A_Kerr_injected"].iloc[0], demo["B_Kerr_injected"].iloc[0]]
    recovered = [demo["A_Kerr_recovered"].iloc[0], demo["B_Kerr_recovered"].iloc[0]]
    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    ax.bar(x - width / 2.0, injected, width, label="injected")
    ax.bar(x + width / 2.0, recovered, width, label="recovered")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Kerr residual shift")
    ax.set_title("Synthetic Kerr residual recovery")
    ax.legend()
    fig.savefig(filename, dpi=180)
    plt.close(fig)


def main() -> None:
    out_dir = Path("outputs") / "rotating_kerr_residuals"
    out_dir.mkdir(parents=True, exist_ok=True)

    demo = synthetic_residual_demo(M, a, ell, m, n, branch)
    demo.to_csv(out_dir / "kerr_residual_demo.csv", index=False)
    sanity = schwarzschild_limit_check(M)
    sanity.to_csv(out_dir / "schwarzschild_limit_check.csv", index=False)

    scan = scan_kerr_baselines(M)
    plot_branch_scan(
        scan,
        "r_ph_over_M",
        "r_ph/M",
        "Equatorial Kerr photon-ring radius",
        out_dir / "kerr_photon_radius_vs_spin.png",
    )
    plot_branch_scan(
        scan,
        "M_Omega_phi_signed",
        "M * Omega_phi (signed)",
        "Equatorial Kerr photon-ring frequency",
        out_dir / "kerr_photon_frequency_vs_spin.png",
    )
    plot_branch_scan(
        scan,
        "M_lambda_K",
        "M * lambda_K",
        "Equatorial Kerr Lyapunov exponent",
        out_dir / "kerr_lyapunov_vs_spin.png",
    )
    plot_synthetic_recovery(demo, out_dir / "synthetic_kerr_residual_recovery.png")

    print("Wrote rotating Kerr residual outputs to:")
    print(f"  {out_dir}")
    print()
    print(demo.to_string(index=False))
    print()
    print("Schwarzschild-limit check:")
    print(sanity.to_string(index=False))


if __name__ == "__main__":
    main()
