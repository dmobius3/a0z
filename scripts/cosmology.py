"""
cosmology.py
============

Flat LambdaCDM expansion history under Planck 2018 parameters, as specified in
§3 of the paper:

    E(z) = sqrt( Omega_m (1+z)^3 + Omega_r (1+z)^4 + Omega_Lambda )      (3.2)
    H(z) = H_0 * E(z)
    t_age(z) = integral_z^infty dz' / [(1+z') H(z')]                     (3.3)

Planck 2018 parameters (paper §3, citation [6]):
    Omega_m      = 0.315
    Omega_r      = 9.2e-5
    Omega_Lambda = 0.685
    H_0          = 67.4 km/s/Mpc

Numerical conventions (per scripts/README.md "Coding standards"):
    - Cosmic-age integral uses scipy.integrate.quad with epsrel = 1e-8.
    - No fitting; deterministic given the parameters above.

Verification target (per scripts/README.md):
    E(z = 2) = 3.0327
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy import integrate


# ---------------------------------------------------------------------------
# Planck 2018 parameters (§3, citation [6])
# ---------------------------------------------------------------------------

OMEGA_M = 0.315
OMEGA_R = 9.2e-5
OMEGA_LAMBDA = 0.685
H0_KM_S_MPC = 67.4

# Unit conversions for the cosmic-age integral.
# 1 Mpc = 3.085677581491367e19 km, so (km/s/Mpc)^-1 in seconds is
# 1 Mpc / (1 km/s) = 3.085677581491367e19 s.
KM_PER_MPC = 3.085677581491367e19
SECONDS_PER_GYR = 3.15576e16  # Julian Gyr (365.25 d/yr)
HUBBLE_TIME_GYR = (KM_PER_MPC / SECONDS_PER_GYR) / H0_KM_S_MPC  # 1/H0 in Gyr


# ---------------------------------------------------------------------------
# Parameter container
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Cosmology:
    """Flat LambdaCDM with radiation. Defaults are Planck 2018 (paper §3, citation [6])."""
    Omega_m: float = OMEGA_M
    Omega_r: float = OMEGA_R
    Omega_Lambda: float = OMEGA_LAMBDA
    H0: float = H0_KM_S_MPC  # km/s/Mpc

    def E(self, z):
        """Dimensionless Hubble parameter E(z) = H(z)/H0, eq. (3.2)."""
        z = np.asarray(z, dtype=float)
        one_plus_z = 1.0 + z
        return np.sqrt(
            self.Omega_m * one_plus_z**3
            + self.Omega_r * one_plus_z**4
            + self.Omega_Lambda
        )

    def H(self, z):
        """H(z) in km/s/Mpc."""
        return self.H0 * self.E(z)

    def t_age(self, z):
        """
        Cosmic age at redshift z in Gyr, via
            t_age(z) = (1/H0) * integral_z^infty dz' / [(1+z') E(z')].
        """
        z_in = np.asarray(z, dtype=float)
        scalar_input = z_in.ndim == 0
        z_arr = np.atleast_1d(z_in)
        out = np.empty_like(z_arr)
        for i, zi in enumerate(z_arr):
            integrand = lambda zp: 1.0 / ((1.0 + zp) * float(self.E(zp)))
            value, _ = integrate.quad(
                integrand, float(zi), np.inf, epsrel=1e-8, limit=200
            )
            out[i] = value
        # Convert (1/H0) [in units (km/s/Mpc)^-1] to Gyr:
        out_gyr = out * HUBBLE_TIME_GYR
        if scalar_input:
            return float(out_gyr[0])
        return out_gyr


# Module-level default instance (Planck 2018).
PLANCK18 = Cosmology()


# Convenience top-level functions wrapping the default cosmology.
def E(z):
    return PLANCK18.E(z)


def H(z):
    return PLANCK18.H(z)


def t_age(z):
    return PLANCK18.t_age(z)


# ---------------------------------------------------------------------------
# Self-check / verification
# ---------------------------------------------------------------------------

def _verification_table(redshifts: Iterable[float] = (0.0, 0.5, 1.0, 2.0, 3.0,
                                                       5.0, 10.0, 15.0)):
    """Reproduce the columns of paper Table 1 that this module is responsible
    for: z, E(z), H(z), t_age(z). The downstream a_0(z) and enhancement
    columns belong to the framework module."""
    rows = []
    for z in redshifts:
        Ez = float(np.asarray(E(z)))
        Hz = float(np.asarray(H(z)))
        tz = float(t_age(z))
        rows.append((z, Ez, Hz, tz))
    return rows


def main():
    # Print the four-column subset of Table 1 covered by this module.
    print(f"{'z':>5} {'E(z)':>10} {'H(z) [km/s/Mpc]':>18} {'t_age [Gyr]':>14}")
    for z, Ez, Hz, t in _verification_table():
        print(f"{z:5.1f} {Ez:10.4f} {Hz:18.2f} {t:14.3f}")

    # Verification target from scripts/README.md: E(2) = 3.0327
    target = 3.0327
    computed = float(E(2.0))
    rounded = round(computed, 4)
    print()
    print(f"Verification target  E(z=2) = {target}")
    print(f"Computed             E(z=2) = {computed:.10f}")
    print(f"Rounded to 4 dp      E(z=2) = {rounded}")
    print(f"Match (4 dp)         : {rounded == target}")


if __name__ == "__main__":
    main()
