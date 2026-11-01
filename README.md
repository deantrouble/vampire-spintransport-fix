# VAMPIRE spin-transport fix and simulation inputs

Modified VAMPIRE source and input files used for the atomistic spin-dynamics simulations in
"Ultralow-current-induced spin-transfer-torque switching in van der Waals magnetic tunnel junctions".

## What this repository contains

| Path | Content |
|---|---|
| `patches/spintransport_fix.patch` | Change relative to upstream VAMPIRE (git diff, one file) |
| `src/spintransport/initialize.cpp` | The modified file in full (drop-in replacement) |
| `make_inputs.py` | Generates the `input`, `multilayer.mat` and `multilayer.ucf` files for the N-layer simulations |
| `examples/N20_Jinter0.3_tau1_J58mT/` | Example input set (N = 20, J_inter = 0.3 erg/cm^2, tau = 1, J = 58.12 mT, Case 1) |
| `regression/run_regression_tests.py` | Regression tests for the cell assignment and resistance (paths must be adjusted) |
| `license/` | Licence files of upstream VAMPIRE (kept unchanged) |

## Base version

Upstream: https://github.com/richard-evans/vampire, commit `525bc27e` ("Updates to the manual for v7.0").
All other source files are identical to upstream.

## What was changed (spin-transport module only)

1. **Resistance unit.** The cell dimensions are stored in angstrom while the resistivity is given in ohm metre.
   The original code computed rho * l / A with the dimensions in angstrom, so the cell resistance was too large
   by a factor of 10^10 (A / l = 1e-10 m). The geometry is now converted to metre before computing the resistance.
2. **Assignment of atoms to transport cells.** The cell index was obtained by truncating coord / cell_size.
   For coordinates that lie exactly on a cell boundary (e.g. layer spacing = cell size = 3.54 A), floating-point
   rounding can give a value just below the integer and the atom is assigned to the cell below, which can leave
   cells empty and break the series circuit. The index is now computed with a small tolerance (`snap_near_integer`),
   and coordinates outside the cell grid raise an error instead of being silently clamped.

The torque, spin-resistance model, and all other modules are unchanged.

## Build

```
git clone https://github.com/richard-evans/vampire.git
cd vampire
git checkout 525bc27e
git apply /path/to/patches/spintransport_fix.patch      # or copy src/spintransport/initialize.cpp
make serial
```

The binary used in this work has SHA-256 starting with `32c6570a1a79...`.

## Run an example

```
cd examples/N20_Jinter0.3_tau1_J58mT
/path/to/vampire-serial
```

`make_inputs.py <directory>` regenerates these files. Edit `N`, `ji`, `J_mT`, `tau` at the bottom of the script.

Output: hysteresis loop (`output`), from +0.5 T to -0.5 T and back, in steps of 10 mT with 1 ns equilibration
and 2 ns per field step. The coercive field is taken from the descending branch. J_c is the STT amplitude J
at which the coercive field vanishes (see Methods of the paper).

## Key parameters (Case 1)

| Parameter | Value |
|---|---|
| Lattice, size | a = 3.54 A, 2 x 2 nm^2 periodic, N layers along z |
| Moment, anisotropy, damping | 1.5 mu_B, K_u = 1e-24 J/atom, alpha = 0.1 |
| Temperature, time step | 0.01 K, 10 fs (Heun) |
| Interlayer exchange | J_ij = J_inter a^2, J_inter = 0.01 - 10 erg/cm^2 |
| External STT | polarization -z, amplitude J exp[-(n-1)/tau], field-like 0.01 T exp[-(n-1)/tau] |
| Internal STT | spin-transport module, V = 0.01 V, rho_ep = 2e-7 Ohm m, rho_sp = 1e-4 Ohm m, prefactors 0.03 T and 0.0012 T |
| Initial state | all spins +z, layer 1 tilted by 1 degree |

## Citation

If you use this code, please cite the paper above and the original VAMPIRE reference:
R. F. L. Evans et al., J. Phys.: Condens. Matter 26, 103202 (2014).

## Licence

The modified source file keeps the original header, which states that VAMPIRE is released under the FreeBSD
(BSD 2-clause) licence (copyright Richard F L Evans). The upstream repository also contains a GPL v2 text
file; both licence files are included unchanged in `license/`. This repository is distributed under the same
terms as upstream VAMPIRE.
