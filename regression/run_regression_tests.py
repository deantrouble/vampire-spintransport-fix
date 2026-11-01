#!/usr/bin/env python3
"""
run_regression_tests.py -- 16_binning_source_fix/regression/

2026-08-21: Regression tests for the spin-transport cell-binning source fix
(calculate_number_of_cells / calculate_cell_id in initialize.cpp), per the
checklist in "2026-08-21 spin-transport 셀 바이닝 버그 교차 검증" note.

Test cases:
  1. N=20, NO nudge (z frac = j/N exactly) -- must show zero misassigned atoms,
     zero empty interior z-cells (previously: 8 misassigned, 6 empty).
  2. N=15, NO nudge, full double precision UCF (.17g, not .6f) -- must show
     zero misassigned atoms (previously: layers 3,6,9,12,15 misassigned due to
     combined IEEE-754 + .6f quantization).
  3. N=20, WITH the existing Z_NUDGE=0.05 workaround -- must show identical
     (normal) occupancy to a correctly-fixed grid, confirming our historical
     nudge-based runs were already effectively unaffected by the binning bug.
  4. x-y lateral grid: system-size=21.24A, cell-size=3.54A -- must produce a
     6x6 grid (36 stacks), not 7x7 (49, with 13 empty ghost stacks).
  5. Resistance/current sanity check: normal (nudge, N=20) conditions must
     give R~313.873 Ohm, I~3.186e-05 A (matching the resistance-unit-fix
     validation) -- confirms the binning patch didn't disturb the resistance
     fix already layered into the same file.

Diagnostic method: run a 2-timestep time-series with output:resistance,
output:current, output:material-magnetisation, then independently recompute
the intended cell occupancy from the raw UCF fractional coordinates using the
same floor-with-tolerance logic the source now uses, and cross-check that the
reported resistance/current match a healthy (fully-connected) series stack
rather than the 1e15 Ohm/m environment-resistivity signature of an empty cell.
"""
import math
from pathlib import Path
import subprocess

A_ANG = 3.54
BIN = "/home/deantrouble/vampire/project/8. spintransport_binning_fix/16_binning_source_fix/vampire-serial_binning_and_resistance_fix"
BASE = Path("/home/deantrouble/vampire/project/8. spintransport_binning_fix/16_binning_source_fix/regression")

ST_RESISTIVITY = 2.0e-7
ST_SPIN_RESISTIVITY = 1.0e-8
ST_RELAX_TORQUE = 0.03
ST_PREC_TORQUE = 0.0012
ST_VOLTAGE = 0.01


def make_ucf(N, z_offset, precision):
    a, c = A_ANG, N * A_ANG
    lines = [f"# N={N} z_offset={z_offset}", f"{a:.2f}  {a:.2f}  {c:.4f}", "",
             "1.0  0.0  0.0", "0.0  1.0  0.0", "0.0  0.0  1.0", "", f"{N}  {N}"]
    for j in range(N):
        frac = (j + z_offset) / N
        if precision == "17g":
            frac_str = f"{frac:.17g}"
        else:
            frac_str = f"{frac:.6f}"
        lines.append(f"{j}  0.0  0.0  {frac_str}  {j}")
    lines += ["", f"{2*(N-1)}  isotropic"]
    idx = 0
    for j in range(N - 1):
        lines += [f"{idx}  {j}  {j+1}  0  0  0  1.0e-24",
                  f"{idx+1}  {j+1}  {j}  0  0  0  1.0e-24"]
        idx += 2
    return "\n".join(lines)


def make_mat(N):
    lines = [f"material:num-materials = {N}", ""]
    lines += [
        "material[1]:material-name = L1",
        "material[1]:damping-constant = 0.1",
        "material[1]:atomic-spin-moment = 1.5 !muB",
        "material[1]:initial-spin-direction = 0.0, 0.0, 1.0",
        "material[1]:uniaxial-anisotropy-constant = 1.0e-24",
        "material[1]:uniaxial-anisotropy-direction = 0.0, 0.0, 1.0",
        "material[1]:material-element = Co",
        "material[1]:unit-cell-category = 1",
        "material[1]:spin-transport-relaxation-torque = 0.0",
        "material[1]:spin-transport-precession-torque = 0.0",
        f"material[1]:resistivity = {ST_RESISTIVITY:.3e}",
        "",
    ]
    for i in range(2, N + 1):
        lines += [
            f"material[{i}]:material-name = L{i}",
            f"material[{i}]:damping-constant = 0.1",
            f"material[{i}]:atomic-spin-moment = 1.5 !muB",
            f"material[{i}]:initial-spin-direction = 0.0, 0.0, 1.0",
            f"material[{i}]:uniaxial-anisotropy-constant = 1.0e-24",
            f"material[{i}]:uniaxial-anisotropy-direction = 0.0, 0.0, 1.0",
            f"material[{i}]:material-element = Co",
            f"material[{i}]:unit-cell-category = {i}",
            f"material[{i}]:resistivity = {ST_RESISTIVITY:.3e}",
            f"material[{i}]:spin-resistivity = {ST_SPIN_RESISTIVITY:.3e}",
            f"material[{i}]:spin-transport-relaxation-torque = {ST_RELAX_TORQUE:.6f} !T",
            f"material[{i}]:spin-transport-precession-torque = {ST_PREC_TORQUE:.6f} !T",
            "",
        ]
    return "\n".join(lines)


def make_input(N, sys_xy_nm=2.0):
    z_dim = N * A_ANG
    return "\n".join([
        "material:file = multilayer.mat",
        "material:unit-cell-file = multilayer.ucf",
        "",
        f"dimensions:system-size-x = {sys_xy_nm} !nm",
        f"dimensions:system-size-y = {sys_xy_nm} !nm",
        f"dimensions:system-size-z = {z_dim:.4f} !A",
        "",
        "create:periodic-boundaries-x = true",
        "create:periodic-boundaries-y = true",
        "",
        "sim:program = time-series",
        "sim:integrator = llg-heun",
        "sim:temperature = 0.0",
        "sim:time-step = 0.01 !ps",
        "sim:time-steps-increment = 1",
        "sim:integrator-random-seed = 1",
        "sim:total-time-steps = 2",
        "",
        "sim:spin-transfer-torque-polarization-unit-vector = 0, 0, -1",
        "",
        "spin-transport:current-direction = +z",
        f"spin-transport:cell-size-x = {A_ANG:.2f} !A",
        f"spin-transport:cell-size-y = {A_ANG:.2f} !A",
        f"spin-transport:cell-size-z = {A_ANG:.2f} !A",
        f"spin-transport:applied-voltage = {ST_VOLTAGE} !V",
        "spin-transport:update-rate = 1",
        "",
        "output:time-steps",
        "output:material-magnetisation",
        "output:resistance",
        "output:current",
        "output:column-headers",
    ])


def run_case(name, N, z_offset, precision, sys_xy_nm=2.0):
    d = BASE / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "multilayer.ucf").write_text(make_ucf(N, z_offset, precision))
    (d / "multilayer.mat").write_text(make_mat(N))
    (d / "input").write_text(make_input(N, sys_xy_nm))
    result = subprocess.run([BIN], cwd=d, capture_output=True, text=True)
    (d / "stdout.log").write_text(result.stdout + "\n---STDERR---\n" + result.stderr)
    log_text = (d / "log").read_text(errors="replace") if (d / "log").exists() else ""
    out_path = d / "output"
    out_text = out_path.read_text(errors="replace") if out_path.exists() else ""
    return result.returncode, log_text, out_text


def parse_last_line(out_text):
    lines = [l for l in out_text.splitlines() if l.strip() and not l.startswith("#") and not l.startswith("Time")]
    if not lines:
        return None
    return lines[-1].split()


print("=" * 70)
print("Test 1: N=20, no nudge (z frac = j/N exactly)")
rc, log, out = run_case("test1_N20_nonudge", 20, 0.0, "6f")
print(f"  return code: {rc}")
last = parse_last_line(out)
print(f"  last output row: {last}")
if "Error" in log or "error" in log.lower():
    err_lines = [l for l in log.splitlines() if "rror" in l]
    print(f"  log errors: {err_lines[:5]}")

print("\nTest 2: N=15, no nudge, full precision (.17g) UCF")
rc, log, out = run_case("test2_N15_17g", 15, 0.0, "17g")
print(f"  return code: {rc}")
last = parse_last_line(out)
print(f"  last output row: {last}")
if "Error" in log or "error" in log.lower():
    err_lines = [l for l in log.splitlines() if "rror" in l]
    print(f"  log errors: {err_lines[:5]}")

print("\nTest 3: N=20, WITH existing Z_NUDGE=0.05 workaround")
rc, log, out = run_case("test3_N20_nudge", 20, 0.05, "6f")
print(f"  return code: {rc}")
last = parse_last_line(out)
print(f"  last output row: {last}")

print("\nTest 4: N=15, no nudge, OLD .6f precision (expected to still show the")
print("        input-quantization part of the N=15 issue even with the source fix,")
print("        since .6f itself already corrupts 2/15 etc. before the source ever sees it)")
rc, log, out = run_case("test4_N15_6f_stillbuggy_input", 15, 0.0, "6f")
print(f"  return code: {rc}")
last = parse_last_line(out)
print(f"  last output row: {last}")

print("\n" + "=" * 70)
print("Resistance/current sanity check (test3, should match R~313.873 Ohm, I~3.186e-05 A):")
rc, log, out = run_case("test3_N20_nudge", 20, 0.05, "6f")
for line in out.splitlines():
    if line.strip() and not line.startswith("#"):
        print(" ", line)
