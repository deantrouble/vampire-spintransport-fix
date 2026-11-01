#!/usr/bin/env python3
"""Input-file generators for the Case 1 (External + Internal STT) hysteresis simulations.
Same functions as used for the production runs; see README.md."""
import math

A_ANG      = 3.54
ALPHA      = 0.1
TILT_DEG   = 1.0
STT_PREC   = 0.01
ERCM2_TO_J = 1e-3 * (A_ANG * 1e-10)**2

N       = 20
TAU     = 1.0
SCALE   = 10000.0

ST_RESISTIVITY           = 2.0e-7
ST_SPIN_RESISTIVITY_BASE = 1.0e-8
ST_RELAX_TORQUE          = 0.03
ST_PREC_TORQUE           = 0.0012
ST_VOLTAGE               = 0.01



def make_ucf(N, J_inter):
    a, c = A_ANG, N * A_ANG
    lines = [f"# N={N} J_inter={J_inter:.3e} (tau=1 ultrafine, scale=10000x)",
             f"{a:.2f}  {a:.2f}  {c:.4f}", "",
             "1.0  0.0  0.0", "0.0  1.0  0.0", "0.0  0.0  1.0", "",
             f"{N}  {N}"]
    for j in range(N):
        frac = j / N
        lines.append(f"{j}  0.0  0.0  {frac:.17g}  {j}")
    lines += ["", f"{2*(N-1)}  isotropic"]
    idx = 0
    for j in range(N - 1):
        lines += [f"{idx}  {j}  {j+1}  0  0  0  {J_inter:.6e}",
                  f"{idx+1}  {j+1}  {j}  0  0  0  {J_inter:.6e}"]
        idx += 2
    return "\n".join(lines)


def make_mat(N, J_scan, tau, spin_res):
    tilt_x = math.sin(math.radians(TILT_DEG))
    lines = [f"material:num-materials = {N}", ""]
    for i in range(1, N + 1):
        w = math.exp(-(i - 1) / tau)
        rl = J_scan * w
        pr = STT_PREC * w
        init_dir = f"{tilt_x:.6f}, 0.0, 1.0" if i == 1 else "0.0, 0.0, 1.0"
        lines += [
            f"material[{i}]:material-name = FGT-L{i}",
            f"material[{i}]:damping-constant = {ALPHA}",
            "material[{i}]:atomic-spin-moment = 1.5 !muB".format(i=i),
            f"material[{i}]:initial-spin-direction = {init_dir}",
            "material[{i}]:uniaxial-anisotropy-constant = 1.0e-24".format(i=i),
            f"material[{i}]:uniaxial-anisotropy-direction = 0.0, 0.0, 1.0",
            f"material[{i}]:material-element = Co",
            f"material[{i}]:unit-cell-category = {i}",
            f"material[{i}]:spin-transfer-relaxation-torque = {rl:.8f} !T",
            f"material[{i}]:spin-transfer-precession-torque = {pr:.8f} !T",
            f"material[{i}]:resistivity = {ST_RESISTIVITY:.3e}",
        ]
        if i == 1:
            lines += [
                "material[1]:spin-transport-relaxation-torque = 0.0",
                "material[1]:spin-transport-precession-torque = 0.0",
            ]
        else:
            lines += [
                f"material[{i}]:spin-resistivity = {spin_res:.3e}",
                f"material[{i}]:spin-transport-relaxation-torque = {ST_RELAX_TORQUE:.6f} !T",
                f"material[{i}]:spin-transport-precession-torque = {ST_PREC_TORQUE:.6f} !T",
            ]
        lines.append("")
    return "\n".join(lines)


def make_input(N):
    z_dim = N * A_ANG
    return "\n".join([
        "material:file = multilayer.mat",
        "material:unit-cell-file = multilayer.ucf",
        "",
        "dimensions:system-size-x = 2 !nm",
        "dimensions:system-size-y = 2 !nm",
        f"dimensions:system-size-z = {z_dim:.4f} !A",
        "",
        "create:periodic-boundaries-x = true",
        "create:periodic-boundaries-y = true",
        "",
        "sim:program = hysteresis-loop",
        "sim:integrator = llg-heun",
        "sim:temperature = 0.01",
        "sim:time-step = 0.01 !ps",
        "sim:time-steps-increment = 100",
        "sim:integrator-random-seed = 1",
        "",
        "sim:equilibration-time-steps = 100000",
        "sim:loop-time-steps = 200000",
        "sim:maximum-applied-field-strength = 0.5000 !T",
        "sim:applied-field-strength-increment = 0.0100 !T",
        "sim:applied-field-unit-vector = 0.0, 0.0, 1.0",
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
        "output:applied-field-strength",
        "output:magnetisation",
        "output:material-magnetisation",
        "output:column-headers",
    ])




if __name__ == "__main__":
    import sys
    from pathlib import Path
    N, ji, J_mT, tau = 20, 0.3, 58.12, 1.0
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "run")
    out.mkdir(exist_ok=True)
    (out / "multilayer.ucf").write_text(make_ucf(N, ji * ERCM2_TO_J))
    (out / "multilayer.mat").write_text(make_mat(N, J_mT * 1e-3, tau, ST_SPIN_RESISTIVITY_BASE * SCALE))
    (out / "input").write_text(make_input(N))
    print("written to", out)
