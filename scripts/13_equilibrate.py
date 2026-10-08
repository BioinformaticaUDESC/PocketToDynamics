#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path.home() / "PocketToDynamics"
CONFIG = ROOT / "config/dynamics.yaml"

DEFAULT_INPUT = Path(
    "/mnt/labio-fast-disk/"
    "PocketToDynamics/dynamics/"
    "minimized/"
    "minimization_manifest.tsv"
)

DEFAULT_LOG = (
    ROOT /
    "logs/M13_equilibration.log"
)


def load_tsv(path: Path):

    if not path.exists():
        return []

    with path.open() as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def run(
    command,
    cwd: Path,
    env=None,
):

    return subprocess.run(
        command,
        cwd=cwd,
        capture_output=True,
        text=True,
        env=env,
    )


def parse_gmx_energy_xvg(
    path: Path,
):

    values = []

    for line in path.read_text(
        errors="replace"
    ).splitlines():

        if not line:
            continue

        if line.startswith(
            ("#", "@")
        ):
            continue

        parts = line.split()

        if len(parts) < 2:
            continue

        try:
            values.append(
                float(parts[1])
            )
        except ValueError:
            continue

    if not values:
        return None

    mean = (
        sum(values)
        /
        len(values)
    )

    var = (
        sum(
            (x - mean) ** 2
            for x in values
        )
        /
        len(values)
    )

    return {
        "mean":
            mean,

        "sd":
            math.sqrt(var),

        "min":
            min(values),

        "max":
            max(values),

        "n":
            len(values),
    }


def extract_energy(
    gmx: Path,
    edr: Path,
    term: str,
    output: Path,
    cwd: Path,
):

    proc = subprocess.run(
        [
            str(gmx),
            "energy",
            "-f",
            str(edr),
            "-o",
            str(output),
        ],
        input=f"{term}\n0\n",
        capture_output=True,
        text=True,
        cwd=cwd,
    )

    return proc


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_INPUT,
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
    )

    parser.add_argument(
        "--log",
        type=Path,
        default=DEFAULT_LOG,
    )

    args = parser.parse_args()

    cfg = yaml.safe_load(
        CONFIG.read_text()
    )

    ecfg = cfg["equilibration"]

    manifest = (
        args.manifest
        .expanduser()
        .resolve()
    )

    rows = load_tsv(
        manifest
    )

    log_path = (
        args.log
        .expanduser()
        .resolve()
    )

    log_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if args.output_root:

        output_root = (
            args.output_root
            .expanduser()
            .resolve()
        )

    else:

        output_root = Path(
            ecfg[
                "heavy_output_root"
            ]
        )

    output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not rows:

        status = (
            "SKIP_M13_NO_MINIMIZED_SYSTEMS"
        )

        log_path.write_text(
            status
            + "\nSystems: 0\n"
        )

        print(status)
        print("Systems: 0")

        return

    gmx = Path(
        ecfg["gmx_executable"]
    )

    if not gmx.exists():
        raise SystemExit(
            "FAIL_M13_GROMACS_NOT_FOUND"
        )

    failures = 0
    results = []

    for row in rows:

        if row.get("status") != "PASS":
            continue

        job_id = (
            row["dynamics_job_id"]
        )

        em_gro = Path(
            row["em_gro"]
        )

        source_job = (
            em_gro.parent
        )

        topology = (
            source_job /
            "system_solvated.top"
        )

        if not (
            em_gro.exists()
            and
            topology.exists()
        ):

            print(
                f"[M13] {job_id} "
                "FAIL missing M12 files"
            )

            failures += 1
            continue

        job_dir = (
            output_root /
            job_id
        )

        if job_dir.exists():
            job_dir.rename(job_dir.with_name(job_dir.name + '.previous-' + str(__import__('time').time_ns())))

        job_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        top_copy = (
            job_dir /
            "system.top"
        )

        shutil.copy2(
            topology,
            top_copy,
        )

        ntmpi = int(
            ecfg["ntmpi"]
        )

        ntomp = int(
            ecfg["ntomp"]
        )

        env = os.environ.copy()

        env[
            "OMP_NUM_THREADS"
        ] = str(
            ntomp
        )

        dt = float(
            ecfg["timestep_ps"]
        )

        nvt_steps = int(
            round(
                float(
                    ecfg["nvt"][
                        "duration_ps"
                    ]
                )
                / dt
            )
        )

        npt_steps = int(
            round(
                float(
                    ecfg["npt"][
                        "duration_ps"
                    ]
                )
                / dt
            )
        )

        try:

            nvt_mdp = (
                job_dir /
                "nvt.mdp"
            )

            nvt_mdp.write_text(
                f"""\
integrator               = md
dt                       = {dt}
nsteps                   = {nvt_steps}

continuation             = no
constraint_algorithm     = lincs
constraints              = {ecfg['constraints']}

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
coulombtype              = PME
vdwtype                  = Cut-off

tcoupl                   = {ecfg['thermostat']}
tc-grps                  = System
tau_t                    = {float(ecfg['nvt']['tau_t_ps'])}
ref_t                    = {float(ecfg['temperature_k'])}

pcoupl                   = no

gen_vel                  = yes
gen_temp                 = {float(ecfg['temperature_k'])}
gen_seed                 = 20261008

nstxout                   = 0
nstvout                   = 0
nstfout                   = 0
nstxout-compressed        = {int(ecfg['nstxout_compressed'])}
nstenergy                = {int(ecfg['nstenergy'])}
nstlog                   = {int(ecfg['nstlog'])}

pbc                      = xyz
"""
            )

            nvt_tpr = (
                job_dir /
                "nvt.tpr"
            )

            grompp_nvt = run(
                [
                    str(gmx),
                    "grompp",
                    "-f",
                    str(nvt_mdp),
                    "-c",
                    str(em_gro),
                    "-p",
                    str(top_copy),
                    "-o",
                    str(nvt_tpr),
                    "-maxwarn",
                    "0",
                ],
                job_dir,
            )

            (
                job_dir /
                "grompp_nvt.stdout"
            ).write_text(
                grompp_nvt.stdout
                or ""
            )

            (
                job_dir /
                "grompp_nvt.stderr"
            ).write_text(
                grompp_nvt.stderr
                or ""
            )

            if (
                grompp_nvt.returncode != 0
                or
                not nvt_tpr.exists()
            ):

                raise RuntimeError(
                    "NVT grompp failed"
                )

            nvt_prefix = (
                job_dir /
                "nvt"
            )

            nvt_run = run(
                [
                    str(gmx),
                    "mdrun",
                    "-deffnm",
                    str(nvt_prefix),
                    "-ntmpi",
                    str(ntmpi),
                    "-ntomp",
                    str(ntomp),
                ],
                job_dir,
                env=env,
            )

            (
                job_dir /
                "mdrun_nvt.stdout"
            ).write_text(
                nvt_run.stdout
                or ""
            )

            (
                job_dir /
                "mdrun_nvt.stderr"
            ).write_text(
                nvt_run.stderr
                or ""
            )

            nvt_gro = (
                job_dir /
                "nvt.gro"
            )

            nvt_cpt = (
                job_dir /
                "nvt.cpt"
            )

            nvt_edr = (
                job_dir /
                "nvt.edr"
            )

            if (
                nvt_run.returncode != 0
                or
                not nvt_gro.exists()
                or
                not nvt_cpt.exists()
                or
                not nvt_edr.exists()
            ):

                raise RuntimeError(
                    "NVT mdrun failed"
                )

            temp_xvg = (
                job_dir /
                "nvt_temperature.xvg"
            )

            temp_proc = (
                extract_energy(
                    gmx,
                    nvt_edr,
                    "Temperature",
                    temp_xvg,
                    job_dir,
                )
            )

            if (
                temp_proc.returncode
                != 0
                or
                not temp_xvg.exists()
            ):

                raise RuntimeError(
                    "NVT temperature "
                    "extraction failed"
                )

            temp_stats = (
                parse_gmx_energy_xvg(
                    temp_xvg
                )
            )

            if temp_stats is None:

                raise RuntimeError(
                    "No NVT temperature data"
                )

            npt_mdp = (
                job_dir /
                "npt.mdp"
            )

            npt_mdp.write_text(
                f"""\
integrator               = md
dt                       = {dt}
nsteps                   = {npt_steps}

continuation             = yes
constraint_algorithm     = lincs
constraints              = {ecfg['constraints']}

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
coulombtype              = PME
vdwtype                  = Cut-off

tcoupl                   = {ecfg['thermostat']}
tc-grps                  = System
tau_t                    = {float(ecfg['npt']['tau_t_ps'])}
ref_t                    = {float(ecfg['temperature_k'])}

pcoupl                   = {ecfg['barostat']}
pcoupltype               = isotropic
tau_p                    = {float(ecfg['npt']['tau_p_ps'])}
ref_p                    = {float(ecfg['pressure_bar'])}
compressibility          = {float(ecfg['npt']['compressibility_bar_inv'])}

gen_vel                  = no

nstxout                   = 0
nstvout                   = 0
nstfout                   = 0
nstxout-compressed        = {int(ecfg['nstxout_compressed'])}
nstenergy                = {int(ecfg['nstenergy'])}
nstlog                   = {int(ecfg['nstlog'])}

pbc                      = xyz
"""
            )

            npt_tpr = (
                job_dir /
                "npt.tpr"
            )

            grompp_npt = run(
                [
                    str(gmx),
                    "grompp",
                    "-f",
                    str(npt_mdp),
                    "-c",
                    str(nvt_gro),
                    "-t",
                    str(nvt_cpt),
                    "-p",
                    str(top_copy),
                    "-o",
                    str(npt_tpr),
                    "-maxwarn",
                    "0",
                ],
                job_dir,
            )

            (
                job_dir /
                "grompp_npt.stdout"
            ).write_text(
                grompp_npt.stdout
                or ""
            )

            (
                job_dir /
                "grompp_npt.stderr"
            ).write_text(
                grompp_npt.stderr
                or ""
            )

            if (
                grompp_npt.returncode != 0
                or
                not npt_tpr.exists()
            ):

                raise RuntimeError(
                    "NPT grompp failed"
                )

            npt_prefix = (
                job_dir /
                "npt"
            )

            npt_run = run(
                [
                    str(gmx),
                    "mdrun",
                    "-deffnm",
                    str(npt_prefix),
                    "-ntmpi",
                    str(ntmpi),
                    "-ntomp",
                    str(ntomp),
                ],
                job_dir,
                env=env,
            )

            (
                job_dir /
                "mdrun_npt.stdout"
            ).write_text(
                npt_run.stdout
                or ""
            )

            (
                job_dir /
                "mdrun_npt.stderr"
            ).write_text(
                npt_run.stderr
                or ""
            )

            npt_gro = (
                job_dir /
                "npt.gro"
            )

            npt_cpt = (
                job_dir /
                "npt.cpt"
            )

            npt_edr = (
                job_dir /
                "npt.edr"
            )

            if (
                npt_run.returncode != 0
                or
                not npt_gro.exists()
                or
                not npt_cpt.exists()
                or
                not npt_edr.exists()
            ):

                raise RuntimeError(
                    "NPT mdrun failed"
                )

            density_xvg = (
                job_dir /
                "npt_density.xvg"
            )

            pressure_xvg = (
                job_dir /
                "npt_pressure.xvg"
            )

            temp2_xvg = (
                job_dir /
                "npt_temperature.xvg"
            )

            dproc = extract_energy(
                gmx,
                npt_edr,
                "Density",
                density_xvg,
                job_dir,
            )

            pproc = extract_energy(
                gmx,
                npt_edr,
                "Pressure",
                pressure_xvg,
                job_dir,
            )

            tproc = extract_energy(
                gmx,
                npt_edr,
                "Temperature",
                temp2_xvg,
                job_dir,
            )

            if (
                dproc.returncode != 0
                or
                pproc.returncode != 0
                or
                tproc.returncode != 0
            ):

                raise RuntimeError(
                    "NPT energy extraction "
                    "failed"
                )

            density_stats = (
                parse_gmx_energy_xvg(
                    density_xvg
                )
            )

            pressure_stats = (
                parse_gmx_energy_xvg(
                    pressure_xvg
                )
            )

            temp2_stats = (
                parse_gmx_energy_xvg(
                    temp2_xvg
                )
            )

            if (
                density_stats is None
                or
                pressure_stats is None
                or
                temp2_stats is None
            ):

                raise RuntimeError(
                    "Incomplete NPT metrics"
                )

            metadata = {
                "dynamics_job_id":
                    job_id,

                "temperature_target_k":
                    float(
                        ecfg[
                            "temperature_k"
                        ]
                    ),

                "pressure_target_bar":
                    float(
                        ecfg[
                            "pressure_bar"
                        ]
                    ),

                "nvt": {
                    "duration_ps":
                        float(
                            ecfg["nvt"][
                                "duration_ps"
                            ]
                        ),

                    "temperature":
                        temp_stats,
                },

                "npt": {
                    "duration_ps":
                        float(
                            ecfg["npt"][
                                "duration_ps"
                            ]
                        ),

                    "temperature":
                        temp2_stats,

                    "pressure":
                        pressure_stats,

                    "density":
                        density_stats,
                },

                "parallelization": {
                    "ntmpi":
                        ntmpi,

                    "ntomp":
                        ntomp,

                    "omp_num_threads":
                        ntomp,
                },

                "files": {
                    "topology":
                        str(
                            top_copy
                        ),

                    "nvt_gro":
                        str(
                            nvt_gro
                        ),

                    "nvt_cpt":
                        str(
                            nvt_cpt
                        ),

                    "npt_gro":
                        str(
                            npt_gro
                        ),

                    "npt_cpt":
                        str(
                            npt_cpt
                        ),

                    "npt_edr":
                        str(
                            npt_edr
                        ),
                },

                "status":
                    "PASS",
            }

            (
                job_dir /
                "equilibration_manifest.json"
            ).write_text(
                json.dumps(
                    metadata,
                    indent=2,
                )
                + "\n"
            )

            results.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "nvt_temp_mean":
                        temp_stats[
                            "mean"
                        ],

                    "npt_temp_mean":
                        temp2_stats[
                            "mean"
                        ],

                    "npt_pressure_mean":
                        pressure_stats[
                            "mean"
                        ],

                    "npt_density_mean":
                        density_stats[
                            "mean"
                        ],

                    "npt_gro":
                        str(
                            npt_gro
                        ),

                    "status":
                        "PASS",
                }
            )

            print(
                f"[M13] {job_id} PASS "
                f"NVT_T="
                f"{temp_stats['mean']:.2f} "
                f"NPT_T="
                f"{temp2_stats['mean']:.2f} "
                f"rho="
                f"{density_stats['mean']:.2f}"
            )

        except Exception as exc:

            failures += 1

            results.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "nvt_temp_mean":
                        "",

                    "npt_temp_mean":
                        "",

                    "npt_pressure_mean":
                        "",

                    "npt_density_mean":
                        "",

                    "npt_gro":
                        "",

                    "status":
                        "FAIL",
                }
            )

            print(
                f"[M13] {job_id} FAIL "
                f"{type(exc).__name__}: "
                f"{exc}"
            )

    table = (
        output_root /
        "equilibration_manifest.tsv"
    )

    fields = [
        "dynamics_job_id",
        "nvt_temp_mean",
        "npt_temp_mean",
        "npt_pressure_mean",
        "npt_density_mean",
        "npt_gro",
        "status",
    ]

    with table.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            results
        )

    passed = sum(
        r["status"] == "PASS"
        for r in results
    )

    if failures:

        status = (
            "FAIL_M13_EQUILIBRATION"
        )

    elif passed == 0:

        status = (
            "SKIP_M13_NO_VALID_SYSTEMS"
        )

    else:

        status = (
            "PASS_M13_EQUILIBRATION"
        )

    log_path.write_text(
        "\n".join(
            [
                status,
                f"Systems: {len(rows)}",
                f"PASS: {passed}",
                f"FAIL: {failures}",
                (
                    "Output root: "
                    f"{output_root}"
                ),
            ]
        )
        + "\n"
    )

    print()
    print(status)
    print(
        "Systems:",
        len(rows),
    )
    print(
        "PASS:",
        passed,
    )
    print(
        "FAIL:",
        failures,
    )

    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
