#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import parmed as pmd
import yaml


ROOT = Path.home() / "PocketToDynamics"
CONFIG = ROOT / "config/dynamics.yaml"

DEFAULT_INPUT = Path(
    "/mnt/labio-fast-disk/"
    "PocketToDynamics/dynamics/"
    "parameterized/"
    "parameterization_manifest.tsv"
)

DEFAULT_LOG = (
    ROOT /
    "logs/M12_minimization.log"
)


def sha256(path: Path) -> str:

    h = hashlib.sha256()

    with path.open("rb") as f:

        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


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


def parse_em_log(path: Path):

    text = path.read_text(
        errors="replace"
    )

    potential = None
    fmax = None
    converged = False

    potential_matches = re.findall(
        r"Potential Energy\s*=\s*"
        r"([-+0-9.eE]+)",
        text,
    )

    if potential_matches:

        potential = float(
            potential_matches[-1]
        )

    fmax_matches = re.findall(
        r"(?:Fmax|Maximum force)\s*=\s*"
        r"([-+0-9.eE]+)",
        text,
        re.IGNORECASE,
    )

    if fmax_matches:

        fmax = float(
            fmax_matches[-1]
        )

    if re.search(
        r"converged to Fmax",
        text,
        re.IGNORECASE,
    ):

        converged = True

    return {
        "potential_energy_kj_mol":
            potential,

        "fmax_kj_mol_nm":
            fmax,

        "converged":
            converged,
    }


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

    setup = cfg["system_setup"]
    emcfg = cfg["energy_minimization"]
    pcfg = cfg["parameterization"]

    manifest_path = (
        args.manifest
        .expanduser()
        .resolve()
    )

    rows = load_tsv(
        manifest_path
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
            emcfg[
                "heavy_output_root"
            ]
        )

    output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not rows:

        status = (
            "SKIP_M12_NO_PARAMETERIZED_SYSTEMS"
        )

        log_path.write_text(
            status
            + "\nSystems: 0\n"
        )

        print(status)
        print("Systems: 0")

        return

    amber_root = Path(
        pcfg["ambertools_root"]
    )

    tleap = (
        amber_root /
        "bin/tleap"
    )

    gmx = Path(
        emcfg[
            "gmx_executable"
        ]
    )

    if not tleap.exists():

        raise SystemExit(
            "FAIL_M12_TLEAP_NOT_FOUND"
        )

    if not gmx.exists():

        raise SystemExit(
            "FAIL_M12_GROMACS_NOT_FOUND"
        )

    results = []
    failures = 0

    for row in rows:

        if row.get("status") != "PASS":
            continue

        job_id = (
            row["dynamics_job_id"]
        )

        source_top = Path(
            row["system_top"]
        )

        source_gro = Path(
            row["system_gro"]
        )

        source_job = (
            source_top.parent
        )

        receptor_clean = (
            source_job /
            "receptor_clean.pdb"
        )

        ligand_gaff2 = (
            source_job /
            "ligand_gaff2.mol2"
        )

        ligand_frcmod = (
            source_job /
            "ligand.frcmod"
        )

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

        solv_prmtop = (
            job_dir /
            "solvated.prmtop"
        )

        solv_inpcrd = (
            job_dir /
            "solvated.inpcrd"
        )

        solv_pdb = (
            job_dir /
            "solvated.pdb"
        )

        gmx_top = (
            job_dir /
            "system_solvated.top"
        )

        gmx_gro = (
            job_dir /
            "system_solvated.gro"
        )

        leap_input = (
            job_dir /
            "tleap_solvate.in"
        )

        em_mdp = (
            job_dir /
            "em.mdp"
        )

        em_tpr = (
            job_dir /
            "em.tpr"
        )

        em_prefix = (
            job_dir /
            "em"
        )

        try:

            required_m11 = [
                receptor_clean,
                ligand_gaff2,
                ligand_frcmod,
            ]

            missing_m11 = [
                str(path)
                for path in required_m11
                if not path.exists()
            ]

            if missing_m11:

                raise FileNotFoundError(
                    "Missing M11 components: "
                    + ", ".join(missing_m11)
                )

            leap_lines = [
                (
                    "source "
                    f"leaprc.protein."
                    f"{pcfg['protein_forcefield']}"
                ),
                (
                    "source "
                    f"leaprc."
                    f"{pcfg['ligand_forcefield']}"
                ),
                (
                    "source "
                    f"leaprc.water."
                    f"{setup['water_model']}"
                ),
                "",
                (
                    "REC = loadPDB "
                    f"{receptor_clean}"
                ),
                (
                    "LIG = loadMol2 "
                    f"{ligand_gaff2}"
                ),
                (
                    "loadAmberParams "
                    f"{ligand_frcmod}"
                ),
                "",
                "COM = combine { REC LIG }",
                "",
            ]

            shape = (
                setup[
                    "box_shape"
                ].lower()
            )

            buffer_a = float(
                setup[
                    "buffer_angstrom"
                ]
            )

            if shape == "octahedral":

                leap_lines.append(
                    "solvateOct COM "
                    f"{setup['water_box']} "
                    f"{buffer_a}"
                )

            elif shape == "rectangular":

                leap_lines.append(
                    "solvateBox COM "
                    f"{setup['water_box']} "
                    f"{buffer_a}"
                )

            else:

                raise ValueError(
                    "Unsupported box_shape: "
                    f"{shape}"
                )

            if setup.get(
                "neutralize",
                True,
            ):

                leap_lines += [
                    (
                        "addIonsRand COM "
                        f"{setup['positive_ion']} "
                        "0"
                    ),
                    (
                        "addIonsRand COM "
                        f"{setup['negative_ion']} "
                        "0"
                    ),
                ]

            leap_lines += [
                "",
                "check COM",
                (
                    "saveAmberParm COM "
                    f"{solv_prmtop} "
                    f"{solv_inpcrd}"
                ),
                (
                    "savePDB COM "
                    f"{solv_pdb}"
                ),
                "",
                "quit",
            ]

            leap_input.write_text(
                "\n".join(
                    leap_lines
                )
                + "\n"
            )

            leap_cmd = [
                str(tleap),
                "-f",
                str(leap_input),
            ]

            leap = run(
                leap_cmd,
                job_dir,
            )

            (
                job_dir /
                "tleap_solvate.stdout"
            ).write_text(
                leap.stdout or ""
            )

            (
                job_dir /
                "tleap_solvate.stderr"
            ).write_text(
                leap.stderr or ""
            )

            if (
                leap.returncode != 0
                or
                not solv_prmtop.exists()
                or
                not solv_inpcrd.exists()
            ):

                raise RuntimeError(
                    "tleap solvation failed"
                )

            amber = pmd.load_file(
                str(solv_prmtop),
                str(solv_inpcrd),
            )

            amber.save(
                str(gmx_top),
                format="gromacs",
                overwrite=True,
            )

            amber.save(
                str(gmx_gro),
                overwrite=True,
            )

            natoms = len(
                amber.atoms
            )

            nresidues = len(
                amber.residues
            )

            charge = sum(
                atom.charge
                for atom in amber.atoms
            )

            em_mdp.write_text(
                f"""\
integrator               = {emcfg['integrator']}
emtol                    = {float(emcfg['emtol_kj_mol_nm'])}
emstep                   = {float(emcfg['emstep_nm'])}
nsteps                   = {int(emcfg['nsteps'])}

cutoff-scheme            = {emcfg['cutoff_scheme']}
nstlist                  = {int(emcfg['nstlist'])}
rlist                    = {float(emcfg['rlist_nm'])}
coulombtype              = {emcfg['coulombtype']}
rcoulomb                 = {float(emcfg['rcoulomb_nm'])}
vdwtype                  = Cut-off
rvdw                     = {float(emcfg['rvdw_nm'])}

constraints              = none
pbc                      = xyz
"""
            )

            grompp_cmd = [
                str(gmx),
                "grompp",

                "-f",
                str(em_mdp),

                "-c",
                str(gmx_gro),

                "-p",
                str(gmx_top),

                "-o",
                str(em_tpr),

                "-maxwarn",
                "0",
            ]

            grompp = run(
                grompp_cmd,
                job_dir,
            )

            (
                job_dir /
                "grompp.stdout"
            ).write_text(
                grompp.stdout or ""
            )

            (
                job_dir /
                "grompp.stderr"
            ).write_text(
                grompp.stderr or ""
            )

            if (
                grompp.returncode != 0
                or
                not em_tpr.exists()
            ):

                raise RuntimeError(
                    "gmx grompp failed"
                )

            mdrun_cmd = [
                str(gmx),
                "mdrun",

                "-deffnm",
                str(em_prefix),

                "-ntmpi",
                str(
                    int(
                        emcfg["ntmpi"]
                    )
                ),

                "-ntomp",
                str(
                    int(
                        emcfg["ntomp"]
                    )
                ),
            ]

            mdrun_env = os.environ.copy()

            mdrun_env[
                "OMP_NUM_THREADS"
            ] = str(
                int(
                    emcfg["ntomp"]
                )
            )

            mdrun = run(
                mdrun_cmd,
                job_dir,
                env=mdrun_env,
            )

            (
                job_dir /
                "mdrun.stdout"
            ).write_text(
                mdrun.stdout or ""
            )

            (
                job_dir /
                "mdrun.stderr"
            ).write_text(
                mdrun.stderr or ""
            )

            em_gro = (
                job_dir /
                "em.gro"
            )

            em_log = (
                job_dir /
                "em.log"
            )

            if (
                mdrun.returncode != 0
                or
                not em_gro.exists()
                or
                not em_log.exists()
            ):

                raise RuntimeError(
                    "gmx mdrun minimization failed"
                )

            em_metrics = (
                parse_em_log(
                    em_log
                )
            )

            potential = (
                em_metrics[
                    "potential_energy_kj_mol"
                ]
            )

            fmax = (
                em_metrics[
                    "fmax_kj_mol_nm"
                ]
            )

            if (
                potential is None
                or
                not math.isfinite(
                    potential
                )
            ):

                raise RuntimeError(
                    "Non-finite or missing "
                    "minimized potential energy"
                )

            metadata = {
                "dynamics_job_id":
                    job_id,

                "water_model":
                    setup[
                        "water_model"
                    ],

                "water_box":
                    setup[
                        "water_box"
                    ],

                "box_shape":
                    shape,

                "buffer_angstrom":
                    buffer_a,

                "neutralize":
                    bool(
                        setup.get(
                            "neutralize",
                            True,
                        )
                    ),

                "total_atoms":
                    natoms,

                "total_residues":
                    nresidues,

                "net_charge_after_ions":
                    charge,

                "minimization": {
                    **em_metrics,

                    "emtol_kj_mol_nm":
                        float(
                            emcfg[
                                "emtol_kj_mol_nm"
                            ]
                        ),

                    "nsteps_max":
                        int(
                            emcfg[
                                "nsteps"
                            ]
                        ),

                    "ntmpi":
                        int(
                            emcfg[
                                "ntmpi"
                            ]
                        ),

                    "ntomp":
                        int(
                            emcfg[
                                "ntomp"
                            ]
                        ),

                    "omp_num_threads":
                        int(
                            emcfg[
                                "ntomp"
                            ]
                        ),
                },

                "commands": {
                    "tleap":
                        leap_cmd,

                    "grompp":
                        grompp_cmd,

                    "mdrun":
                        mdrun_cmd,
                },

                "files": {
                    "solvated_prmtop":
                        str(
                            solv_prmtop
                        ),

                    "solvated_inpcrd":
                        str(
                            solv_inpcrd
                        ),

                    "gromacs_top":
                        str(
                            gmx_top
                        ),

                    "gromacs_gro":
                        str(
                            gmx_gro
                        ),

                    "em_tpr":
                        str(
                            em_tpr
                        ),

                    "em_gro":
                        str(
                            em_gro
                        ),

                    "em_log":
                        str(
                            em_log
                        ),
                },

                "sha256": {
                    "solvated_prmtop":
                        sha256(
                            solv_prmtop
                        ),

                    "solvated_inpcrd":
                        sha256(
                            solv_inpcrd
                        ),

                    "em_gro":
                        sha256(
                            em_gro
                        ),
                },

                "status":
                    "PASS",
            }

            (
                job_dir /
                "minimization_manifest.json"
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

                    "total_atoms":
                        natoms,

                    "net_charge":
                        f"{charge:.6f}",

                    "potential_energy":
                        potential,

                    "fmax":
                        (
                            ""
                            if fmax is None
                            else fmax
                        ),

                    "converged":
                        int(
                            em_metrics[
                                "converged"
                            ]
                        ),

                    "em_gro":
                        str(
                            em_gro
                        ),

                    "status":
                        "PASS",
                }
            )

            print(
                f"[M12] {job_id} PASS "
                f"atoms={natoms} "
                f"charge={charge:.3f} "
                f"Epot={potential:.3f}"
            )

        except Exception as exc:

            failures += 1

            results.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "total_atoms":
                        "",

                    "net_charge":
                        "",

                    "potential_energy":
                        "",

                    "fmax":
                        "",

                    "converged":
                        0,

                    "em_gro":
                        "",

                    "status":
                        "FAIL",
                }
            )

            print(
                f"[M12] {job_id} FAIL "
                f"{type(exc).__name__}: "
                f"{exc}"
            )

    result_table = (
        output_root /
        "minimization_manifest.tsv"
    )

    fields = [
        "dynamics_job_id",
        "total_atoms",
        "net_charge",
        "potential_energy",
        "fmax",
        "converged",
        "em_gro",
        "status",
    ]

    with result_table.open(
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
            "FAIL_M12_MINIMIZATION"
        )

    elif passed == 0:

        status = (
            "SKIP_M12_NO_VALID_SYSTEMS"
        )

    else:

        status = (
            "PASS_M12_MINIMIZATION"
        )

    log_path.write_text(
        "\n".join(
            [
                status,
                f"Systems: {len(rows)}",
                f"PASS: {passed}",
                f"FAIL: {failures}",
                f"Output root: {output_root}",
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
