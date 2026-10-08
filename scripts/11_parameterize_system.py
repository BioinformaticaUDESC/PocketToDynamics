#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import yaml
import parmed as pmd


ROOT = Path.home() / "PocketToDynamics"

CONFIG = ROOT / "config/dynamics.yaml"

DEFAULT_INPUT = (
    ROOT /
    "work/dynamics/prepared/"
    "pose_extraction_manifest.tsv"
)

DEFAULT_LOG = (
    ROOT /
    "logs/M11_parameterization.log"
)


def sha256(path: Path) -> str:

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

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
):

    proc = subprocess.run(
        command,
        cwd=cwd,
        capture_output=True,
        text=True,
    )

    return proc


def prepare_receptor_pdb(
    source: Path,
    output: Path,
    strip_waters: bool,
    water_names: set[str],
):
    """
    Prepare receptor coordinates for force-field parameterization.

    Only crystallographic waters explicitly listed in water_names are
    removed when strip_waters=True. Other HETATM records are preserved
    so that cofactors are never silently discarded.
    """

    kept = []
    removed_waters = 0
    removed_atoms = 0

    removed_residues = set()

    for line in source.read_text(
        errors="replace"
    ).splitlines():

        if not (
            line.startswith("ATOM")
            or
            line.startswith("HETATM")
        ):
            kept.append(line)
            continue

        resname = line[17:20].strip()

        if (
            strip_waters
            and
            resname in water_names
        ):

            removed_atoms += 1

            chain = line[21:22].strip()
            resid = line[22:26].strip()

            removed_residues.add(
                (
                    resname,
                    chain,
                    resid,
                )
            )

            continue

        kept.append(line)

    output.write_text(
        "\n".join(kept)
        + "\n"
    )

    removed_waters = len(
        removed_residues
    )

    return {
        "removed_water_residues":
            removed_waters,

        "removed_water_atoms":
            removed_atoms,
    }


def formal_charge_for_ligand(
    ligand_id: str,
) -> int:

    table = (
        ROOT /
        "work/chemoinformatics/"
        "ligand_descriptors.tsv"
    )

    rows = load_tsv(table)

    for row in rows:

        if row["ligand_id"] == ligand_id:

            return int(
                float(
                    row["formal_charge"]
                )
            )

    # Safe technical default only when descriptor data are unavailable.
    return 0


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

    pcfg = cfg["parameterization"]

    receptor_cfg = cfg.get(
        "receptor_preparation",
        {},
    )

    strip_waters = bool(
        receptor_cfg.get(
            "strip_crystallographic_waters",
            True,
        )
    )

    water_names = set(
        receptor_cfg.get(
            "water_residue_names",
            ["HOH", "WAT", "SOL"],
        )
    )

    amber_root = Path(
        pcfg["ambertools_root"]
    )

    antechamber = (
        amber_root /
        "bin/antechamber"
    )

    parmchk2 = (
        amber_root /
        "bin/parmchk2"
    )

    tleap = (
        amber_root /
        "bin/tleap"
    )

    required_tools = {
        "antechamber": antechamber,
        "parmchk2": parmchk2,
        "tleap": tleap,
    }

    missing = [
        name
        for name, path
        in required_tools.items()
        if not path.exists()
    ]

    if missing:

        print(
            "FAIL_M11_MISSING_TOOLS:",
            ",".join(missing),
        )

        raise SystemExit(1)

    manifest = (
        args.manifest.expanduser().resolve()
    )

    rows = load_tsv(
        manifest
    )

    if args.output_root:

        output_root = (
            args.output_root
            .expanduser()
            .resolve()
        )

    else:

        output_root = Path(
            pcfg[
                "heavy_output_root"
            ]
        )

    output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    log_path = (
        args.log.expanduser().resolve()
    )

    log_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not rows:

        status = (
            "SKIP_M11_NO_EXTRACTED_POSES"
        )

        log_path.write_text(
            status
            + "\nSystems: 0\n"
        )

        print(status)
        print("Systems: 0")

        return

    result_rows = []
    failures = 0

    for row in rows:

        if row.get("status") != "PASS":
            continue

        job_id = (
            row["dynamics_job_id"]
        )

        ligand_id = (
            row["ligand_id"]
        )

        receptor = Path(
            row["receptor_pdb"]
        )

        ligand_input = Path(
            row["selected_pose_mol2"]
        )

        if not receptor.is_absolute():
            receptor = ROOT / receptor

        if not ligand_input.is_absolute():
            ligand_input = ROOT / ligand_input

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

        receptor_copy = (
            job_dir /
            "receptor_clean.pdb"
        )

        ligand_source = (
            job_dir /
            "ligand_input.mol2"
        )

        receptor_qc = (
            prepare_receptor_pdb(
                source=receptor,
                output=receptor_copy,
                strip_waters=strip_waters,
                water_names=water_names,
            )
        )

        shutil.copy2(
            ligand_input,
            ligand_source,
        )

        print(
            f"[M11] {job_id} "
            f"receptor waters_removed="
            f"{receptor_qc['removed_water_residues']} "
            f"water_atoms_removed="
            f"{receptor_qc['removed_water_atoms']}"
        )

        charge = (
            formal_charge_for_ligand(
                ligand_id
            )
        )

        ligand_gaff = (
            job_dir /
            "ligand_gaff2.mol2"
        )

        frcmod = (
            job_dir /
            "ligand.frcmod"
        )

        leap_input = (
            job_dir /
            "tleap.in"
        )

        prmtop = (
            job_dir /
            "complex.prmtop"
        )

        inpcrd = (
            job_dir /
            "complex.inpcrd"
        )

        complex_pdb = (
            job_dir /
            "complex_amber.pdb"
        )

        gmx_top = (
            job_dir /
            "system.top"
        )

        gmx_gro = (
            job_dir /
            "system.gro"
        )

        manifest_json = (
            job_dir /
            "parameterization_manifest.json"
        )

        try:

            ac_cmd = [
                str(antechamber),

                "-i",
                str(ligand_source),

                "-fi",
                "mol2",

                "-o",
                str(ligand_gaff),

                "-fo",
                "mol2",

                "-at",
                pcfg[
                    "ligand_forcefield"
                ],

                "-c",
                pcfg[
                    "charge_method"
                ],

                "-nc",
                str(charge),

                "-s",
                "2",
            ]

            ac = run(
                ac_cmd,
                job_dir,
            )

            (
                job_dir /
                "antechamber.stdout"
            ).write_text(
                ac.stdout or ""
            )

            (
                job_dir /
                "antechamber.stderr"
            ).write_text(
                ac.stderr or ""
            )

            if (
                ac.returncode != 0
                or
                not ligand_gaff.exists()
            ):
                raise RuntimeError(
                    "antechamber failed"
                )

            pc_cmd = [
                str(parmchk2),

                "-i",
                str(ligand_gaff),

                "-f",
                "mol2",

                "-o",
                str(frcmod),

                "-s",
                pcfg[
                    "ligand_forcefield"
                ],
            ]

            pc = run(
                pc_cmd,
                job_dir,
            )

            (
                job_dir /
                "parmchk2.stdout"
            ).write_text(
                pc.stdout or ""
            )

            (
                job_dir /
                "parmchk2.stderr"
            ).write_text(
                pc.stderr or ""
            )

            if (
                pc.returncode != 0
                or
                not frcmod.exists()
            ):
                raise RuntimeError(
                    "parmchk2 failed"
                )

            leap_text = f"""\
source leaprc.protein.{pcfg['protein_forcefield']}
source leaprc.{pcfg['ligand_forcefield']}

REC = loadPDB {receptor_copy}
LIG = loadMol2 {ligand_gaff}

loadAmberParams {frcmod}

COM = combine {{ REC LIG }}

check COM

saveAmberParm COM {prmtop} {inpcrd}
savePDB COM {complex_pdb}

quit
"""

            leap_input.write_text(
                leap_text
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
                "tleap.stdout"
            ).write_text(
                leap.stdout or ""
            )

            (
                job_dir /
                "tleap.stderr"
            ).write_text(
                leap.stderr or ""
            )

            if (
                leap.returncode != 0
                or
                not prmtop.exists()
                or
                not inpcrd.exists()
            ):
                raise RuntimeError(
                    "tleap failed"
                )

            amber = pmd.load_file(
                str(prmtop),
                str(inpcrd),
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

            if (
                not gmx_top.exists()
                or
                not gmx_gro.exists()
            ):
                raise RuntimeError(
                    "ParmEd GROMACS "
                    "conversion failed"
                )

            metadata = {
                "dynamics_job_id":
                    job_id,

                "target_id":
                    row["target_id"],

                "ligand_id":
                    ligand_id,

                "formal_charge":
                    charge,

                "receptor_preparation": {
                    "strip_crystallographic_waters":
                        strip_waters,

                    "water_residue_names":
                        sorted(
                            water_names
                        ),

                    "removed_water_residues":
                        receptor_qc[
                            "removed_water_residues"
                        ],

                    "removed_water_atoms":
                        receptor_qc[
                            "removed_water_atoms"
                        ],
                },

                "protein_forcefield":
                    pcfg[
                        "protein_forcefield"
                    ],

                "ligand_forcefield":
                    pcfg[
                        "ligand_forcefield"
                    ],

                "charge_method":
                    pcfg[
                        "charge_method"
                    ],

                "antechamber_command":
                    ac_cmd,

                "parmchk2_command":
                    pc_cmd,

                "tleap_command":
                    leap_cmd,

                "files": {
                    "ligand_gaff2":
                        str(ligand_gaff),

                    "frcmod":
                        str(frcmod),

                    "prmtop":
                        str(prmtop),

                    "inpcrd":
                        str(inpcrd),

                    "complex_pdb":
                        str(complex_pdb),

                    "gromacs_top":
                        str(gmx_top),

                    "gromacs_gro":
                        str(gmx_gro),
                },

                "sha256": {
                    "prmtop":
                        sha256(prmtop),

                    "inpcrd":
                        sha256(inpcrd),

                    "gromacs_top":
                        sha256(gmx_top),

                    "gromacs_gro":
                        sha256(gmx_gro),
                },

                "status":
                    "PASS",
            }

            manifest_json.write_text(
                json.dumps(
                    metadata,
                    indent=2,
                )
                + "\n"
            )

            result_rows.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "target_id":
                        row["target_id"],

                    "ligand_id":
                        ligand_id,

                    "system_top":
                        str(gmx_top),

                    "system_gro":
                        str(gmx_gro),

                    "status":
                        "PASS",
                }
            )

            print(
                f"[M11] {job_id} "
                "PASS"
            )

        except Exception as exc:

            failures += 1

            result_rows.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "target_id":
                        row["target_id"],

                    "ligand_id":
                        ligand_id,

                    "system_top":
                        "",

                    "system_gro":
                        "",

                    "status":
                        "FAIL",
                }
            )

            print(
                f"[M11] {job_id} "
                f"FAIL "
                f"{type(exc).__name__}: "
                f"{exc}"
            )

    result_table = (
        output_root /
        "parameterization_manifest.tsv"
    )

    fields = [
        "dynamics_job_id",
        "target_id",
        "ligand_id",
        "system_top",
        "system_gro",
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
            result_rows
        )

    passed = sum(
        r["status"] == "PASS"
        for r in result_rows
    )

    if failures:

        status = (
            "FAIL_M11_PARAMETERIZATION"
        )

    elif passed == 0:

        status = (
            "SKIP_M11_NO_VALID_SYSTEMS"
        )

    else:

        status = (
            "PASS_M11_PARAMETERIZATION"
        )

    log_path.write_text(
        "\n".join(
            [
                status,
                (
                    "Systems requested: "
                    f"{len(rows)}"
                ),
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
        "Systems requested:",
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
