#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path.home() / "PocketToDynamics"

CONFIG = (
    ROOT /
    "config/dynamics.yaml"
)

DEFAULT_INPUT = Path(
    "/mnt/labio-fast-disk/"
    "PocketToDynamics/dynamics/"
    "equilibrated/"
    "equilibration_manifest.tsv"
)

DEFAULT_LOG = (
    ROOT /
    "logs/M14_production_md.log"
)


def sha256(
    path: Path,
) -> str:

    h = hashlib.sha256()

    with path.open("rb") as f:

        for chunk in iter(
            lambda: f.read(
                1024 * 1024
            ),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def load_tsv(
    path: Path,
):

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


def validate_positive(
    name: str,
    value: float,
):

    if (
        not math.isfinite(value)
        or
        value <= 0
    ):

        raise ValueError(
            f"{name} must be "
            f"positive: {value}"
        )


def build_mdp(
    path: Path,
    cfg: dict,
    segment_steps: int,
    init_step: int,
    nstxout_compressed: int,
    nstenergy: int,
    nstlog: int,
):

    dt = float(
        cfg["timestep_ps"]
    )

    path.write_text(
        f"""\
integrator               = md
dt                       = {dt}
nsteps                   = {segment_steps}

tinit                    = 0
init-step                = {init_step}

continuation             = yes

constraint_algorithm     = lincs
constraints              = {cfg['constraints']}

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
coulombtype              = PME
vdwtype                  = Cut-off

tcoupl                   = {cfg['thermostat']}
tc-grps                  = System
tau_t                    = {float(cfg['tau_t_ps'])}
ref_t                    = {float(cfg['temperature_k'])}

pcoupl                   = {cfg['barostat']}
pcoupltype               = isotropic
tau_p                    = {float(cfg['tau_p_ps'])}
ref_p                    = {float(cfg['pressure_bar'])}
compressibility          = {float(cfg['compressibility_bar_inv'])}

gen_vel                  = no

nstxout                   = 0
nstvout                   = 0
nstfout                   = 0
nstxout-compressed        = {nstxout_compressed}

nstenergy                = {nstenergy}
nstlog                   = {nstlog}

pbc                      = xyz
"""
    )


def segment_is_complete(
    segment_dir: Path,
):

    marker = (
        segment_dir /
        "segment_complete.json"
    )

    if not marker.exists():
        return False

    try:

        data = json.loads(
            marker.read_text()
        )

    except Exception:
        return False

    required = [
        segment_dir / "md.gro",
        segment_dir / "md.cpt",
        segment_dir / "md.edr",
        segment_dir / "md.log",
        segment_dir / "md.tpr",
        segment_dir / "md.xtc",
    ]

    if not all(
        p.exists()
        and
        p.stat().st_size > 0
        for p in required
    ):
        return False

    if data.get('status') != 'PASS':
        return False
    expected = data.get('sha256', {})
    return all(expected.get(ext) == sha256(segment_dir / ('md.' + ext))
               for ext in ['gro', 'cpt', 'edr', 'xtc', 'tpr', 'log'])



def remove_partial_outputs(
    segment_dir: Path,
):

    names = [
        "md.gro",
        "md.cpt",
        "md.edr",
        "md.log",
        "md.xtc",
        "md.trr",
        "md_prev.cpt",
        "segment_complete.json",
    ]

    for name in names:

        path = (
            segment_dir /
            name
        )

        if path.exists():
            path.unlink()


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

    parser.add_argument(
        "--profile",
        choices=[
            "production",
            "fixture",
        ],
        default="production",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    args = parser.parse_args()

    cfg = yaml.safe_load(
        CONFIG.read_text()
    )

    pcfg = dict(
        cfg["production"]
    )

    if args.profile == "fixture":

        fixture_cfg = cfg.get(
            "production_fixture",
            {},
        )

        pcfg.update(
            fixture_cfg
        )

    input_manifest = (
        args.manifest
        .expanduser()
        .resolve()
    )

    rows = load_tsv(
        input_manifest
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
            pcfg[
                "heavy_output_root"
            ]
        )

    output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not rows:

        status = (
            "SKIP_M14_NO_EQUILIBRATED_SYSTEMS"
        )

        log_path.write_text(
            status
            + "\nSystems: 0\n"
        )

        print(status)
        print("Systems: 0")

        return

    gmx = Path(
        pcfg[
            "gmx_executable"
        ]
    )

    if not gmx.exists():

        raise SystemExit(
            "FAIL_M14_GROMACS_NOT_FOUND"
        )

    dt = float(
        pcfg[
            "timestep_ps"
        ]
    )

    total_ps = float(
        pcfg[
            "total_duration_ps"
        ]
    )

    segment_ps = float(
        pcfg[
            "segment_duration_ps"
        ]
    )

    validate_positive(
        "timestep_ps",
        dt,
    )

    validate_positive(
        "total_duration_ps",
        total_ps,
    )

    validate_positive(
        "segment_duration_ps",
        segment_ps,
    )

    if (
        segment_ps
        >
        total_ps
    ):

        raise SystemExit(
            "FAIL_M14_SEGMENT_LONGER_THAN_TOTAL"
        )

    total_steps = int(
        round(
            total_ps /
            dt
        )
    )

    nominal_segment_steps = int(
        round(
            segment_ps /
            dt
        )
    )

    if total_steps < 1:

        raise SystemExit(
            "FAIL_M14_INVALID_TOTAL_STEPS"
        )

    if nominal_segment_steps < 1:

        raise SystemExit(
            "FAIL_M14_INVALID_SEGMENT_STEPS"
        )

    ntmpi = int(
        pcfg["ntmpi"]
    )

    ntomp = int(
        pcfg["ntomp"]
    )

    env = os.environ.copy()

    env[
        "OMP_NUM_THREADS"
    ] = str(
        ntomp
    )

    results = []
    failures = 0

    for row in rows:

        if (
            row.get("status")
            != "PASS"
        ):
            continue

        job_id = (
            row[
                "dynamics_job_id"
            ]
        )

        npt_gro = Path(
            row[
                "npt_gro"
            ]
        )

        source_job = (
            npt_gro.parent
        )

        npt_cpt = (
            source_job /
            "npt.cpt"
        )

        source_top = (
            source_job /
            "system.top"
        )

        required_source = [
            npt_gro,
            npt_cpt,
            source_top,
        ]

        missing = [
            str(p)
            for p in required_source
            if (
                not p.exists()
                or
                p.stat().st_size == 0
            )
        ]

        if missing:

            failures += 1

            print(
                f"[M14] {job_id} FAIL "
                "missing equilibration files"
            )

            results.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "profile":
                        args.profile,

                    "target_duration_ps":
                        total_ps,

                    "completed_duration_ps":
                        0.0,

                    "segments_total":
                        0,

                    "segments_completed":
                        0,

                    "final_gro":
                        "",

                    "final_cpt":
                        "",

                    "status":
                        "FAIL",
                }
            )

            continue

        job_dir = (
            output_root /
            job_id
        )

        if (
            args.overwrite
            and
            job_dir.exists()
        ):

            shutil.rmtree(
                job_dir
            )

        job_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        top_copy = (
            job_dir /
            "system.top"
        )

        if (
            not top_copy.exists()
            or
            args.overwrite
        ):

            shutil.copy2(
                source_top,
                top_copy,
            )

        remaining_steps = (
            total_steps
        )

        segment_specs = []

        init_step = 0
        index = 1

        while remaining_steps > 0:

            seg_steps = min(
                nominal_segment_steps,
                remaining_steps,
            )

            segment_specs.append(
                {
                    "index":
                        index,

                    "steps":
                        seg_steps,

                    "init_step":
                        init_step,

                    "duration_ps":
                        seg_steps * dt,

                    "start_ps":
                        init_step * dt,

                    "end_ps":
                        (
                            init_step
                            +
                            seg_steps
                        )
                        *
                        dt,
                }
            )

            init_step += (
                seg_steps
            )

            remaining_steps -= (
                seg_steps
            )

            index += 1

        segment_records = []

        previous_gro = (
            npt_gro
        )

        previous_cpt = (
            npt_cpt
        )

        job_failed = False

        try:
            provenance = {
                'config': pcfg, 'profile': args.profile,
                'inputs': {str(p): sha256(p) for p in required_source},
                'engine_sha256': sha256(gmx.resolve()),
            }
            signature = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
            provenance_path = job_dir / 'production_inputs.json'
            if provenance_path.exists():
                if json.loads(provenance_path.read_text()).get('signature') != signature:
                    raise RuntimeError('Production inputs/config changed; use a new output root or explicit --overwrite')
            elif any(job_dir.glob('segment_*/md.cpt')):
                raise RuntimeError('Legacy checkpoints lack input provenance; use a new output root to revalidate')
            if sha256(top_copy) != sha256(source_top):
                raise RuntimeError('Production topology differs from source')
            provenance_path.write_text(json.dumps({'signature': signature, 'inputs': provenance}, indent=2)+'\n')

            for spec in segment_specs:

                seg_index = (
                    spec["index"]
                )

                seg_name = (
                    f"segment_"
                    f"{seg_index:04d}"
                )

                segment_dir = (
                    job_dir /
                    seg_name
                )

                segment_dir.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                mdp = (
                    segment_dir /
                    "production.mdp"
                )

                tpr = (
                    segment_dir /
                    "md.tpr"
                )

                prefix = (
                    segment_dir /
                    "md"
                )

                complete = (
                    segment_is_complete(
                        segment_dir
                    )
                )

                marker = segment_dir / 'segment_complete.json'
                if marker.exists() and not complete and not args.overwrite:
                    raise RuntimeError(f'{seg_name}: completed artifact integrity failed; use a new output root')

                if (
                    complete
                    and
                    bool(
                        pcfg.get(
                            "resume",
                            True,
                        )
                    )
                    and
                    not args.overwrite
                ):

                    print(
                        f"[M14] {job_id} "
                        f"{seg_name} "
                        "SKIP_COMPLETED"
                    )

                    previous_gro = (
                        segment_dir /
                        "md.gro"
                    )

                    previous_cpt = (
                        segment_dir /
                        "md.cpt"
                    )

                    segment_records.append(
                        {
                            **spec,

                            "segment":
                                seg_name,

                            "status":
                                "SKIP_COMPLETED",

                            "gro":
                                str(
                                    previous_gro
                                ),

                            "cpt":
                                str(
                                    previous_cpt
                                ),
                        }
                    )

                    continue

                build_mdp(
                    mdp,
                    pcfg,
                    spec["steps"],
                    spec["init_step"],
                    int(
                        pcfg[
                            "nstxout_compressed"
                        ]
                    ),
                    int(
                        pcfg[
                            "nstenergy"
                        ]
                    ),
                    int(
                        pcfg[
                            "nstlog"
                        ]
                    ),
                )

                grompp_cmd = [
                    str(gmx),
                    "grompp",

                    "-f",
                    str(mdp),

                    "-c",
                    str(
                        previous_gro
                    ),

                    "-t",
                    str(
                        previous_cpt
                    ),

                    "-p",
                    str(
                        top_copy
                    ),

                    "-o",
                    str(tpr),

                    "-maxwarn",
                    "0",
                ]

                grompp = run(
                    grompp_cmd,
                    segment_dir,
                )

                (
                    segment_dir /
                    "grompp.stdout"
                ).write_text(
                    grompp.stdout
                    or ""
                )

                (
                    segment_dir /
                    "grompp.stderr"
                ).write_text(
                    grompp.stderr
                    or ""
                )

                if (
                    grompp.returncode
                    != 0
                    or
                    not tpr.exists()
                    or
                    tpr.stat().st_size
                    == 0
                ):

                    raise RuntimeError(
                        f"{seg_name} "
                        "grompp failed"
                    )

                checkpoint = (
                    segment_dir /
                    "md.cpt"
                )

                partial_checkpoint = (
                    checkpoint.exists()
                    and
                    checkpoint.stat().st_size
                    > 0
                )

                if (
                    partial_checkpoint
                    and
                    bool(
                        pcfg.get(
                            "resume",
                            True,
                        )
                    )
                    and
                    not args.overwrite
                ):

                    mode = "RESUME"

                    mdrun_cmd = [
                        str(gmx),
                        "mdrun",

                        "-deffnm",
                        str(prefix),

                        "-s",
                        str(tpr),

                        "-cpi",
                        str(
                            checkpoint
                        ),

                        "-append",

                        "-ntmpi",
                        str(ntmpi),

                        "-ntomp",
                        str(ntomp),
                    ]

                else:

                    mode = "START"

                    remove_partial_outputs(
                        segment_dir
                    )

                    mdrun_cmd = [
                        str(gmx),
                        "mdrun",

                        "-deffnm",
                        str(prefix),

                        "-s",
                        str(tpr),

                        "-ntmpi",
                        str(ntmpi),

                        "-ntomp",
                        str(ntomp),
                    ]

                print(
                    f"[M14] {job_id} "
                    f"{seg_name} "
                    f"{mode} "
                    f"{spec['start_ps']:.3f}"
                    "→"
                    f"{spec['end_ps']:.3f} ps"
                )

                mdrun = run(
                    mdrun_cmd,
                    segment_dir,
                    env=env,
                )

                (
                    segment_dir /
                    "mdrun.stdout"
                ).write_text(
                    mdrun.stdout
                    or ""
                )

                (
                    segment_dir /
                    "mdrun.stderr"
                ).write_text(
                    mdrun.stderr
                    or ""
                )

                gro = (
                    segment_dir /
                    "md.gro"
                )

                cpt = (
                    segment_dir /
                    "md.cpt"
                )

                edr = (
                    segment_dir /
                    "md.edr"
                )

                log = (
                    segment_dir /
                    "md.log"
                )

                xtc = (
                    segment_dir /
                    "md.xtc"
                )

                required_outputs = [
                    gro,
                    cpt,
                    edr,
                    log,
                    xtc,
                ]

                if (
                    mdrun.returncode
                    != 0
                    or
                    not all(
                        p.exists()
                        and
                        p.stat().st_size
                        > 0
                        for p
                        in required_outputs
                    )
                ):

                    raise RuntimeError(
                        f"{seg_name} "
                        "mdrun failed or "
                        "outputs incomplete"
                    )

                completion = {
                    "dynamics_job_id":
                        job_id,

                    "segment":
                        seg_name,

                    "segment_index":
                        seg_index,

                    "start_ps":
                        spec[
                            "start_ps"
                        ],

                    "end_ps":
                        spec[
                            "end_ps"
                        ],

                    "duration_ps":
                        spec[
                            "duration_ps"
                        ],

                    "steps":
                        spec[
                            "steps"
                        ],

                    "init_step":
                        spec[
                            "init_step"
                        ],

                    "mode":
                        mode,

                    "sha256": {
                        "tpr": sha256(tpr),
                        "log": sha256(segment_dir / "md.log"),
                        "gro":
                            sha256(gro),

                        "cpt":
                            sha256(cpt),

                        "edr":
                            sha256(edr),

                        "xtc":
                            sha256(xtc),
                    },

                    "status":
                        "PASS",
                }

                (
                    segment_dir /
                    "segment_complete.json"
                ).write_text(
                    json.dumps(
                        completion,
                        indent=2,
                    )
                    + "\n"
                )

                segment_records.append(
                    {
                        **spec,

                        "segment":
                            seg_name,

                        "status":
                            "PASS",

                        "gro":
                            str(gro),

                        "cpt":
                            str(cpt),
                    }
                )

                previous_gro = gro
                previous_cpt = cpt

            completed_ps = sum(
                float(
                    x[
                        "duration_ps"
                    ]
                )
                for x
                in segment_records
                if (
                    x[
                        "status"
                    ]
                    in {
                        "PASS",
                        "SKIP_COMPLETED",
                    }
                )
            )

            if not math.isclose(
                completed_ps,
                total_ps,
                rel_tol=0.0,
                abs_tol=max(
                    dt,
                    1e-6,
                ),
            ):

                raise RuntimeError(
                    "Accumulated production "
                    "duration mismatch: "
                    f"{completed_ps} "
                    f"!= {total_ps}"
                )

            final_gro = (
                previous_gro
            )

            final_cpt = (
                previous_cpt
            )

            production_manifest = {
                "dynamics_job_id":
                    job_id,

                "profile":
                    args.profile,

                "timestep_ps":
                    dt,

                "target_duration_ps":
                    total_ps,

                "segment_duration_ps":
                    segment_ps,

                "total_steps":
                    total_steps,

                "segments_total":
                    len(
                        segment_specs
                    ),

                "segments_completed":
                    len(
                        segment_records
                    ),

                "completed_duration_ps":
                    completed_ps,

                "temperature_k":
                    float(
                        pcfg[
                            "temperature_k"
                        ]
                    ),

                "pressure_bar":
                    float(
                        pcfg[
                            "pressure_bar"
                        ]
                    ),

                "parallelization": {
                    "ntmpi":
                        ntmpi,

                    "ntomp":
                        ntomp,

                    "omp_num_threads":
                        ntomp,
                },

                "resume_enabled":
                    bool(
                        pcfg.get(
                            "resume",
                            True,
                        )
                    ),

                "segments":
                    segment_records,

                "files": {
                    "topology":
                        str(
                            top_copy
                        ),

                    "final_gro":
                        str(
                            final_gro
                        ),

                    "final_cpt":
                        str(
                            final_cpt
                        ),
                },

                "sha256": {
                    "final_gro":
                        sha256(
                            final_gro
                        ),

                    "final_cpt":
                        sha256(
                            final_cpt
                        ),
                },

                "status":
                    "PASS",
            }

            (
                job_dir /
                "production_manifest.json"
            ).write_text(
                json.dumps(
                    production_manifest,
                    indent=2,
                )
                + "\n"
            )

            results.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "profile":
                        args.profile,

                    "target_duration_ps":
                        total_ps,

                    "completed_duration_ps":
                        completed_ps,

                    "segments_total":
                        len(
                            segment_specs
                        ),

                    "segments_completed":
                        len(
                            segment_records
                        ),

                    "final_gro":
                        str(
                            final_gro
                        ),

                    "final_cpt":
                        str(
                            final_cpt
                        ),

                    "status":
                        "PASS",
                }
            )

            print(
                f"[M14] {job_id} PASS "
                f"{completed_ps:.3f}/"
                f"{total_ps:.3f} ps "
                f"segments="
                f"{len(segment_records)}"
            )

        except Exception as exc:

            job_failed = True
            failures += 1

            results.append(
                {
                    "dynamics_job_id":
                        job_id,

                    "profile":
                        args.profile,

                    "target_duration_ps":
                        total_ps,

                    "completed_duration_ps":
                        sum(
                            float(
                                x.get(
                                    "duration_ps",
                                    0.0,
                                )
                            )
                            for x
                            in segment_records
                        ),

                    "segments_total":
                        len(
                            segment_specs
                        ),

                    "segments_completed":
                        len(
                            segment_records
                        ),

                    "final_gro":
                        "",

                    "final_cpt":
                        "",

                    "status":
                        "FAIL",
                }
            )

            print(
                f"[M14] {job_id} FAIL "
                f"{type(exc).__name__}: "
                f"{exc}"
            )

    table = (
        output_root /
        "production_manifest.tsv"
    )

    fields = [
        "dynamics_job_id",
        "profile",
        "target_duration_ps",
        "completed_duration_ps",
        "segments_total",
        "segments_completed",
        "final_gro",
        "final_cpt",
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
            "FAIL_M14_PRODUCTION_MD"
        )

    elif passed == 0:

        status = (
            "SKIP_M14_NO_VALID_SYSTEMS"
        )

    else:

        status = (
            "PASS_M14_PRODUCTION_MD"
        )

    log_path.write_text(
        "\n".join(
            [
                status,
                (
                    f"Profile: "
                    f"{args.profile}"
                ),
                (
                    f"Systems: "
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
