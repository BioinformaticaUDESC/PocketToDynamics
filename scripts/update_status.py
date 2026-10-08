#!/usr/bin/env python3

from pathlib import Path
import re


ROOT = Path.home() / "PocketToDynamics"
LOG_DIR = ROOT / "logs"
OUT = ROOT / "manifest" / "STATUS.tsv"


SPECIAL = {
    "M05_engine_preflight.log": "M05A",
    "M05_prepare_inputs.log": "M05B",
    "M05C_screening.log": "M05C",
}


def module_id(path: Path):

    if path.name in SPECIAL:
        return SPECIAL[path.name]

    match = re.match(
        r"(M\d+)",
        path.name,
    )

    if match:
        return match.group(1)

    return None


rows = []

for path in sorted(
    LOG_DIR.glob("M*.log")
):

    text = path.read_text(
        errors="replace"
    ).strip()

    if not text:
        continue

    first = (
        text.splitlines()[0].strip()
    )

    module = module_id(path)

    if not module:
        continue

    if first.startswith("PASS_"):
        status = "PASS"

    elif first.startswith("WARN_"):
        status = "WARN"

    elif first.startswith("FAIL_"):
        status = "FAIL"

    elif first.startswith("SKIP_"):
        status = "SKIP"

    else:
        status = "UNKNOWN"

    rows.append(
        (
            module,
            status,
            first,
            str(
                path.relative_to(ROOT)
            ),
        )
    )


def sort_key(row):

    module = row[0]

    m = re.match(
        r"M(\d+)([A-Z]?)$",
        module,
    )

    if not m:
        return (999, "")

    return (
        int(m.group(1)),
        m.group(2),
    )


rows.sort(
    key=sort_key
)


OUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)


with OUT.open("w") as f:

    f.write(
        "module\tstatus\tmessage\tlog_file\n"
    )

    for row in rows:

        f.write(
            "\t".join(row)
            + "\n"
        )


summary = {}

for _, status, _, _ in rows:

    summary[status] = (
        summary.get(status, 0)
        + 1
    )


print(
    "===== PIPELINE STATUS ====="
)

for (
    module,
    status,
    message,
    _,
) in rows:

    print(
        f"{module:<5} "
        f"{status:<8} "
        f"{message}"
    )


print()
print(
    "Summary:",
    summary,
)

print(
    "Status file:",
    OUT,
)
