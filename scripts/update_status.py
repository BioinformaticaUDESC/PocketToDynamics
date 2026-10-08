#!/usr/bin/env python3

from pathlib import Path
import re

ROOT = Path.home() / "PocketToDynamics"
LOG_DIR = ROOT / "logs"
OUT = ROOT / "manifest" / "STATUS.tsv"

rows = []

for path in sorted(LOG_DIR.glob("M*.log")):

    text = path.read_text(errors="replace").strip()

    if not text:
        continue

    first = text.splitlines()[0].strip()

    m = re.match(r"(M\d+)", path.name)

    if not m:
        continue

    module = m.group(1)

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
            str(path.relative_to(ROOT)),
        )
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
        f.write("\t".join(row) + "\n")

summary = {}

for _, status, _, _ in rows:
    summary[status] = summary.get(status, 0) + 1

print("===== PIPELINE STATUS =====")

for module, status, message, _ in rows:
    print(
        f"{module:<5} "
        f"{status:<8} "
        f"{message}"
    )

print()
print("Summary:", summary)
print("Status file:", OUT)
