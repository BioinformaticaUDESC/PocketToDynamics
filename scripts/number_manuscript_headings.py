#!/usr/bin/env python3
"""Create a numbered-section copy of the PocketToDynamics Markdown manuscript."""
from pathlib import Path
import argparse
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
source = args.input.read_text().splitlines()
section = 0
subsection = 0
result = []
for line in source:
    match = re.match(r'^(#{2,3}) (.+)$', line)
    if not match:
        result.append(line)
        continue
    level, title = len(match.group(1)), match.group(2)
    if level == 2:
        if title.strip().lower() in {'abstract', 'references'}:
            result.append(line)
            continue
        section += 1
        subsection = 0
        result.append(f'{match.group(1)} {section}. {title}')
    else:
        subsection += 1
        result.append(f'{match.group(1)} {section}.{subsection}. {title}')
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text('\n'.join(result) + '\n')
print(args.output)
