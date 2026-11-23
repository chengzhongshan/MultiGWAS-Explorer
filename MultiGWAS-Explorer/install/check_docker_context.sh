#!/usr/bin/env bash
# Fast COPY-only build: prove the actual Docker context contains every file
# required by the pinned SAS runtime, before installing system dependencies.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "${script_dir}/.." && pwd)"
dockerfile="$(mktemp)"
trap 'rm -f "$dockerfile"' EXIT
python3 - "$root" > "$dockerfile" <<'PY'
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
manifest=json.loads((root/'docs/SAS_RUNTIME_MACRO_QUERY_MANIFEST.json').read_text())
print('FROM scratch')
for item in manifest['files']:
    path=Path(item['path'])
    if path.is_absolute() or '..' in path.parts:
        raise SystemExit('Unsafe manifest source path')
    print('COPY '+json.dumps([path.as_posix(),'/runtime/'+path.as_posix()]))
PY
docker build --progress=plain --file "$dockerfile" "$root"
printf 'Docker runtime source context: PASS\n'
