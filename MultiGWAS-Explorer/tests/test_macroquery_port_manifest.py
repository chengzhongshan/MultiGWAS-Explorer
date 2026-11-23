"""Verify the incremental port map without SAS, credentials or network."""
import hashlib
import json
from pathlib import Path

root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'docs/SAS_RUNTIME_MACRO_QUERY_MANIFEST.json').read_text())
assert manifest['upstream_commit']=='58c2fd4378b27b65e9d8d2c684d3b0c208790a66'
for item in manifest['files']:
    file=root/item['path']
    assert root in file.resolve().parents, item['path']
    assert hashlib.sha256(file.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==item['sha256'],item['path']
print('PASS: incremental source manifest:',len(manifest['files']),'files.')
