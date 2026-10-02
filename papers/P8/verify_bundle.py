"""Read-only closed-inventory and SHA-256 check, not a proof or signature."""
import hashlib
import json
from pathlib import Path, PurePosixPath


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root/'MANIFEST.json').read_bytes())
    rows = manifest['members']
    names = [r['path'] for r in rows]
    if names != sorted(set(names)) or 'MANIFEST.json' in names:
        raise ValueError('MANIFEST_INVENTORY')
    actual = sorted(p.relative_to(root).as_posix() for p in root.rglob('*')
                    if p.is_file() and p != root/'MANIFEST.json')
    if actual != names:
        raise ValueError('EXTRA_OR_MISSING_FILES')
    for row in rows:
        name = row['path']
        rel = PurePosixPath(name)
        path = root/name
        if rel.is_absolute() or '..' in rel.parts or '\\' in name or path.is_symlink():
            raise ValueError('UNSAFE_PATH')
        if not path.resolve().is_relative_to(root):
            raise ValueError('PATH_ESCAPE')
        data = path.read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ValueError('MEMBER_MISMATCH:'+name)
    print('MATH_PAPER_PACKAGE_INTEGRITY_PASS', manifest['paper'], len(rows))


if __name__ == '__main__':
    main()
