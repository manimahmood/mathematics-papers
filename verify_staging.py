"""Read-only export integrity check. Copyright 2026 Luqman Mahmood.

Software rights reserved. This is not a mathematical proof or digital signature.
"""
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import zipfile


def require(condition, label):
    if not condition:
        raise ValueError(label)


def safe_path(root, name):
    rel = PurePosixPath(name)
    require(not rel.is_absolute() and ".." not in rel.parts and "\\" not in name,
            "UNSAFE_MEMBER")
    path = root / name
    require(path.resolve().is_relative_to(root) and not path.is_symlink(),
            "UNSAFE_MEMBER")
    return path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "STAGING_MANIFEST.json").read_bytes())
    require(manifest["schema"] == "math-papers.github-private-staging.v1", "SCHEMA")
    require(manifest["self_excluded"] is True, "SELF_EXCLUSION")
    rows = manifest["members"]
    names = [row["path"] for row in rows]
    require(names == sorted(set(names)), "MEMBER_ORDER")
    actual = sorted(p.relative_to(root).as_posix() for p in root.rglob("*")
                    if p.is_file() and p.relative_to(root).parts[0] != ".git"
                    and p != root / "STAGING_MANIFEST.json")
    require(actual == names, "INVENTORY_MISMATCH")
    for row in rows:
        data = safe_path(root, row["path"]).read_bytes()
        require(len(data) == row["bytes"] and digest(data) == row["sha256"],
                "MEMBER_MISMATCH:" + row["path"])
    catalogue = json.loads((root / "PAPERS.json").read_bytes())
    expected_ids = [f"P{i}" for i in range(1, 9)] + [
        "SM1", "SM2", "SM3", "SM4", "LENS", "SILVER", "EP"]
    require([p["id"] for p in catalogue["papers"]] == expected_ids,
            "PAPER_INVENTORY")
    for paper in catalogue["papers"]:
        folder = safe_path(root, paper["directory"])
        package = json.loads((folder / "MANIFEST.json").read_bytes())
        require(package["paper"] == paper["id"], "PACKAGE_ID")
        for path, field in [(safe_path(root, paper["pdf"]), "pdf_sha256"),
                            (safe_path(root, paper["archive"]), "archive_sha256"),
                            (folder / "MANIFEST.json", "package_manifest_sha256")]:
            require(digest(path.read_bytes()) == paper[field],
                    "CATALOGUE_DIGEST:" + paper["id"] + "/" + field)
        expected = {r["path"] for r in package["members"]} | {"MANIFEST.json"}
        archive = safe_path(root, paper["archive"])
        with zipfile.ZipFile(archive) as z:
            members = z.namelist()
            require(len(members) == len(set(members)) and set(members) == expected,
                    "ARCHIVE_INVENTORY:" + paper["id"])
            for name in members:
                require(z.read(name) == safe_path(folder, name).read_bytes(),
                        "ARCHIVE_BYTES:" + paper["id"] + "/" + name)
        run = subprocess.run([sys.executable, "-B", str(folder / "verify_bundle.py")],
                             cwd=folder, capture_output=True, text=True, timeout=60)
        expected_token = (f"MATH_PAPER_PACKAGE_INTEGRITY_PASS {paper['id']} "
                          f"{len(package['members'])}")
        require(run.returncode == 0 and not run.stderr and run.stdout.strip() == expected_token,
                "PACKAGE_CHECK:" + paper["id"])
        print(expected_token)
    sums = (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    expected_sums = [f"{r['sha256']}  {r['path']}" for r in rows
                     if r["path"] != "SHA256SUMS"]
    require(sums == expected_sums, "SHA256SUMS_MISMATCH")
    print(f"PRIVATE_PAPERS_RELEASE_STAGING_PASS 15 papers / {len(rows)} bound files")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("PRIVATE_PAPERS_STAGING_REFUSED " + str(error), file=sys.stderr)
        sys.exit(1)
