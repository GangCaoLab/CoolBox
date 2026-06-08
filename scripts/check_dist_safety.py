#!/usr/bin/env python3
"""Reject distribution artifacts with unexpected executable payload hooks."""

import argparse
import sys
import tarfile
import zipfile
from pathlib import Path, PurePosixPath


SUSPICIOUS_SUFFIXES = (".pth", ".js", ".mjs", ".cjs")
SUSPICIOUS_NAMES = {
    ".npmrc",
    "bun.lock",
    "bun.lockb",
    "package.json",
    "package-lock.json",
    "pnpm-lock.yaml",
    "yarn.lock",
}
SUSPICIOUS_PARTS = {"node_modules", "__pycache__"}
SUSPICIOUS_CONTENT = (
    "oven-sh/bun",
    "bun/releases/download",
    "bun-v1.",
)


def normalize_member_name(name):
    return name.replace("\\", "/")


def is_path_unsafe(name):
    path = PurePosixPath(normalize_member_name(name))
    return path.is_absolute() or ".." in path.parts


def is_suspicious_name(name):
    path = PurePosixPath(normalize_member_name(name))
    lowered = path.name.lower()
    suffix = path.suffix.lower()

    return (
        suffix in SUSPICIOUS_SUFFIXES
        or lowered in SUSPICIOUS_NAMES
        or any(part.lower() in SUSPICIOUS_PARTS for part in path.parts)
    )


def content_matches(content):
    try:
        text = content.decode("utf-8", errors="ignore")
    except AttributeError:
        return False
    return any(pattern in text for pattern in SUSPICIOUS_CONTENT)


def scan_zip(path):
    findings = []
    with zipfile.ZipFile(path) as archive:
        for member in archive.infolist():
            if member.is_dir():
                continue
            if is_path_unsafe(member.filename):
                findings.append((member.filename, "unsafe archive path"))
                continue
            if is_suspicious_name(member.filename):
                findings.append((member.filename, "unexpected executable or JS-related file"))
                continue
            with archive.open(member) as fh:
                if content_matches(fh.read()):
                    findings.append((member.filename, "suspicious payload marker in file content"))
    return findings


def scan_tar(path):
    findings = []
    with tarfile.open(path) as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            if is_path_unsafe(member.name):
                findings.append((member.name, "unsafe archive path"))
                continue
            if is_suspicious_name(member.name):
                findings.append((member.name, "unexpected executable or JS-related file"))
                continue
            extracted = archive.extractfile(member)
            if extracted and content_matches(extracted.read()):
                findings.append((member.name, "suspicious payload marker in file content"))
    return findings


def scan_artifact(path):
    if zipfile.is_zipfile(path):
        return scan_zip(path)
    if tarfile.is_tarfile(path):
        return scan_tar(path)
    return [(path.name, "unsupported distribution archive format")]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args(argv)

    all_findings = []
    for artifact in args.artifacts:
        if not artifact.exists():
            all_findings.append((str(artifact), "artifact does not exist"))
            continue
        for member, reason in scan_artifact(artifact):
            all_findings.append((f"{artifact}: {member}", reason))

    if all_findings:
        print("Unsafe distribution contents detected:")
        for target, reason in all_findings:
            print(f"- {target}: {reason}")
        return 1

    print("Distribution safety check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
