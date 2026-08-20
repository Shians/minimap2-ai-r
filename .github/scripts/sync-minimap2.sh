#!/usr/bin/env bash
# Sync src/minimap2/ to a given upstream lh3/minimap2 release tag.
# Usage: sync-minimap2.sh vX.Y
set -euo pipefail

tag="${1:?usage: sync-minimap2.sh <tag>}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT

echo "Downloading minimap2 ${tag}..."
curl -sSL "https://github.com/lh3/minimap2/archive/refs/tags/${tag}.tar.gz" -o "$work_dir/minimap2.tar.gz"
tar -xzf "$work_dir/minimap2.tar.gz" -C "$work_dir"

src_dir="$(find "$work_dir" -mindepth 1 -maxdepth 1 -type d -name 'minimap2-*')"
if [ -z "$src_dir" ]; then
    echo "Could not find extracted minimap2 source directory" >&2
    exit 1
fi

# Mirror the upstream tree into src/minimap2/, dropping the pieces the R
# package doesn't need to build (python bindings, upstream's own CI/test
# suite, packaging metadata for PyPI, docs sources).
rsync -a --delete \
    --exclude='.git' \
    --exclude='.github' \
    --exclude='.gitmodules' \
    --exclude='Makefile.simde' \
    --exclude='MANIFEST.in' \
    --exclude='pyproject.toml' \
    --exclude='setup.py' \
    --exclude='code_of_conduct.md' \
    --exclude='python/' \
    --exclude='test/' \
    --exclude='tex/' \
    --exclude='misc/' \
    "$src_dir"/ "$repo_root/src/minimap2/"

echo "$tag" > "$repo_root/.minimap2-version"

version="${tag#v}"
sed -E -i.bak "s/^Version: .*/Version: ${version}/" "$repo_root/DESCRIPTION"
rm -f "$repo_root/DESCRIPTION.bak"

echo "Synced src/minimap2/ to ${tag} (package Version set to ${version})"
