#!/usr/bin/env bash
set -euo pipefail

# Run from anywhere: bash update_github_upload.sh "Describe the change"
# Preview without changing the real Git index: bash update_github_upload.sh --dry-run
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$ROOT_DIR"

REMOTE_NAME="${REMOTE_NAME:-origin}"
BRANCH_NAME="${BRANCH_NAME:-main}"
COMMIT_MSG='Update MultiGWAS-Explorer scripts'
DRY_RUN=0
MESSAGE_SET=0

usage() {
  cat <<'USAGE'
Usage: bash update_github_upload.sh [--dry-run] ["commit message"]

Stages source and documentation changes, excluding SAS ODA and gnuplot run
artifacts. --dry-run previews the staged files without changing Git or GitHub.
REMOTE_NAME and BRANCH_NAME can override the default origin/main destination.
USAGE
}

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    -h|--help) usage; exit 0 ;;
    --*) echo "Unknown option: $arg" >&2; usage >&2; exit 2 ;;
    *)
      if (( MESSAGE_SET )); then
        echo 'Provide the commit message as one quoted argument.' >&2
        exit 2
      fi
      COMMIT_MSG="$arg"
      MESSAGE_SET=1
      ;;
  esac
done

[[ -n "$COMMIT_MSG" ]] || { echo 'Commit message cannot be empty.' >&2; exit 2; }
git rev-parse --is-inside-work-tree >/dev/null
[[ -z "$(git rev-parse --show-prefix)" ]] || {
  echo "The script must live at the Git repository root: $ROOT_DIR" >&2
  exit 2
}
current_branch="$(git symbolic-ref --quiet --short HEAD)" || {
  echo 'Cannot upload from a detached HEAD.' >&2
  exit 2
}
[[ "$current_branch" == "$BRANCH_NAME" ]] || {
  echo "Current branch is $current_branch; expected $BRANCH_NAME. Set BRANCH_NAME explicitly if intended." >&2
  exit 2
}
git remote get-url "$REMOTE_NAME" >/dev/null

is_generated_file() {
  local path="$1"
  # Tracked benchmark fixtures are deliberate source material.
  [[ "$path" == MultiGWAS-Explorer/benchmark/* ]] && return 1
  case "$path" in
    MultiGWAS-Explorer/cache/*|MultiGWAS-Explorer/local/*|\
    MultiGWAS-Explorer/run_local*/*|MultiGWAS-Explorer/run_manhattan_*/*|\
    MultiGWAS-Explorer/run_single_snp_with_gtf_*/*|MultiGWAS-Explorer/upload_*/*|\
    MultiGWAS-Explorer/tmp*/*|MultiGWAS-Explorer/debug_single_local_gtf_*/*|\
    MultiGWAS-Explorer/configs/auto_*_diff_merged_*|\
    MultiGWAS-Explorer/auto_gtf_import_single_snp.*.sas|\
    MultiGWAS-Explorer/auto_wide_import_single_snp.*.sas|\
    MultiGWAS-Explorer/run_sas_oda_single_snp_with_gtf.*.sas|\
    MultiGWAS-Explorer/run_sas_oda_local_top_hits_manhattan.*.sas|\
    MultiGWAS-Explorer/run_sas_local_debug_*.sas|\
    MultiGWAS-Explorer/*.local_debug.sas|\
    MultiGWAS-Explorer/*_top_hits_*.csv|\
    MultiGWAS-Explorer/gnuplot_fallback_*.status.json)
      return 0 ;;
  esac
  case "${path##*/}" in
    *.png|*.html|*.gp|*.plot.tsv|*.manifest.tsv|*.status.json|\
    *.info.txt|*.log|*.gz|*.bgz|*.tbi|*.csi|*.sas7bdat|*.request.md5)
      return 0 ;;
  esac
  return 1
}

is_new_source_file() {
  local path="$1"
  is_generated_file "$path" && return 1
  case "${path##*/}" in
    .gitignore|Makefile|*.pl|*.pm|*.sh|*.sas|*.py|*.R|*.Rmd|\
    *.md|*.txt|*.json|*.yaml|*.yml|*.toml|*.ini|*.cfg|*.sql|\
    *.js|*.ts|*.ipynb)
      return 0 ;;
  esac
  return 1
}

if (( DRY_RUN )); then
  index_path="$(git rev-parse --git-path index)"
  preview_index="$(mktemp "${TMPDIR:-/tmp}/multigwas-index.XXXXXXXX")"
  [[ ! -f "$index_path" ]] || cp "$index_path" "$preview_index"
  export GIT_INDEX_FILE="$preview_index"
  trap 'rm -f -- "$preview_index"' EXIT
fi

# Include modifications and deletions of tracked files. New data files and
# run products are never picked up merely because they are in the workspace.
git add -u -- .
while IFS= read -r -d '' path; do
  if is_new_source_file "$path"; then
    git add -- "$path"
  fi
done < <(git ls-files --others --exclude-standard -z)

# Protect against results that were manually staged, or were tracked by an
# older commit. This only changes Git's index; local plot files stay in place.
excluded_count=0
while IFS= read -r -d '' path; do
  if is_generated_file "$path"; then
    git restore --staged -- "$path"
    (( excluded_count += 1 ))
  fi
done < <(git diff --cached --name-only --no-renames -z)

echo "Repository: $ROOT_DIR"
echo "Destination: $REMOTE_NAME/$BRANCH_NAME"
echo "Generated files kept out of commit: $excluded_count"
if git diff --cached --quiet; then
  echo 'No source or documentation changes to upload.'
  exit 0
fi

git diff --cached --check
echo 'Files to commit:'
git diff --cached --name-status
if (( DRY_RUN )); then
  echo 'Preview only; Git index and GitHub were not changed.'
  exit 0
fi

remote_ref="refs/remotes/$REMOTE_NAME/$BRANCH_NAME"
git fetch "$REMOTE_NAME" "refs/heads/$BRANCH_NAME:$remote_ref"
if ! git show-ref --verify --quiet "$remote_ref"; then
  echo "Cannot find $remote_ref after fetch; no commit or push was made." >&2
  exit 1
fi
if ! git merge-base --is-ancestor "$remote_ref" HEAD; then
  echo "The remote branch has changes missing locally. No commit or push was made." >&2
  echo "Review and integrate $REMOTE_NAME/$BRANCH_NAME, then rerun this script." >&2
  exit 1
fi

git commit -m "$COMMIT_MSG"
git push "$REMOTE_NAME" "HEAD:refs/heads/$BRANCH_NAME"
echo "Uploaded source changes to $REMOTE_NAME/$BRANCH_NAME."
