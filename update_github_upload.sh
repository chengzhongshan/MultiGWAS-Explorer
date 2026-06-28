#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

COMMIT_MSG="${1:-Update MultiGWAS-Explorer scripts}"
REMOTE_NAME="${REMOTE_NAME:-origin}"
BRANCH_NAME="${BRANCH_NAME:-main}"

echo "Working directory: $ROOT_DIR"
echo "Remote: $REMOTE_NAME"
echo "Branch: $BRANCH_NAME"
echo "Commit message: $COMMIT_MSG"
echo

git status
echo
git diff --stat || true
echo

git add .

if git diff --cached --quiet; then
  echo "No staged changes to commit."
else
  git commit -m "$COMMIT_MSG"
fi

git pull --rebase "$REMOTE_NAME" "$BRANCH_NAME"
git push "$REMOTE_NAME" "$BRANCH_NAME"

echo
echo "GitHub upload package is up to date."
