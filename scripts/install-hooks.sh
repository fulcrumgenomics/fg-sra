#!/bin/bash
#
# Install git hooks for fg-sra development.
#
# Usage: ./scripts/install-hooks.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
HOOKS_SOURCE="${SCRIPT_DIR}/hooks"
HOOKS_DEST="${REPO_ROOT}/.git/hooks"

echo "Installing git hooks for fg-sra..."

if [[ ! -f "${REPO_ROOT}/Cargo.toml" ]]; then
    echo "Error: Must be run from the fg-sra repository"
    exit 1
fi

if [[ ! -d "${REPO_ROOT}/.git" ]]; then
    echo "Error: .git directory not found"
    exit 1
fi

mkdir -p "${HOOKS_DEST}"

for hook in "${HOOKS_SOURCE}"/*; do
    if [[ -f "${hook}" ]]; then
        hook_name="$(basename "${hook}")"
        dest="${HOOKS_DEST}/${hook_name}"

        if [[ -e "${dest}" ]]; then
            if [[ -L "${dest}" ]]; then
                echo "Updating symlink: ${hook_name}"
                rm "${dest}"
            else
                echo "Backing up existing ${hook_name} to ${hook_name}.backup"
                mv "${dest}" "${dest}.backup"
            fi
        else
            echo "Installing: ${hook_name}"
        fi

        ln -s "${hook}" "${dest}"
        chmod +x "${hook}"
    fi
done

echo ""
echo "Git hooks installed successfully!"
echo ""
echo "The following hooks are now active:"
echo "  - pre-commit: Runs cargo ci-fmt and ci-lint before each commit"
