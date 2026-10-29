#!/usr/bin/env bash
set -euo pipefail

CONFIG="${ISO_SPATIAL_CONFIG:-config/run.yaml}"
CORES="${CORES:-1}"
ENV_PREFIX="${ISO_SPATIAL_ENV:-}"
RUNTIME_ROOT="${ISO_SPATIAL_RUNTIME_ROOT:-tmp/snakemake_runtime}"
DRY_RUN=0
TARGETS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --target) TARGETS+=("$2"); shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

mkdir -p "$RUNTIME_ROOT/tmp" "$RUNTIME_ROOT/xdg-cache"
export TMPDIR="${ISO_SPATIAL_TMPDIR:-$RUNTIME_ROOT/tmp}"
export TMP="$TMPDIR" TEMP="$TMPDIR"
export XDG_CACHE_HOME="${ISO_SPATIAL_CACHE:-$RUNTIME_ROOT/xdg-cache}"
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME"

ARGS=(--snakefile pipeline/Snakefile --configfile "$CONFIG" --cores "$CORES" --rerun-incomplete --default-resources "tmpdir='$TMPDIR'")
if [[ "$DRY_RUN" == 1 ]]; then ARGS+=(--dry-run); fi
if [[ "${#TARGETS[@]}" -gt 0 ]]; then ARGS+=("${TARGETS[@]}"); fi

if [[ -n "$ENV_PREFIX" && -x "$ENV_PREFIX/bin/snakemake" ]]; then
  exec "$ENV_PREFIX/bin/snakemake" "${ARGS[@]}"
fi
exec snakemake "${ARGS[@]}"
