#!/bin/bash

# Main script to generate Table 5: Monte Carlo Results
#
# Usage:
#   ./table_5.sh [--partial|--full]        Run all remaining experiments
#   ./table_5.sh --experiment CONFIG       Run a specific experiment
#   ./table_5.sh --status                  Show completion status
#   ./table_5.sh --save                    Generate LaTeX table only

set -e

# Parse command line arguments
MODE="partial"
ACTION="run"
SPECIFIC_EXP=""

while [ $# -gt 0 ]; do
  case "$1" in
    --full) MODE="full"; shift ;;
    --partial) MODE="partial"; shift ;;
    --experiment) ACTION="run-one"; SPECIFIC_EXP="$2"; shift 2 ;;
    --status) ACTION="status"; shift ;;
    --save) ACTION="save"; shift ;;
    *)
      echo "Usage: $0 [OPTIONS]"
      echo ""
      echo "Options:"
      echo "  --partial           Use -partial control files (default)"
      echo "  --full              Use regular control files"
      echo "  --experiment CONFIG Run specific experiment (e.g., mc-02-0.0)"
      echo "  --status            Show what's complete"
      echo "  --save              Generate table without running experiments"
      exit 1
      ;;
  esac
done

# Set suffix based on mode
SUFFIX=""
if [ "$MODE" = "partial" ]; then
    SUFFIX="-partial"
fi

# Check for required binary
if [ "$ACTION" != "status" ] && [ "$ACTION" != "save" ]; then
  if [ ! -x ./mcnp ]; then
      echo "mcnp binary not found. Building..."
      make mcnp
      echo ""
  fi
fi

# Function to check if experiment is complete
is_complete() {
  local CONFIG=$1
  local CTL="control/${CONFIG}${SUFFIX}.ctl"
  local BASE="${CONFIG}${SUFFIX}"

  [ ! -f "$CTL" ] && return 1

  if ls results/${BASE}_*.txt &>/dev/null 2>&1; then
    NREPS=$(cat results/${BASE}_*.txt 2>/dev/null | grep -v '^#' | wc -l | tr -d ' ')
    EXPECTED=$(grep "^nmc" "$CTL" | awk '{print $2}')
    [ "$NREPS" -eq "$EXPECTED" ] && return 0
  fi

  return 1
}

# Function to run experiment (with skip if complete)
run_experiment() {
  local CONFIG=$1
  local THREADS=$2
  local CTL="control/${CONFIG}${SUFFIX}.ctl"

  if [ ! -f "$CTL" ]; then
    echo "Warning: Control file $CTL not found, skipping..."
    return 1
  fi

  if is_complete "$CONFIG"; then
    echo "Skipping $CONFIG (already complete)"
    return 0
  fi

  echo ""
  echo "Running $CONFIG..."
  ./run_parallel.sh "$CTL" "$THREADS"
}

# Main execution
case "$ACTION" in
  status)
    echo "Experiment Status ($MODE mode):"
    echo ""
    for CONFIG in mc-02-0.0 mc-02-1.0 mc-04-0.0 mc-04-1.0 mc-06-0.0 mc-06-1.0 mc-08-0.0 mc-08-1.0; do
      if is_complete "$CONFIG"; then
        echo "  [DONE] $CONFIG"
      else
        echo "  [TODO] $CONFIG"
      fi
    done
    ;;

  save)
    echo "Generating LaTeX table..."
    ./table_5_save.sh --$MODE
    ;;

  run-one)
    echo "Running experiment: $SPECIFIC_EXP (mode: $MODE)"
    echo ""

    case "$SPECIFIC_EXP" in
      mc-02-*) THREADS=2 ;;
      mc-04-*) THREADS=4 ;;
      mc-06-*) THREADS=4 ;;
      mc-08-*) THREADS=8 ;;
      *)
        echo "Error: Unknown experiment: $SPECIFIC_EXP"
        exit 1
        ;;
    esac

    run_experiment "$SPECIFIC_EXP" "$THREADS"
    echo ""
    echo "To generate table: $0 --save"
    ;;

  run)
    echo "Table 5: Quality Ladder Model Monte Carlo Results"
    echo "================================================="
    echo "Mode: $MODE"
    echo ""

    run_experiment "mc-02-0.0" 2
    run_experiment "mc-02-1.0" 2
    run_experiment "mc-04-0.0" 4
    run_experiment "mc-04-1.0" 4
    run_experiment "mc-06-0.0" 4
    run_experiment "mc-06-1.0" 4
    run_experiment "mc-08-0.0" 8
    run_experiment "mc-08-1.0" 8

    echo ""
    echo "All experiments complete!"
    echo "Generating LaTeX table..."
    echo ""

    ./table_5_save.sh --$MODE
    ;;
esac
