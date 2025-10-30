#!/bin/bash

# Main script to generate Table 5: Monte Carlo Results
# Runs all experiments and generates the LaTeX table

set -e

# Function to run a single Monte Carlo experiment
run_experiment() {
  local CONFIG=$1
  local THREADS=$2
  local CTL="control/${CONFIG}.ctl"

  if [ ! -f "$CTL" ]; then
    echo "Warning: Control file $CTL not found, skipping..." >&2
    return 1
  fi

  ./run_parallel.sh "$CTL" "$THREADS"
}

echo "Table 5: Monte Carlo Experiments"
echo "================================"
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

./make_mc_table.sh
