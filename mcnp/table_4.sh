#!/bin/bash

# Main script to generate Table 4: Timing Results
# Runs timing experiments and generates the LaTeX table

set -e

echo "Table 4: Quality Ladder Model Timing Results"
echo "============================================"
echo ""

# Check for required binary and build if necessary
if [ ! -x ./mcnp ]; then
    echo "mcnp binary not found. Building..."
    make mcnp
    echo ""
fi

# Function to check if timing experiment is complete
is_complete() {
  local NN=$1
  local LOG="logs/time-${NN}.log"

  # Check if log exists and contains elapsed time (successful completion)
  if [ -f "$LOG" ] && grep -q "Elapsed time:" "$LOG" 2>/dev/null; then
    return 0
  fi

  return 1
}

# Function to run a single timing experiment
run_timing() {
  local NN=$1
  local CTL="control/time-${NN}.ctl"
  local LOG="logs/time-${NN}.log"

  if [ ! -f "$CTL" ]; then
    echo "Warning: Control file $CTL not found, skipping..." >&2
    return 1
  fi

  if is_complete "$NN"; then
    echo "Skipping time-${NN} (already complete)"
    return 0
  fi

  echo "Running time-${NN}..."

  # Create logs directory
  mkdir -p logs

  # Run experiment and save log
  ./mcnp "$CTL" > "$LOG" 2>&1

  echo "  Log saved to $LOG"
}


# Run all timing experiments
run_timing "02"
run_timing "04"
run_timing "06"
run_timing "08"
run_timing "10"
run_timing "12"
run_timing "14"
run_timing "16"
run_timing "18"
run_timing "20"
run_timing "22"
run_timing "24"
run_timing "26"
run_timing "28"
run_timing "30"

echo ""
echo "Timing experiments complete!"
echo "Generating LaTeX table..."
echo ""

./make_timing_table.sh
