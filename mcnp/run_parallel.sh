#!/bin/bash

# Parallel Monte Carlo execution script
# Usage: ./run_parallel.sh control/mc-<nn>-<delta>.ctl [num_threads]
#
# num_threads: OpenMP threads per process (default: 2)
# num_processes: Computed as total_cores / num_threads
# total_cores is detected with the nprocs command, but can be
# overridden by setting the OMP_NUM_PROCS environment variable.

set -e

if [ $# -lt 1 ]; then
  echo "Usage: $0 <control_file> [num_threads]"
  echo "Example: $0 control/mc-08-0.0.ctl 8"
  echo ""
  echo "num_threads controls OMP_NUM_THREADS (inner parallelism)"
  echo "num_processes is computed as: available_cores / num_threads"
  exit 1
fi

CTL="$1"
NUM_THREADS=${2:-2}

# Detect available processors
if [ -n "$OMP_NUM_PROCS" ]; then
  TOTAL_CORES=$OMP_NUM_PROCS
elif command -v nproc &> /dev/null; then
  TOTAL_CORES=$(nproc)
elif command -v sysctl &> /dev/null; then
  TOTAL_CORES=$(sysctl -n hw.ncpu)
else
  TOTAL_CORES=2
  echo "Warning: Could not detect CPU count, assuming $TOTAL_CORES cores"
fi

# Calculate number of parallel processes
NPROC=$((TOTAL_CORES / NUM_THREADS))
if [ $NPROC -lt 1 ]; then
  NPROC=1
fi

if [ ! -f "$CTL" ]; then
  echo "Error: Control file not found: $CTL"
  exit 1
fi

# Read parameters from control file
NMC=$(grep "^nmc" "$CTL" | awk '{print $2}')
if [ -z "$NMC" ]; then
  echo "Error: Could not read nmc from control file"
  exit 1
fi

NN=$(grep "^nn" "$CTL" | awk '{print $2}')
DELTA=$(grep "^delta" "$CTL" | awk '{print $2}')

# Determine sampling type
if [ "$DELTA" = "0.00" ] || [ "$DELTA" = "0.0" ]; then
  SAMPLING="Continuous"
else
  SAMPLING="Discrete (Δ = $DELTA)"
fi

# Cap number of processes at number of replications
if [ $NPROC -gt $NMC ]; then
  NPROC=$NMC
fi

# Calculate replications per process with better load balancing
BASE_PER_PROC=$((NMC / NPROC))
REMAINDER=$((NMC % NPROC))

# Create directories
mkdir -p results logs

BASE=$(basename "$CTL" .ctl)

echo "Parallel Monte Carlo Execution"
echo "=============================="
echo ""
echo "Experiment:"
echo "  Number of firms: $NN"
echo "  Data sampling: $SAMPLING"
echo "  Replications: $NMC"
echo ""
echo "Parallelization:"
echo "  Total cores: $TOTAL_CORES"
echo "  Threads per process: $NUM_THREADS"
echo "  Number of processes: $NPROC"
echo "  Base replications per process: $BASE_PER_PROC"
if [ $REMAINDER -gt 0 ]; then
  echo "  First $REMAINDER processes get +1 extra replication"
fi
echo ""

# Launch all processes in background
PIDS=()
CURRENT_POS=1

for ((i=1; i<=NPROC; i++)); do
  START=$CURRENT_POS

  # First REMAINDER processes get one extra replication
  if [ $i -le $REMAINDER ]; then
    COUNT=$((BASE_PER_PROC + 1))
  else
    COUNT=$BASE_PER_PROC
  fi

  END=$((START + COUNT - 1))
  CURRENT_POS=$((END + 1))

  PROC_NUM=$(printf "%02d" $i)
  RESULTS_FILE="results/${BASE}_${PROC_NUM}.txt"
  LOG_FILE="logs/${BASE}_${PROC_NUM}.log"

  echo "Launching process $i: replications $START-$END"

  # Run mcnp with environment variables
  OMP_NUM_THREADS=$NUM_THREADS MC_START=$START MC_END=$END MC_RESULTS_FILE="$RESULTS_FILE" \
    ./mcnp "$CTL" > "$LOG_FILE" 2>&1 &

  PIDS+=($!)
done

WAIT_TIME=60
echo ""
echo "All processes launched. PIDs: ${PIDS[*]}"
echo "Started at: $(date)"
echo ""
echo "Monitoring progress (refreshes every $WAIT_TIME sec.)..."

# Wait for all jobs and track completion
START_TIME=$(date +%s)
COMPLETED=0
while [ $COMPLETED -lt $NPROC ]; do
  sleep $WAIT_TIME

  # Count completed processes
  COMPLETED=0
  for PID in "${PIDS[@]}"; do
    if ! kill -0 $PID 2>/dev/null; then
      COMPLETED=$((COMPLETED + 1))
    fi
  done

  # Count replications completed in each log file
  TOTAL_REPS_DONE=0
  for ((i=1; i<=NPROC; i++)); do
    PROC_NUM=$(printf "%02d" $i)
    LOG_FILE="logs/${BASE}_${PROC_NUM}.log"

    if [ -f "$LOG_FILE" ]; then
      REPS_DONE=$(grep "ll(theta)" "$LOG_FILE" 2>/dev/null | wc -l)
      REPS_DONE=$(echo $REPS_DONE)  # Strip whitespace
      TOTAL_REPS_DONE=$((TOTAL_REPS_DONE + REPS_DONE))
    fi
  done

  # Calculate progress percentage
  if [ $NMC -gt 0 ]; then
    PERCENT=$((TOTAL_REPS_DONE * 100 / NMC))
  else
    PERCENT=0
  fi

  # Calculate time estimate
  ELAPSED=$(($(date +%s) - START_TIME))
  TIMESTAMP=$(date)
  if [ $TOTAL_REPS_DONE -gt 0 ]; then
    AVG_TIME=$((ELAPSED / TOTAL_REPS_DONE))
    REMAINING=$((NMC - TOTAL_REPS_DONE))
    ETA=$((AVG_TIME * REMAINING))
    ETA_MIN=$((ETA / 60))
    printf "%s Completed %d/%d reps (~%dm remaining)..." "$TIMESTAMP" "$TOTAL_REPS_DONE" "$NMC" "$ETA_MIN"
  else
    printf "%s Completed %d/%d reps..." "$TIMESTAMP" "$TOTAL_REPS_DONE" "$NMC"
  fi
  echo ""
done

echo ""
echo "All processes completed at: $(date)"
echo ""

# Check for errors
echo "Checking for errors..."
FAILED=0
for ((i=1; i<=NPROC; i++)); do
  PROC_NUM=$(printf "%02d" $i)
  LOG_FILE="logs/${BASE}_${PROC_NUM}.log"

  if grep -qi "error\|segmentation\|aborted" "$LOG_FILE"; then
    echo "  Process $i FAILED - check $LOG_FILE"
    FAILED=$((FAILED + 1))
  else
    echo "  Process $i completed successfully"
  fi
done

if [ $FAILED -gt 0 ]; then
  echo ""
  echo "ERROR: $FAILED process(es) failed!"
  exit 1
fi

echo ""
echo "Parallel execution complete!"
echo ""
echo "Results written to:"
echo "  Partial results: results/${BASE}_*.txt"
echo "  Logs:            logs/${BASE}_*.log"
