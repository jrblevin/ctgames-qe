#!/bin/bash

# Generate LaTeX table from timing results.
# Extracts timing data from log files.

set -e

# Function to extract value from log file
extract_value() {
  local LOG=$1
  local PATTERN=$2
  local FIELD=$3

  if [ ! -f "$LOG" ]; then
    echo "?"
    return 1
  fi

  grep "$PATTERN" "$LOG" | head -1 | awk "{print \$$FIELD}"
}

# Function to extract elapsed time
extract_time() {
  local LOG=$1

  if [ ! -f "$LOG" ]; then
    echo "?"
    return 1
  fi

  # Extract "Elapsed time: X.XXX sec."
  grep "Obtaining value function" -A 1 "$LOG" | grep "Elapsed time:" | awk '{print $3}'
}

# Function to extract parameter from control file
get_control_param() {
  local NN=$1
  local PARAM=$2
  local CTL="control/time-${NN}.ctl"

  if [ -f "$CTL" ]; then
    grep "^${PARAM}" "$CTL" | awk '{print $2}'
  else
    echo "?"
  fi
}

echo "Processing timing results..." >&2
echo "" >&2

# Generate LaTeX table
{
cat <<'EOF'
\begin{table}[tbp]
 \centering
 \caption{Quality Ladder Model: Timing Results}
 \label{tab:ladder:timing}
 \begin{tabular}{rrrrr}
   \hline
   \hline
   $N$ & $\bar{\omega}$ & $K$ & $\bar{M}$ & Obtain $V$ (sec.) \\
   \hline
EOF

# Process each timing experiment
for NN in 02 04 06 08 10 12 14 16 18 20 22 24 26 28 30; do
  LOG="logs/time-${NN}.log"

  if [ ! -f "$LOG" ]; then
    echo "  Warning: Log file $LOG not found" >&2
    continue
  fi

  # Extract from control file
  N=$(get_control_param "$NN" "nn")
  NW=$(get_control_param "$NN" "nw")
  MKTSIZE=$(get_control_param "$NN" "mktsize")

  # Number of quality levels
  if [ "$NW" != "?" ] && [ -n "$NW" ]; then
    OMEGA=$(printf "%.0f" "$NW")
  else
    OMEGA="?"
  fi

  # Market size
  if [ "$MKTSIZE" != "?" ] && [ -n "$MKTSIZE" ]; then
    M=$(awk "BEGIN {printf \"%.2f\", $MKTSIZE}")
  else
    M="?"
  fi

  # Extract from log
  K=$(extract_value "$LOG" "Total states:" 3)
  TIME=$(extract_time "$LOG")

  # Format K with commas using awk
  if [ "$K" != "?" ] && [ -n "$K" ]; then
    K=$(echo "$K" | awk '{printf "%'"'"'d\n", $1}' 2>/dev/null || echo "$K")
  fi

  # Format time to 2 decimal places
  if [ "$TIME" != "?" ] && [ -n "$TIME" ]; then
    TIME=$(awk "BEGIN {printf \"%.2f\", $TIME}")
  else
    TIME="?"
  fi

  echo "    $N & $OMEGA & $K & $M & $TIME \\\\"
done

cat <<'EOF'
   \hline
 \end{tabular}
\end{table}
EOF
} | tee results/table_4.tex

echo "" >&2
echo "Table written to results/table_4.tex" >&2
