#!/bin/bash

# Generate LaTeX table from Monte Carlo results
#
# This script:
# 1. Discovers available results files in results/
# 2. Combines partial results for each configuration
# 3. Computes statistics (mean, sd)
# 4. Generates LaTeX table

set -e

# Combine partial results and compute statistics
process_config() {
  local BASE=$1  # e.g., "mc-02-0.0"

  # Check if any partial results exist
  if ! ls results/${BASE}_*.txt &>/dev/null; then
    return 1
  fi

  # Combine partial results
  COMBINED="results/${BASE}-combined.txt"
  cat results/${BASE}_*.txt | grep -v '^#' | sort -n > "$COMBINED"

  # Compute statistics
  STATS="results/${BASE}-stats.txt"
  awk '
  BEGIN {
    n = 0
  }
  {
    n++
    for (i=1; i<=6; i++) {
      sum[i] += $i
      sumsq[i] += $i * $i
    }
  }
  END {
    if (n == 0) exit 1

    params[1] = "lambda_L"
    params[2] = "lambda_H"
    params[3] = "gamma"
    params[4] = "kappa"
    params[5] = "eta"
    params[6] = "fc"

    for (i=1; i<=6; i++) {
      mean = sum[i] / n
      if (n > 1) {
        var = (sumsq[i] - sum[i]*sum[i]/n) / (n-1)
        stddev = sqrt(var > 0 ? var : 0)
      } else {
        stddev = 0
      }
      printf "%-15s %15.6f %15.6f\n", params[i], mean, stddev
    }
  }
  ' "$COMBINED" > "$STATS"

  return 0
}

# Function to extract K from log file
get_K() {
  local BASE=$1
  local LOG=$(ls logs/${BASE}_*.log 2>/dev/null | head -1)

  if [ -f "$LOG" ]; then
    K=$(grep "Total states:" "$LOG" | head -1 | awk '{print $NF}')
    if [ -n "$K" ]; then
      # Add commas for readability
      echo "$K" | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'
      return 0
    fi
  fi

  # Fallback: use hardcoded values based on nn
  case "$BASE" in
    mc-02-*) echo "56" ;;
    mc-04-*) echo "840" ;;
    mc-06-*) echo "5,544" ;;
    mc-08-*) echo "24,024" ;;
    *) echo "?" ;;
  esac
  return 0
}

# Function to get nn from control file
get_nn() {
  local BASE=$1
  local CTL="control/${BASE}.ctl"

  if [ -f "$CTL" ]; then
    grep "^nn" "$CTL" | awk '{print $2}'
  else
    echo "?"
  fi
}

echo "Processing Monte Carlo results..." >&2
echo "" >&2

# Discover and process all configurations
for NN in 02 04 06 08; do
  for DELTA in 0.0 1.0; do
    BASE="mc-${NN}-${DELTA}"

    if process_config "$BASE"; then
      echo "  Processed $BASE" >&2
    fi
  done
done

echo "" >&2
echo "Generating LaTeX table..." >&2

# Generate LaTeX table (output to both console and file)
{
cat <<'EOF'
\begin{table}[tbp]
 \centering
 \caption{Quality Ladder Model Monte Carlo Results}
 \label{tab:ladder:nfxp}
 \begin{tabular}{rrllrrrrrr}
   \hline
   \hline
   $\Nplayers$ & $K$    & Sampling       &      & $\lambda_{\text{L}}$ & $\lambda_{\text{H}}$ & $\gamma$ & $\kappa$ & $\eta$ & $\mu$ \\
   \hline
EOF

# DGP parameters (true values)
echo "               &         & DGP            & True & 1.000 & 1.200 & 0.400 & 0.800 & 4.000 & 0.900 \\\\"
echo "   \\hline"

# Process each configuration
for NN in 02 04 06 08; do
  NN_VAL=$(get_nn "mc-${NN}-0.0")
  K=$(get_K "mc-${NN}-0.0")

  # Process continuous (delta=0.0)
  STATS_CONT="results/mc-${NN}-0.0-stats.txt"
  if [ -f "$STATS_CONT" ]; then
    LAMBDA_L_MEAN=$(grep "^lambda_L" "$STATS_CONT" | awk '{printf "%.3f", $2}')
    LAMBDA_L_STD=$(grep "^lambda_L" "$STATS_CONT" | awk '{printf "%.3f", $3}')
    LAMBDA_H_MEAN=$(grep "^lambda_H" "$STATS_CONT" | awk '{printf "%.3f", $2}')
    LAMBDA_H_STD=$(grep "^lambda_H" "$STATS_CONT" | awk '{printf "%.3f", $3}')
    GAMMA_MEAN=$(grep "^gamma" "$STATS_CONT" | awk '{printf "%.3f", $2}')
    GAMMA_STD=$(grep "^gamma" "$STATS_CONT" | awk '{printf "%.3f", $3}')
    KAPPA_MEAN=$(grep "^kappa" "$STATS_CONT" | awk '{printf "%.3f", $2}')
    KAPPA_STD=$(grep "^kappa" "$STATS_CONT" | awk '{printf "%.3f", $3}')
    ETA_MEAN=$(grep "^eta" "$STATS_CONT" | awk '{printf "%.3f", $2}')
    ETA_STD=$(grep "^eta" "$STATS_CONT" | awk '{printf "%.3f", $3}')
    FC_MEAN=$(grep "^fc" "$STATS_CONT" | awk '{printf "%.3f", $2}')
    FC_STD=$(grep "^fc" "$STATS_CONT" | awk '{printf "%.3f", $3}')

    echo "    $NN_VAL          & $K      & Continuous     & Mean & $LAMBDA_L_MEAN & $LAMBDA_H_MEAN & $GAMMA_MEAN & $KAPPA_MEAN & $ETA_MEAN & $FC_MEAN \\\\"
    echo "               &         &                & S.D. & $LAMBDA_L_STD & $LAMBDA_H_STD & $GAMMA_STD & $KAPPA_STD & $ETA_STD & $FC_STD \\\\"
  fi

  # Process discrete (delta=1.0)
  STATS_DISC="results/mc-${NN}-1.0-stats.txt"
  if [ -f "$STATS_DISC" ]; then
    LAMBDA_L_MEAN=$(grep "^lambda_L" "$STATS_DISC" | awk '{printf "%.3f", $2}')
    LAMBDA_L_STD=$(grep "^lambda_L" "$STATS_DISC" | awk '{printf "%.3f", $3}')
    LAMBDA_H_MEAN=$(grep "^lambda_H" "$STATS_DISC" | awk '{printf "%.3f", $2}')
    LAMBDA_H_STD=$(grep "^lambda_H" "$STATS_DISC" | awk '{printf "%.3f", $3}')
    GAMMA_MEAN=$(grep "^gamma" "$STATS_DISC" | awk '{printf "%.3f", $2}')
    GAMMA_STD=$(grep "^gamma" "$STATS_DISC" | awk '{printf "%.3f", $3}')
    KAPPA_MEAN=$(grep "^kappa" "$STATS_DISC" | awk '{printf "%.3f", $2}')
    KAPPA_STD=$(grep "^kappa" "$STATS_DISC" | awk '{printf "%.3f", $3}')
    ETA_MEAN=$(grep "^eta" "$STATS_DISC" | awk '{printf "%.3f", $2}')
    ETA_STD=$(grep "^eta" "$STATS_DISC" | awk '{printf "%.3f", $3}')
    FC_MEAN=$(grep "^fc" "$STATS_DISC" | awk '{printf "%.3f", $2}')
    FC_STD=$(grep "^fc" "$STATS_DISC" | awk '{printf "%.3f", $3}')

    echo "               &         & \$\\Delta = 1.0\$ & Mean & $LAMBDA_L_MEAN & $LAMBDA_H_MEAN & $GAMMA_MEAN & $KAPPA_MEAN & $ETA_MEAN & $FC_MEAN \\\\"
    echo "               &         &                & S.D. & $LAMBDA_L_STD & $LAMBDA_H_STD & $GAMMA_STD & $KAPPA_STD & $ETA_STD & $FC_STD \\\\"
  fi
done

cat <<'EOF'
   \hline
 \end{tabular}
\end{table}
EOF
} | tee results/table_5.tex

echo "" >&2
echo "Table written to results/table_5.tex" >&2
