#!/bin/bash

# Generate Tables 1 and 2

set -e

echo "Table 2: Rust (1987) Model Estimates"
echo "===================================="
echo ""

# Create logs directory
mkdir -p logs
mkdir -p results

# Run three model variants
echo "Running ABBE model (lambda_L = lambda_H = 1.0)..."
#./rustct-abbe > logs/rustct-abbe.log 2>&1
echo "Done. Log saved to logs/rustct-abbe.log"
echo ""

echo "Running homogeneous model (lambda_L = lambda_H)..."
#./rustct-homogeneous > logs/rustct-homogeneous.log 2>&1
echo "Done. Log saved to logs/rustct-homogeneous.log"
echo ""

echo "Running heterogeneous model (lambda_L != lambda_H)..."
#./rustct-heterogeneous > logs/rustct-heterogeneous.log 2>&1
echo "Done. Log saved to logs/rustct-heterogeneous.log"
echo ""

# Extract results and generate table
echo "Extracting results and generating LaTeX table..."
echo "-------------------------------------------------"
echo ""

# Function to extract log likelihood
get_ll() {
  local LOG=$1
  grep "^Best starting value:" -A 1 "$LOG" | tail -1 | awk '{print $3}'
}

# Function to extract number of observations
get_nobs() {
  local LOG=$1
  grep "Loaded.*observations in total" "$LOG" | awk '{print $3}'
}

# Function to extract estimate value from table
# Format: "    Estimates  & value1  &   value2  & ..."
# Fields: Estimates & value1 & value2 & ...
# So value at column N is at field (2*N + 3)
get_estimate() {
  local LOG=$1
  local COL=$2
  grep "^    Estimates" "$LOG" | awk -v col="$COL" '{print $(2*col+3)}'
}

# Function to extract standard error from table
get_se() {
  local LOG=$1
  local COL=$2
  grep "^    S.E." "$LOG" | awk -v col="$COL" '{print $(2*col+3)}'
}

# Extract all values and format log-likelihoods to 2 decimal places
LL_ABBE=$(get_ll "logs/rustct-abbe.log" | awk '{printf "%.2f", $1}')
LL_HOMO=$(get_ll "logs/rustct-homogeneous.log" | awk '{printf "%.2f", $1}')
LL_HETE=$(get_ll "logs/rustct-heterogeneous.log" | awk '{printf "%.2f", $1}')

NOBS=$(get_nobs "logs/rustct-abbe.log")

# ABBE model: lambda=1.0 (fixed), q1, beta, c
ABBE_Q1=$(get_estimate "logs/rustct-abbe.log" 1)
ABBE_BETA=$(get_estimate "logs/rustct-abbe.log" 2)
ABBE_C=$(get_estimate "logs/rustct-abbe.log" 3)
ABBE_Q1_SE=$(get_se "logs/rustct-abbe.log" 1)
ABBE_BETA_SE=$(get_se "logs/rustct-abbe.log" 2)
ABBE_C_SE=$(get_se "logs/rustct-abbe.log" 3)

# Homogeneous model: lambda, q1, beta, c
HOMO_LAMBDA=$(get_estimate "logs/rustct-homogeneous.log" 0)
HOMO_Q1=$(get_estimate "logs/rustct-homogeneous.log" 1)
HOMO_BETA=$(get_estimate "logs/rustct-homogeneous.log" 2)
HOMO_C=$(get_estimate "logs/rustct-homogeneous.log" 3)
HOMO_LAMBDA_SE=$(get_se "logs/rustct-homogeneous.log" 0)
HOMO_Q1_SE=$(get_se "logs/rustct-homogeneous.log" 1)
HOMO_BETA_SE=$(get_se "logs/rustct-homogeneous.log" 2)
HOMO_C_SE=$(get_se "logs/rustct-homogeneous.log" 3)

# Heterogeneous model: lambda1, lambda2, q1, beta, c
HETE_LAMBDA1=$(get_estimate "logs/rustct-heterogeneous.log" 0)
HETE_LAMBDA2=$(get_estimate "logs/rustct-heterogeneous.log" 1)
HETE_Q1=$(get_estimate "logs/rustct-heterogeneous.log" 2)
HETE_BETA=$(get_estimate "logs/rustct-heterogeneous.log" 3)
HETE_C=$(get_estimate "logs/rustct-heterogeneous.log" 4)
HETE_LAMBDA1_SE=$(get_se "logs/rustct-heterogeneous.log" 0)
HETE_LAMBDA2_SE=$(get_se "logs/rustct-heterogeneous.log" 1)
HETE_Q1_SE=$(get_se "logs/rustct-heterogeneous.log" 2)
HETE_BETA_SE=$(get_se "logs/rustct-heterogeneous.log" 3)
HETE_C_SE=$(get_se "logs/rustct-heterogeneous.log" 4)

# Compute likelihood ratio statistics and p-values
echo "Computing likelihood ratio statistics..."

# Chi-squared p-value calculation using R
chi2_pvalue() {
  local LR=$1
  local DF=$2

  if command -v Rscript &> /dev/null; then
    Rscript --vanilla -e "cat(sprintf('%.5f', pchisq($LR, $DF, lower.tail=FALSE)))" 2>/dev/null
  else
    echo "--"
  fi
}

# Test 1a: H0: lambda = 1 (ABBE vs homogeneous), df=1
LR1a=$(awk -v ll_homo="$LL_HOMO" -v ll_abbe="$LL_ABBE" 'BEGIN {printf "%.2f", 2*(ll_homo - ll_abbe)}')
PVAL1a=$(chi2_pvalue "$LR1a" 1)

# Test 1b: H0: lambda_L = lambda_H = 1 (ABBE vs heterogeneous), df=2
LR1b=$(awk -v ll_hete="$LL_HETE" -v ll_abbe="$LL_ABBE" 'BEGIN {printf "%.2f", 2*(ll_hete - ll_abbe)}')
PVAL1b=$(chi2_pvalue "$LR1b" 2)

# Test 2: H0: lambda_L = lambda_H (homogeneous vs heterogeneous), df=1
LR2=$(awk -v ll_hete="$LL_HETE" -v ll_homo="$LL_HOMO" 'BEGIN {printf "%.2f", 2*(ll_hete - ll_homo)}')
PVAL2=$(chi2_pvalue "$LR2" 1)

echo "LR test 1a (ABBE vs homogeneous): LR = $LR1a, p = $PVAL1a"
echo "LR test 1b (ABBE vs heterogeneous): LR = $LR1b, p = $PVAL1b"
echo "LR test 2 (homogeneous vs heterogeneous): LR = $LR2, p = $PVAL2"
echo ""

# Generate LaTeX table
{
cat <<'EOF'
\begin{table}[tbh]
  \centering
  \begin{tabular}{lrcrcrc}
    \hline
    \hline
    & \multicolumn{2}{c}{Fixed $\lambda = 1$} & \multicolumn{2}{c}{Variable $\lambda$} & \multicolumn{2}{c}{Heterogeneous $\lambda$} \\
                                           & Est.      & S.E.    & Est.      & S.E.    & Est.      & S.E.    \\
    \hline
EOF

# Decision rates and mileage parameters
printf "    Decision rate   (\$\\lambda\$)            &     1.000 &      -- & %9s & (%s) &        -- &      -- \\\\\\\\\n" "$HOMO_LAMBDA" "$HOMO_LAMBDA_SE"
printf "    Decision rate 1 (\$\\lambda_{\\\\text{L}}\$) &        -- &      -- &        -- &      -- & %9s & (%s) \\\\\\\\\n" "$HETE_LAMBDA1" "$HETE_LAMBDA1_SE"
printf "    Decision rate 2 (\$\\lambda_{\\\\text{H}}\$) &        -- &      -- &        -- &      -- & %9s & (%s) \\\\\\\\\n" "$HETE_LAMBDA2" "$HETE_LAMBDA2_SE"
printf "    Mileage increase (\$\\gamma\$)            & %9s & (%s) & %9s & (%s) & %9s & (%s) \\\\\\\\\n" "$ABBE_Q1" "$ABBE_Q1_SE" "$HOMO_Q1" "$HOMO_Q1_SE" "$HETE_Q1" "$HETE_Q1_SE"
printf "    Mileage cost (\$\\\\beta\$)                 & %9s & (%s) & %9s & (%s) & %9s & (%s) \\\\\\\\\n" "$ABBE_BETA" "$ABBE_BETA_SE" "$HOMO_BETA" "$HOMO_BETA_SE" "$HETE_BETA" "$HETE_BETA_SE"
printf "    Replacement cost (\$\\mu\$)               & %9s & (%s) & %9s & (%s) & %9s & (%s) \\\\\\\\\n" "$ABBE_C" "$ABBE_C_SE" "$HOMO_C" "$HOMO_C_SE" "$HETE_C" "$HETE_C_SE"

cat <<EOF
    \hline
    Log likelihood                         & $LL_ABBE &         & $LL_HOMO &         & $LL_HETE &         \\\\
    Observations                           &     $NOBS &         &     $NOBS &         &     $NOBS &         \\\\
    \hline
    \\multicolumn{7}{l}{Test for \$H_0: \\lambda_{\\text{L}} = \\lambda_{\\text{H}}=1\$}                            \\\\
    \\quad LR                              &         -- &         &     $LR1a &         &     $LR1b &         \\\\
    \\quad \$p\$-value                       &         -- &         &   $PVAL1a &         &   $PVAL1b &         \\\\
    \\multicolumn{7}{l}{Test for \$H_0: \\lambda_{\\text{L}} = \\lambda_{\\text{H}}\$}                              \\\\
    \\quad LR                              &         -- &         &        -- &         &      $LR2 &         \\\\
    \\quad \$p\$-value                       &         -- &         &        -- &         &   $PVAL2 &         \\\\
    \hline
  \end{tabular}
  \caption{Model Estimates Based on Data from \cite{rust87optimal}}
  \label{tab:mc1p:rustct:estimates}
\end{table}
EOF
} | tee results/table_2.tex

echo ""
echo "Table 1 written to results/table_1.tex"
echo "Table 2 written to results/table_2.tex"
