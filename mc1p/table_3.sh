#!/bin/bash

# Create directories for logs and results
LOGDIR=logs
RESULTSDIR=results
mkdir -p $LOGDIR
mkdir -p $RESULTSDIR

# Output file for Table 3
OUTFILE=$RESULTSDIR/table_3.tex

echo "Running Monte Carlo experiments for Table 3..."
echo ""

# Run all experiments and save logs
for nm in 200 800 3200; do
    for delta in 0.00 1.00 8.00; do
        echo "Running experiment: M=$nm, delta=$delta"
        ./mc1p control/mc-$delta-$nm.ctl | tee $LOGDIR/mc-$delta-$nm.log
    done
done

echo ""
echo "Generating Table 3..."

# Generate LaTeX table and save to file
{
    echo "\begin{table}[htbp]"
    echo " \centering"
    echo " \caption{Single Agent Renewal Model Monte Carlo Results}"
    echo " \label{tab:mc1p}"
    echo " \begin{tabular}{rllrrrrrr}"
    echo "   \hline"
    echo "   \hline"
    echo "       \$M$   & Sampling        &         &  \$\\lambda_1$ &  \$\\lambda_2$ &    \$\\gamma$  &     \$\\beta$  &       \$\\mu$  &  \$\\mu/\\beta$ \\\\"
    echo "   \hline"
    echo "   \$\\infty$  & DGP             & True    &       0.050  &       0.100  &       0.500  &      -2.000  &      -9.000  &       4.500  \\\\"
    echo "   \hline"

    for nm in 200 800 3200; do
        for delta in 0.00 1.00 8.00; do
            tail -5 $LOGDIR/mc-$delta-$nm.log | head -2
        done
        echo "     \\hline"
    done

    echo " \end{tabular}"
    echo " \medskip"
    echo " \begin{footnotesize}"
    echo "   The mean and standard deviation are reported for"
    echo "   100 replications under several sampling regimes."
    echo "   For each replication, \$M\$ markets were simulated over a fixed time"
    echo "   interval \$[0,T]\$ with \$T = 120\$."
    echo " \end{footnotesize}"
    echo "\end{table}"
} > $OUTFILE

# Display the table to console
echo "LaTeX Table"
echo "==========="
cat $OUTFILE

echo ""
echo "Monte Carlo log files saved to: $LOGDIR/"
echo "Table 3 saved to: $RESULTSDIR/table_3.tex"
