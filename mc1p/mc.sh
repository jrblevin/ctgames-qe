#!/bin/bash

LOGDIR=logs
mkdir -p $LOGDIR

for nm in 200 800 3200; do
    for delta in 0.00 1.00 8.00; do
        ./mc1p control/mc-$delta-$nm.ctl | tee $LOGDIR/mc-$delta-$nm.log
    done
done

echo "\begin{table}[htbp]"
echo " \centering"
echo " \caption{Single Agent Renewal Model Monte Carlo Results}"
echo " \label{tab:mc1p}"
echo " \begin{tabular}{rllrrrrrr}"
echo "   \hline"
echo "   \hline"
echo "   \$M$   & Sampling        &      & $\lambda_1$ & $\lambda_2$ & $\gamma$  & $\beta$ & $\mu$    & $\mu/\beta$ \\\\"
echo "   \hline"
echo "   $\infty$ & DGP          & True & 0.050 & 0.100 & 0.500 & -2.000 & -9.000 & 4.500 \\\\"
echo "   \hline"

for nm in 200 800 3200; do
    for delta in 0.00 1.00 8.00; do
        tail -5 $LOGDIR/mc-$delta-$nm.log | head -2
    done
    echo "     \\hline"
done

echo " \end{tabular}                                                              \\\\"
echo " \medskip"
echo " \begin{footnotesize}"
echo "   The mean and standard deviation are reported for"
echo "   100 replications under several sampling regimes."
echo "   For each replication, \$M$ markets were simulated over a fixed time"
echo "   interval \$[0,T]$ with \$T = 120$."
echo " \end{footnotesize}"
echo "\end{table}"
