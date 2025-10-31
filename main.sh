#!/bin/sh

# For users with access to HPC resources, this script automates the
# complete replication of all tables.  This script runs each table in
# order (Tables 1-5), with parallelization within each table as appropriate.
#
# Given the computational requirements, particularly for Table 5,
# this script should only be used on HPC systems with substantial resources
# (32+ GB RAM, 40+ cores).
#
# For most users, we recommend running the individual table scripts
# separately, which provides more flexibility and control over the
# replication process.

# Exit on any error
set -e

#
# Single agent model results
#

cd mc1p

# Table 1 and Table 2
./table_1_2.sh

# Table 3
./table_3.sh

cd ..

#
# Quality ladder model results
#

cd mcnp

# Table 4
./table_4.sh

# Table 5
./table_5.sh

cd..

echo ""
echo "All tables generated successfully!"
echo ""
echo "Results are saved in:"
echo "  - mc1p/results/table_1.tex"
echo "  - mc1p/results/table_2.tex"
echo "  - mc1p/results/table_3.tex"
echo "  - mcnp/results/table_4.tex"
echo "  - mcnp/results/table_5.tex"
echo ""
