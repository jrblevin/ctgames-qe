#!/bin/sh

# Exit on any error
set -e

#
# Build binaries
#

# Function to check if a binary exists
check_binary() {
    if [ ! -x "$1" ]; then
        echo "Error: Binary $1 not found or not executable!"
        echo "Build may have failed. Please check the build output above."
        exit 1
    else
        echo "Binary $1 found and is executable"
    fi
}

# Build rustct variants and mc1p
echo "Building mc1p and rustct binaries in mc1p/ directory"
echo "===================================================="
echo ""
cd mc1p
make mc1p rustct-heterogeneous rustct-homogeneous rustct-abbe
echo ""
check_binary "rustct-heterogeneous"
check_binary "rustct-homogeneous"
check_binary "rustct-abbe"
check_binary "mc1p"
cd ..

# Build mcnp
echo ""
echo "Building mcnp binary in mcnp/ directory"
echo "======================================="
echo ""
cd mcnp
make mcnp
echo ""
check_binary "mcnp"
cd ..

#
# Single agent model results
#

# Table 1 and Table 2
echo "Generating Table 1 and Table 2"
echo "=============================="
cd mc1p
if [ -x ./rustct-abbe ] && [ -x ./rustct-homogeneous ] && [ -x ./rustct-heterogeneous ]; then
    echo "Running table_1_2.sh for Table 1 (sample characteristics) and Table 2 (estimates)..."
    ./table_1_2.sh
else
    echo "Error: rustct binaries not found in mc1p/"
    echo "Please ensure rustct-abbe, rustct-homogeneous, and rustct-heterogeneous are built."
    exit 1
fi
cd ..

# Table 3
echo ""
echo "Generating Table 3"
echo "=================="
cd mc1p
if [ -x ./mc1p ]; then
    echo "Running Table 3 experiments..."
    ./table_3.sh
else
    echo "Error: mc1p binary not found in mc1p/"
    exit 1
fi
cd ..

#
# Quality ladder model results
#

# Table 4
echo ""
echo "Generating Table 4"
echo "=================="
cd mcnp
if [ -x ./mcnp ]; then
    echo "Running Table 4 timing experiments..."
    ./table_4.sh
else
    echo "Error: mcnp binary not found in mcnp/"
    exit 1
fi
cd ..

# Table 5
echo ""
echo "Generating Table 5"
echo "=================="
cd mcnp
if [ -x ./mcnp ]; then
    echo "Running Table 5 Monte Carlo experiments..."
    ./table_5.sh
else
    echo "Error: mcnp binary not found in mcnp/"
    exit 1
fi
cd ..

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
