#!/bin/bash
# Batch script to run ABACUS calculations on multiple conformations using Docker
#
# Optimized defaults for Mac Mini M4:
#   - M4 base: 10 cores (4 performance + 6 efficiency)
#   - M4 Pro:  14 cores (10 performance + 4 efficiency)
#   - M4 Max:  16 cores (12 performance + 4 efficiency)
#
# Recommended settings:
#   MPI_PROCS=8, OMP_NUM_THREADS=1 (best for most ABACUS calculations)
#   Total cores used = MPI × OMP = 8 (leaves headroom for system)

set -e

# If not already running under caffeinate, re-execute with caffeinate
if [ -z "$CAFFEINATE_ACTIVE" ]; then
    export CAFFEINATE_ACTIVE=1
    exec caffeinate -i -d "$0" "$@"
fi

# Get the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_SCRIPT="$SCRIPT_DIR/4_run_abacus_docker.sh"

# Parse command line arguments
BASE_DIR="${1:-data/abacus_inputs}"
MAX_CONFORMATIONS="${2:-}"
MPI_ARG="${3:-}"
OMP_ARG="${4:-}"

# Set MPI processes (argument > env var > default)
# Default: 8 for Mac Mini M4 (uses 8 of 10 cores)
if [ -n "$MPI_ARG" ]; then
    export MPI_PROCS="$MPI_ARG"
else
    export MPI_PROCS=${MPI_PROCS:-4}
fi

# Set OpenMP threads (argument > env var > default)
# Default: 1 (MPI scales better than OpenMP for ABACUS)
if [ -n "$OMP_ARG" ]; then
    export OMP_NUM_THREADS="$OMP_ARG"
    OMP_THREADS="$OMP_ARG"
else
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
    OMP_THREADS=${OMP_NUM_THREADS:-2}
fi

# Limit number of conformations if specified
if [ -n "$MAX_CONFORMATIONS" ]; then
    if ! [[ "$MAX_CONFORMATIONS" =~ ^[0-9]+$ ]]; then
        echo "Error: MAX_CONFORMATIONS must be a number"
        echo "Usage: $0 [base_directory] [max_conformations] [mpi_procs] [omp_threads]"
        echo "Example: $0 data/abacus_inputs 50 8 4"
        exit 1
    fi
fi

# Check if base directory exists
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Directory '$BASE_DIR' does not exist"
    echo "Usage: $0 [base_directory] [max_conformations] [mpi_procs] [omp_threads]"
    echo "Example: $0 data/abacus_inputs 50 8 4"
    exit 1
fi

# Check if Docker script exists
if [ ! -f "$DOCKER_SCRIPT" ]; then
    echo "Error: Docker script not found: $DOCKER_SCRIPT"
    exit 1
fi

echo "=========================================="
echo "Batch ABACUS Calculations with Docker"
echo "=========================================="
echo "Base directory: $BASE_DIR"
echo "MPI processes: $MPI_PROCS"
echo "OpenMP threads: $OMP_THREADS"
if [ -n "$MAX_CONFORMATIONS" ]; then
    echo "Max conformations: $MAX_CONFORMATIONS"
fi
echo "=========================================="
echo ""

# Find all conformation directories (directories with INPUT file, excluding OUT.ABACUS)
CONF_DIRS=($(find "$BASE_DIR" -type f -name "INPUT" -not -path "*/OUT.ABACUS/*" -exec dirname {} \; | sort))

if [ ${#CONF_DIRS[@]} -eq 0 ]; then
    echo "Error: No INPUT files found in $BASE_DIR"
    echo "Expected structure:"
    echo "  $BASE_DIR/"
    echo "    fragment_XXX_conf_000/"
    echo "      INPUT"
    echo "      STRU"
    echo "      KPT"
    echo "      ..."
    exit 1
fi

# Limit if specified
if [ -n "$MAX_CONFORMATIONS" ] && [ ${#CONF_DIRS[@]} -gt $MAX_CONFORMATIONS ]; then
    CONF_DIRS=("${CONF_DIRS[@]:0:$MAX_CONFORMATIONS}")
fi

TOTAL=${#CONF_DIRS[@]}
SUCCESS=0
FAILED=0
SKIPPED=0

echo "Found $TOTAL conformations to process"
echo ""
echo "Note: Running with caffeinate to prevent system sleep"
echo "Note: Already completed simulations will be skipped"
echo ""

# Function to check if simulation is already complete
is_simulation_complete() {
    local conf_dir="$1"
    local log_file="$conf_dir/OUT.ABACUS/running_md.log"
    local md_dump="$conf_dir/OUT.ABACUS/MD_dump"
    
    # Check if MD_dump exists AND log file shows completion
    if [ -f "$md_dump" ] && [ -f "$log_file" ]; then
        # Check for "Finish Time" or "Total  Time" in the last 50 lines
        if tail -50 "$log_file" 2>/dev/null | grep -q "Finish Time\|Total  Time"; then
            return 0  # Complete
        fi
    fi
    return 1  # Not complete
}

# Run ABACUS on each conformation
for i in "${!CONF_DIRS[@]}"; do
    conf_dir="${CONF_DIRS[$i]}"
    conf_name=$(basename "$conf_dir")
    current=$((i + 1))
    
    echo "[$current/$TOTAL] Processing $conf_name..."
    
    # Check if INPUT file exists
    if [ ! -f "$conf_dir/INPUT" ]; then
        echo "  ✗ Skipping: No INPUT file found"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Check if simulation is already complete
    if is_simulation_complete "$conf_dir"; then
        echo "  ⏭ Skipping: Already completed (MD_dump + Finish Time found)"
        SKIPPED=$((SKIPPED + 1))
        SUCCESS=$((SUCCESS + 1))
        continue
    fi
    
    # Run ABACUS using the Docker script
    "$DOCKER_SCRIPT" "$conf_dir" > /tmp/abacus_batch_${current}.log 2>&1
    DOCKER_EXIT_CODE=$?
    
    # Check if output directory was created and has required files
    if [ -d "$conf_dir/OUT.ABACUS" ] && [ -f "$conf_dir/OUT.ABACUS/MD_dump" ]; then
        SUCCESS=$((SUCCESS + 1))
        echo "  ✓ Completed: $conf_name (MD_dump found)"
    elif [ -d "$conf_dir/OUT.ABACUS" ] && ls "$conf_dir/OUT.ABACUS"/running_*.log > /dev/null 2>&1; then
        # Has log but no MD_dump - might be incomplete or failed
        FAILED=$((FAILED + 1))
        echo "  ✗ Failed: $conf_name (no MD_dump file - check log)"
        echo "    Docker exit code: $DOCKER_EXIT_CODE"
        if [ -f /tmp/abacus_batch_${current}.log ]; then
            echo "    Last 5 lines of log:"
            tail -5 /tmp/abacus_batch_${current}.log | sed 's/^/      /'
        fi
    else
        FAILED=$((FAILED + 1))
        echo "  ✗ Failed: $conf_name (no output directory)"
        echo "    Docker exit code: $DOCKER_EXIT_CODE"
        if [ -f /tmp/abacus_batch_${current}.log ]; then
            echo "    Last 10 lines of log:"
            tail -10 /tmp/abacus_batch_${current}.log | sed 's/^/      /'
        fi
    fi
    
    echo ""
done

# Summary
echo "=========================================="
echo "Summary:"
echo "  Total: $TOTAL"
echo "  Success: $SUCCESS"
echo "  Skipped (already done): $SKIPPED"
echo "  Failed: $FAILED"
echo "  New completions: $((SUCCESS - SKIPPED))"
echo "=========================================="

if [ $FAILED -eq 0 ]; then
    echo ""
    echo "✓ All calculations completed successfully!"
    exit 0
else
    echo ""
    echo "⚠ Some calculations failed. Check logs above for details."
    exit 1
fi

