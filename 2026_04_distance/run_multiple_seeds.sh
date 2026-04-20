#!/bin/bash

# Script to run result3 analysis with multiple seeds (1-10)
# Each run will have its own output directory within a user-specified run directory
# This allows running multiple parameter configurations simultaneously

# Base directory
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG_FILE="${BASE_DIR}/config/analysis_config.R"
MAIN_SCRIPT="${BASE_DIR}/main.R"

# Check if files exist
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found: $CONFIG_FILE"
    exit 1
fi

if [ ! -f "$MAIN_SCRIPT" ]; then
    echo "Error: Main script not found: $MAIN_SCRIPT"
    exit 1
fi

# Parse command line arguments
RUN_DIR=""
if [ $# -eq 1 ]; then
    RUN_DIR="$1"
    echo "Using specified run directory: $RUN_DIR"
else
    # Create a timestamp-based run directory
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    RUN_DIR="${BASE_DIR}/run_${TIMESTAMP}"
    echo "No run directory specified. Creating new run directory: $RUN_DIR"
fi

# Create run directory if it doesn't exist
mkdir -p "$RUN_DIR"
echo "Run directory: $RUN_DIR"

# Create a unique backup of the current config file in the run directory
CONFIG_BACKUP="${RUN_DIR}/analysis_config.R.backup"
echo "Creating config backup in run directory: $CONFIG_BACKUP"
cp "$CONFIG_FILE" "$CONFIG_BACKUP"

echo "Starting batch run of result3 analysis with seeds 1-10"
echo "======================================================"
echo "Run directory: $RUN_DIR"
echo "Config backup: $CONFIG_BACKUP"
echo ""

# Loop through seeds 1 to 10
for SEED in {1..10}; do
    echo ""
    echo "Processing seed: $SEED"
    echo "----------------------"
    
    # Create a temporary copy of the config file for this seed
    TEMP_CONFIG="${RUN_DIR}/analysis_config_seed${SEED}.R"
    echo "Creating temporary config for seed $SEED: $TEMP_CONFIG"
    cp "$CONFIG_FILE" "$TEMP_CONFIG"
    
    # Update the seed in the temporary config file
    echo "Updating seed to $SEED in temporary config..."
    sed -i "s/seed = [0-9]*,/seed = ${SEED},/" "$TEMP_CONFIG"
    
    # Also update the set.seed line at the end of the file
    sed -i 's/set.seed(cfg\$seed)/set.seed('"$SEED"')/' "$TEMP_CONFIG"
    
    # Create output directory for this seed within the run directory
    OUTPUT_DIR="${RUN_DIR}/result3_seed${SEED}"
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
    
    # Save the temporary config to the output directory for reference
    cp "$TEMP_CONFIG" "${OUTPUT_DIR}/analysis_config_used.R"
    
    # Run the analysis using the temporary config file
    echo "Running: Rscript main.R result3 (using temporary config)"
    cd "$BASE_DIR"
    
    # First, replace the original config with our temporary one
    mv "$CONFIG_FILE" "${CONFIG_FILE}.original"
    cp "$TEMP_CONFIG" "$CONFIG_FILE"
    
    # Run the analysis and capture output
    Rscript main.R result3 2>&1 | tee "${OUTPUT_DIR}/run_log.txt"
    
    # Restore the original config file immediately
    mv "${CONFIG_FILE}.original" "$CONFIG_FILE"
    
    # Check if the run was successful
    if [ $? -eq 0 ]; then
        echo "Analysis completed successfully for seed $SEED"
        
        # Copy results to the seed-specific directory
        echo "Copying results to ${OUTPUT_DIR}/"
        
        # Copy output directories if they exist
        if [ -d "outputs" ]; then
            cp -r outputs/* "$OUTPUT_DIR/" 2>/dev/null || true
            echo "  Copied outputs directory contents"
        fi
        
        if [ -d "logs" ]; then
            cp -r logs/* "$OUTPUT_DIR/" 2>/dev/null || true
            echo "  Copied logs directory contents"
        fi
        
        # Also copy any RDS files
        if [ -d "outputs/rds" ]; then
            find outputs/rds -name "*.rds" -exec cp {} "$OUTPUT_DIR/" \; 2>/dev/null || true
            echo "  Copied RDS files"
        fi
        
        # Create a summary file
        echo "Seed: $SEED" > "${OUTPUT_DIR}/seed_info.txt"
        echo "Run directory: $RUN_DIR" >> "${OUTPUT_DIR}/seed_info.txt"
        echo "Config used: analysis_config_used.R" >> "${OUTPUT_DIR}/seed_info.txt"
        echo "Completed: $(date)" >> "${OUTPUT_DIR}/seed_info.txt"
    else
        echo "Warning: Analysis failed for seed $SEED"
        echo "Check ${OUTPUT_DIR}/run_log.txt for details"
        
        # Create error info file
        echo "Seed: $SEED - FAILED" > "${OUTPUT_DIR}/error_info.txt"
        echo "Error time: $(date)" >> "${OUTPUT_DIR}/error_info.txt"
    fi
    
    # Clean up temporary config file
    rm -f "$TEMP_CONFIG"
    
    echo "Completed processing seed $SEED"
    echo "==============================="
done

# Create a summary of the entire run
SUMMARY_FILE="${RUN_DIR}/run_summary.txt"
echo "Run Summary" > "$SUMMARY_FILE"
echo "===========" >> "$SUMMARY_FILE"
echo "Run directory: $RUN_DIR" >> "$SUMMARY_FILE"
echo "Original config backup: $(basename $CONFIG_BACKUP)" >> "$SUMMARY_FILE"
echo "Start time: $(date)" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "Seed directories:" >> "$SUMMARY_FILE"
for SEED in {1..10}; do
    SEED_DIR="${RUN_DIR}/result3_seed${SEED}"
    if [ -d "$SEED_DIR" ]; then
        if [ -f "${SEED_DIR}/seed_info.txt" ]; then
            STATUS="SUCCESS"
        else
            STATUS="FAILED or INCOMPLETE"
        fi
        echo "  result3_seed${SEED}: $STATUS" >> "$SUMMARY_FILE"
    else
        echo "  result3_seed${SEED}: NOT FOUND" >> "$SUMMARY_FILE"
    fi
done

echo ""
echo "Batch run completed!"
echo "===================="
echo "Run directory: $RUN_DIR"
echo "Contains:"
echo "  - analysis_config.R.backup: Backup of original config"
echo "  - result3_seed1 through result3_seed10: Results for each seed"
echo "  - run_summary.txt: Summary of the entire run"
echo ""
echo "To compare results across seeds, you can use:"
echo "  ls -la ${RUN_DIR}/result3_seed*/result3_summary.csv"
echo "  or write an R script to aggregate the results"
echo ""
echo "To run with different parameters:"
echo "  1. Modify ${CONFIG_FILE}"
echo "  2. Run: ./run_multiple_seeds_v2.sh /path/to/new/run/directory"
echo "  3. Or run without arguments to auto-create a timestamped directory"
