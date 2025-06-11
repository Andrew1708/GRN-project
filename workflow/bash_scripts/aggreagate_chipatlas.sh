#!/bin/bash

# Usage:
# ./aggregate_beds.sh /path/to/bed_directory /path/to/output_file.bed

# Exit if any command fails
set -e

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/bed_directory /path/to/output_file.bed"
    exit 1
fi

# Assign arguments to variables
BED_DIR=$1
OUTPUT_FILE=$2

# Run the aggregation
cat "${BED_DIR}"/*.bed | python /home/andrem/GRN-project/resources/greta/tfb/aggregate.py > "${OUTPUT_FILE}"
