#!/usr/bin/env bash

# Usage:
# ./chipatlas_fetch_and_process.sh <url> <output> <input> <max_psize>

set -euo pipefail  # safer bash

# Parse command-line arguments
url="$1"
output="$2"
input="$3"
max_psize="$4"

# Check if URL is accessible
if wget --spider "$url" 2>/dev/null; then
    echo "Fetching data from: $url"
    
    # Download, process with python, merge with bedtools
    wget --no-verbose --tries=2 "$url" -O - | \
    python /home/andrem/GRN-project/resources/greta/tfb/chipatlas_tf.py "$output" "$input" "$max_psize" | \
    bedtools merge -c 4,5 -o distinct,distinct > "$output"

    echo "Finished processing. Output written to: $output"
else
    echo "URL not accessible: $url"
    touch "$output"  # create empty output to satisfy downstream
    echo "Created empty output: $output"
fi
