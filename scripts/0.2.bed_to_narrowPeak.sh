#!/bin/bash

# Input parameters
peaks=$1
out_peaks=$2
flank_size=$3


# Process the peaks file with filtering and formatting
awk -v flank_size="$flank_size" 'BEGIN {OFS="\t"} 
    $2 > flank_size { 
        col4=$1 ":" $2 "-" $3; 
        col5=0; 
        col6="."; 
        col7=0; 
        col8=-1; 
        col9=-1; 
        col10=int((($3 - $2) / 2) + 0.5); 
        print $1, $2, $3, col4, col5, col6, col7, col8, col9, col10 
    }' "$peaks" > "$out_peaks"

echo "Processing complete. Output written to $out_peaks"
