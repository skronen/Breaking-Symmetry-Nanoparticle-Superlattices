#!/bin/bash

input_file="pairs.txt"
output_file="filtered_pairs.txt"
cutoffval=39

awk -v cutoffval="$cutoffval" '
BEGIN { skip_block = 0; }
/^ITEM: TIMESTEP/ {
    if (!skip_block && timestep_seen) {
        print timestep_line
        print timestep
        for (i in data_block) print data_block[i]
    }
    delete data_block
    timestep_seen = 1
    skip_block = 0
    getline timestep
    timestep_line = $0
    next
}
/^ITEM: NUMBER OF ENTRIES/ {
    getline num_entries
    if (num_entries == 0) {
        skip_block = 1
    }
    next
}
/^ITEM: BOX BOUNDS/ {
    for (i = 0; i < 3; i++) getline
    next
}
/^ITEM:/ { next }
{
    if (!skip_block) {
        if (NF == 5) {
            first = $1
            second_last = $(NF-1)
            last = $NF
            # keep only if first <= cutoffval AND last two values are (5,10) or (10,5)
            if (first <= cutoffval && 
                ((second_last == 5 && last == 10) || (second_last == 10 && last == 5))) {
                data_block[length(data_block)] = $0
            }
        } else {
            # keep any non-5-field lines untouched
            data_block[length(data_block)] = $0
        }
    }
}
END {
    if (!skip_block && timestep_seen) {
        print timestep_line
        print timestep
        for (i in data_block) print data_block[i]
    }
}
' "$input_file" > "$output_file"

