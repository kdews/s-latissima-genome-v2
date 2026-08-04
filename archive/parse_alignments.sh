#!/bin/bash

cactus_dir="cactus_pangenome_01-28-26_6039802"
mg_log="$cactus_dir/chrom-subproblems/minigraph.split.log"

assign_ptn="^Assigned contig to "
ambig_ptn="^Query contig is ambiguous"


grep -P "$ambig_ptn" "$mg_log" | head
grep -P "$assign_ptn" "$mg_log" \
  | sed "s/$assign_ptn//g" \
  | awk -F ":" '{print $1, $2}' \
  | sed "s/id=\w+\|//g"


grep -P "$ambig_ptn" "$mg_log" | sed "s/$ambig_ptn//g" | awk -F ":" '{print $1, $2}'

grep -h "seq00000000" "$mg_log" | awk -F ":" '{print $2}' | sed ''
