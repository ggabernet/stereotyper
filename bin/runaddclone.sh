#!/bin/bash

REP=("HD07" "HD09" "HD10" "HD13" "F09" "F15" "F48" "F93" "F94" "F67" "HD-6" "HD-8")
SEQ = ("1" "1" "1" "1" "1" "1" "1" "1" "1" "2" "1" "1")
for i in "${!REP[@]}"
do
    echo "Processing ${REP[$i]}..."
    ~/Projects/Pipelines/stereotyper/bin/addclone_to_rep.R --repertoire "../../aggregated_healthy_controls/${REP[i]}__clone-pass.tsv" \\
    --clone "../add_clone_repertoire/${REP[i]}/${REP[i]}_naive_seqs.${SEQ[i]}_translated.tsv" --target_sequence QVQLVQSGAEVKKPGSSVKVSCKASGGTFNSYAISWVRQAPGQGLEWMGGIIPIFGTAKYAQKFQGRVTITADESTGTAYMQLNSLRSDDTAVFYCARVGGELGYCRSIDCFLFYGMDAWGQGTTVTVSS --abund 0.01 --outname "${REP[i]}_repertoire_with_naive_seqs.${SEQ[i]}_clone_sampled_"
done
