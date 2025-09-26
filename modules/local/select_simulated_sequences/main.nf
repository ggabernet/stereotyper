process SELECT_SIMULATED_SEQUENCES {
    tag "$meta.id"
    publishDir "$params.outdir/simulation/${meta.id}_f${fuziness}_a${abundance}_s${repertoire_sample}_w${witness}", mode: 'copy'
    container "docker.io/ggabernet/stereotyper:dev"

    input:
    tuple val(meta), path(repertoire), path(simulated_sequences), path(repertoire_embedding), path(simulated_embedding)

    output:
    tuple val(meta), path("*_repertoire_with_simulated_meta.tsv"), emit: rep_sim_sequences
    tuple val(meta), path("*_repertoire_with_simulated_embedding.tsv"), emit: rep_sim_embedding
    path("*.png"), emit: figures
    path "versions.yml", emit: versions


    script:
    """
    select_simulated_sequences.py --dir . \
    --repertoire_id ${meta.id} \
    --repertoire_embedding ${repertoire_embedding} \
    --simulated_embedding ${simulated_embedding} \
    --repertoire ${repertoire} \
    --simulated ${simulated_sequences} \
    --target_aa ${params.target_seq_aa} \
    --abundance ${meta.abundance} \
    --repertoire_sample ${meta.repertoire_sample} \
    --fuzziness_param ${meta.fuzziness} \
    --witness ${meta.witness} \
    --random_seed ${params.random_seed}

    # Save versions
    cat <<EOF > versions.yml
    python: \$(python --version 2>&1 | head -n 1)
    pandas: \$(python -c "import pandas; print(pandas.__version__)" | head -n 1)
    r: \$(R --version | head -n 1)
    rpy2: \$(python -c "import rpy2; print(rpy2.__version__)" | head -n 1)
    seaborn: \$(python -c "import seaborn; print(seaborn.__version__)" | head -n 1)
    matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)" | head -n 1)
    anndata: \$(python -c "import anndata; print(anndata.__version__)" | head -n 1)
    scanpy: \$(python -c "import scanpy; print(scanpy.__version__)" | head -n 1)
    scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)" | head -n 1)
    umap-learn: \$(python -c "import umap; print(umap.__version__)" | head -n 1)
    numpy: \$(python -c "import numpy; print(numpy.__version__)" | head -n 1)
    EOF
    """
}
