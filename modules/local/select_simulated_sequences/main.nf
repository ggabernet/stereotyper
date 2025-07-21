process SELECT_SIMULATED_SEQUENCES {
    tag "$meta.id"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'
    container "community.wave.seqera.io/library/anndata_matplotlib_numpy_pandas_pruned:7f590d3f128c4839"

    input:
    tuple val(meta), path(repertoire), path(repertoire_embedding)
    tuple val(meta), path(simulated_sequences), path(simulated_embedding)

    output:
    tuple val(meta), path("${meta.id}_selected_simulated_sequences.tsv"), emit: selected_simulated_sequences
    path "versions.yml", emit: versions


    script:
    """
    select_children.py --dir . \
    --repertoire_embedding ${repertoire_embedding} \
    --simulated_embedding ${simulated_embedding} \
    --repertoire ${repertoire} \
    --simulated ${simulated_sequences} \
    --embedding_model ${params.embedding_model} \
    --target_aa ${params.target_sequence_aa} \
    --abundance ${params.abundance_fraction} \
    --repertoire_sample ${params.subsample_size} \
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
