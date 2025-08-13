#!/usr/bin/env python
# Select children based on distance to target and to other centroids in the repertoire
# Usage: python select_children.py --dir <input_directory> --repertoire_embedding <repertoire_embedding.tsv> --simulated_embedding <simulated_embedding.tsv> --repertoire <repertoire.tsv> --simulated <simulated_sequences.tsv> --embedding_model <embedding_model> --target_aa <target_seq_aa> --abundance <abundance_fraction> --repertoire_sample <number_of_repertoire_samples> --random_seed <random_seed> --germline_vcall <germline_V_call> --germline_jcall <germline_J_call>

import pandas as pd
# uses latest unreleased amulety version
from rpy2.robjects.packages import importr
alakazam = importr('alakazam')


#from umap import UMAP
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import anndata as ad
import os
from pathlib import Path
import sklearn
import umap
import numpy as np
import logging
import argparse
from scipy.spatial.distance import cdist


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Select simulated sequences based on centroids and target embedding.")
parser.add_argument("--dir", type=str, default=".", help="results directory")
parser.add_argument("--repertoire_embedding", type=str, help="Path to repertoire embedded file")
parser.add_argument("--simulated_embedding", type=str, help="Path to simulated embedded file")
parser.add_argument("--repertoire", type=str, help="Path to repertoire AIRR file.")
parser.add_argument("--simulated", type=str, help="Path to simulated sequence and metadata file.")
parser.add_argument("--target_aa", type=str, default="GFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSENTLYLQMNSLRAEDTAVYYCARGEIQPYYYYGMDVWGQGTTVTVSS", help="Target amino acid sequence")
#parser.add_argument("--embedding_model", type=str, default="antiberty", help="Embedding model to use")
parser.add_argument("--abundance", type=float, default=0.01, help="Abundance fraction")
parser.add_argument("--repertoire_sample", type=int, default=100000, help="Number of repertoire samples")
parser.add_argument("--fuzziness_param", type=float, default=2.0, help="Parameter for distance weighting (should be in range (1, inf))")
parser.add_argument("--random_seed", type=int, default=42, help="Random seed")
args = parser.parse_args()


dir = args.dir
figures_dir = dir
input_path = Path(dir)
sc.settings.figdir = figures_dir

if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

repertoire_embedding = args.repertoire_embedding
simulated_embedding = args.simulated_embedding
repertoire = args.repertoire
simulated = args.simulated
target_aa = args.target_aa
fuzziness_param = args.fuzziness_param
#embedding_model = args.embedding_model
abundance = args.abundance
repertoire_sample = args.repertoire_sample
random_seed = args.random_seed


# Input data
input_path = Path(dir)

os.listdir(dir)

# TODO update how to get subject ID
subj_id = str(repertoire).split("/")[-1].split("_subsampled")[0]

# Read input files
logger.info("Reading input files...")
rep_embed = pd.read_csv(repertoire_embedding, sep="\t", header=0)
sim_embed = pd.read_csv(simulated_embedding, sep="\t", header=0)
rep = pd.read_csv(repertoire, sep="\t", header=0)
sim = pd.read_csv(simulated, sep="\t", header=0)

# # Embed target
# def compute_embedding(sequences, model='antiberta2'):
#     if model == 'antiberta2':
#         embedding = amulety.bcr_embeddings.antiberta2(sequences=sequences)
#     elif model == 'esm2':
#         embedding = amulety.bcr_embeddings.esm2(sequences=sequences)
#     elif model == 'antiberty':
#         embedding = amulety.bcr_embeddings.antiberty(sequences=sequences)
#     elif model == 'balm_paired':
#         embedding = amulety.bcr_embeddings.balm_paired(sequences=sequences)
#     else:
#         raise ValueError(f"Unsupported model: {model}")
#     return embedding

# target_embed = compute_embedding(pd.Series([target_aa]), model=embedding_model)

# Select target embedding from simulated sequences embedding
# Get target from simulation (dist = 0)
target_id = sim[sim["dist_to_target"] == 0]["sequence_id"].values[0]
target_embed = sim_embed[sim_embed["sequence_id"] == target_id].drop(columns=["sequence_id"])

logger.info("Repertoire embedded shape: %s", rep_embed.shape)
logger.info("Simulation embedded shape: %s", sim_embed.shape)
logger.info("Repertoire shape: %s", rep.shape)
logger.info("Simulation shape: %s", sim.shape)
logger.info("Processing repertoires...")

# Drop sequence ID duplicates
rep_embed = rep_embed.drop_duplicates(subset=["sequence_id"])
# Get meta of only sequences that are in the embedded data
rep_meta = rep[rep["sequence_id"].isin(rep_embed["sequence_id"])]
rep_meta = rep_meta.drop_duplicates(subset=["sequence_id"])

# Sample repertoire sequences to desired depth
if repertoire_sample < rep_meta.shape[0]:
    logger.info("Sampling %d sequences from repertoire", repertoire_sample)
    random_seed = random_seed if random_seed is not None else 42
    # Sample repertoire metadata
    rep_meta = rep_meta.sample(n=repertoire_sample, random_state=random_seed)
    rep_embed = rep_embed[rep_embed["sequence_id"].isin(rep_meta["sequence_id"])]


rep_meta["subject_id"] = subj_id
rep_meta["sequence_id"] = subj_id + "|" + rep_meta["sequence_id"]
rep_meta["simulated"] = False
rep_embed["sequence_id"] = subj_id + "|" + rep_embed["sequence_id"]

sim_embed = sim_embed.drop_duplicates(subset=["sequence_id"])
sim_embed["sequence_id"] = subj_id + "|" + sim_embed["sequence_id"]

# Deduplicate sim_meta
sim["sequence_id"] = subj_id + "|" + sim["sequence_id"]
sim_meta = sim.drop_duplicates(subset=["sequence_id"])
sim_meta = sim_meta.drop_duplicates(subset=["sequence_vdj_aa"])

# remove sequences with zero distance to target
sim_meta = sim_meta[sim_meta["dist_to_target"] > 0]

sim_meta = sim_meta[sim_meta["sequence_id"].isin(sim_embed["sequence_id"])]

# Take only embeddings from non-duplicated sequences
sim_embed = sim_embed[sim_embed["sequence_id"].isin(sim_meta["sequence_id"])]

# TODO if checking different germlines, we will need to annotate the V gene of the germline
sim_meta["simulated"] = True
sim_meta["subject_id"] = subj_id
sim_meta["v_call"] = "unknown"
sim_meta["j_call"] = "unknown"

# Set sequence_id as index
rep_embed.set_index("sequence_id", inplace=True)
rep_meta.set_index("sequence_id", inplace=True)
sim_embed.set_index("sequence_id", inplace=True)
sim_meta.set_index("sequence_id", inplace=True)

# Reorder rep_meta to match index of rep_embed
rep_meta = rep_meta.reindex(rep_embed.index)
# Reorder sim_meta to match index of sim_embed
sim_meta = sim_meta.reindex(sim_embed.index)

assert rep_meta.index.equals(rep_embed.index), "Indices do not match!"
assert sim_meta.index.equals(sim_embed.index), "Indices do not match!"

logger.info("After filtering:")
logger.info(f"rep_meta shape: {rep_meta.shape}")
logger.info(f"rep_embed shape: {rep_embed.shape}")
logger.info(f"sim_meta shape: {sim_meta.shape}")
logger.info(f"sim_embed shape: {sim_embed.shape}")

### Find centroids
### -----------------------------------------------
logger.info("Finding centroids...")


# Define clusters within repertoire
annd = ad.AnnData(rep_embed)
annd.obs['subject_id'] = rep_meta['subject_id']
annd.obs['v_gene'] = rep_meta['v_call'].str.split("*").str[0]
annd.obs['j_gene'] = rep_meta['j_call'].str.split("*").str[0]
annd.obs['vj_gene'] = annd.obs['v_gene'] + "_" + annd.obs['j_gene']
annd.obs['simulated'] = rep_meta['simulated']
annd.layers['counts'] = annd.X.copy()

sc.tl.pca(annd, n_comps = 50)
sc.pp.neighbors(annd)
sc.tl.umap(annd)

sc.pl.umap(
    annd,
    color=["v_gene"],
    size=2,
    ncols=2,
    show=False,
    frameon=False,
    title=["V gene"],
    save="_clusters_v_gene_scanpy.png",
    #legend_loc="none",
)

resolutions = [round(x * 0.1, 2) for x in range(1, 21)]

for res in resolutions:
    sc.tl.leiden(
        annd, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

res_vars = [key for key in annd.obs.keys() if key.startswith("leiden_res_")]
color_vars = res_vars + ["simulated", "v_gene"]

sc.pl.umap(
    annd,
    color=color_vars,
    palette=list(matplotlib.colors.CSS4_COLORS.values()),
    size=2,
    ncols=3,
    legend_loc="none",
    show=False,
    frameon=False,
    save="_clustering_resolutions.png",
)

# Plot number of clusters per resolution
plt.figure(figsize=(10, 6))
num_clusters = [annd.obs[f'leiden_res_{res:4.2f}'].nunique() for res in resolutions]
plt.plot(resolutions, num_clusters, marker='o')
plt.xlabel("Resolution")
plt.ylabel("Number of Clusters")
plt.title("Number of Clusters per Resolution")
plt.savefig(figures_dir + "num_clusters_per_resolution.png")


# Pick resolution
picked_resolution = 1.0

# Get the cluster centroids for the picked resolution
centroids = sc.get.aggregate(
    annd, by=f"leiden_res_{picked_resolution:4.2f}", func=["mean"]
)
centroids.X = centroids.layers['mean'].copy()

# Plot data with centroids

annd_combined = ad.concat([annd, centroids], axis=0)
annd_combined.layers['counts'] = annd_combined.X.copy()

annd_combined.obs['centroid'] = annd_combined.obs.index.isin(centroids.obs.index)
annd_combined.obs['centroid_color'] = pd.Series(annd_combined.obs['centroid']).map({True: 'centroid', False: 'repertoire'})
annd_combined.obs['centroid_size'] = pd.Series(annd_combined.obs['centroid']).map({True: 50, False: 2})
annd_combined.obs['cluster_color'] = list(annd.obs[f'leiden_res_{picked_resolution:4.2f}']) + list(centroids.obs.index)
annd_combined.obs['alpha'] = pd.Series(annd_combined.obs['centroid']).map({True: 1, False: 0.5})


# Perform PCA and UMAP on combined data with centroids
# Using Sklearn to be able to just transform the simulated sequences
pca = sklearn.decomposition.PCA(n_components=50, random_state=random_seed)
pca.fit(annd_combined.X)

pca_orig = pca.transform(annd_combined.X)

reduce_umap = umap.UMAP(n_neighbors=15,
                        min_dist=0.5,
                        n_components=2,
                        metric='euclidean',
                        random_state=random_seed)

umap_orig = reduce_umap.fit_transform(pca_orig)


annd_combined.obsm["X_umap"] = umap_orig

sc.pl.umap(
    annd_combined,
    color=['centroid_color'],
    alpha=annd_combined.obs['alpha'],
    size=annd_combined.obs['centroid_size'],
    ncols=3,
    #legend_loc="none",
    show=False,
    title=["Clusters with centroids"],
    frameon=False,
    save="_centroids.png",
)

sc.pl.umap(
    annd_combined,
    color=['cluster_color'],
    alpha=annd_combined.obs['alpha'],
    size=annd_combined.obs['centroid_size'],
    ncols=3,
    legend_loc="none",
    show=False,
    title=["Clusters with centroids"],
    frameon=False,
    save="_centroids_cluster_colors.png",
)

# Function to pick simulated sequences based on distance to target and centroids
def pick_simulated_sequences(simulated_sequences_embeddings,
                                centroids_embeddings,
                                target_embedding,
                                abundance,
                                repertoire_sample,
                                fuzziness_param=2):
    """
    Pick sequences from simulated sequences according to their weighted distance to the target.

    Parameters:
    - simulated_sequences_embeddings: DataFrame with simulated sequences embeddings.
    - centroids_embeddings: DataFrame with centroids embeddings.
    - target_embedding: 1D numpy array with target embedding.
    - abundance: float, fraction of sequences to sample.
    - repertoire_sample: int, number of sequences to sample from the repertoire.
    - fuzziness_param: float, parameter for distance weighting. Should take values (1, inf). The higher the value, the fuzzier the clustering is (including sequences further away from the target).

    """
    n_centroids = centroids_embeddings.shape[0]
    logger.info(f"Start simulated sequences selection shape: {simulated_sequences_embeddings.shape}")

    # dist_to_target np 1D array (n_simulated_seqs, )
    dist_to_target = np.linalg.norm(simulated_sequences_embeddings - target_embedding, axis=1)
    logger.info(f"Dist to target shape: {dist_to_target.shape}")

    # dist_to_centroids np 2D array (n_simulated_seqs, n_centroids)
    dist_to_centroids = cdist(simulated_sequences_embeddings, centroids_embeddings, metric='euclidean')
    logger.info(f"Dist to centroids shape: {dist_to_centroids.shape}")


    # dist_to_target is a 1D numpy array of shape (n_simulated_seqs,), repeat to number of centroids
    dist_to_target_tiled = np.tile(dist_to_target[:, np.newaxis], (1, n_centroids))
    logger.info(f"Dist to target tiled shape: {dist_to_target_tiled.shape}")

    # dist_to_target_over_centroids is a 2D numpy array of shape (n_simulated_seqs, n_centroids)
    dist_to_target_over_centroids_sq = (dist_to_target_tiled / dist_to_centroids)**(2/(fuzziness_param-1))


    # compute sum of normalized distances to target over centroids (n_simulated_seqs, )
    dist_to_target_over_centroids_sum = dist_to_target_over_centroids_sq.sum(axis=1)

    # weighted dist to centroids from target (n_simulated_seqs, )
    weighted_dist_to_centroids = 1 / dist_to_target_over_centroids_sum
    logger.info(f"Weighted dist to centroids shape: {weighted_dist_to_centroids.shape}")


    # transformed weighted dist to centroids by logistic function (n_simulated_seqs, )
    weighted_dist_to_centroids_transformed = 1 / (1 + np.exp(-0.5 * weighted_dist_to_centroids))

    # Rescaling weighted distances
    weighted_dist_to_centroids_rescaled = (weighted_dist_to_centroids_transformed - min(weighted_dist_to_centroids_transformed)) / (max(weighted_dist_to_centroids_transformed) - min(weighted_dist_to_centroids_transformed))

    logger.info("Weighted distances to centroids rescaled shape: %s", weighted_dist_to_centroids_rescaled.shape)

    # Pick sequences according to probabilities
    sample_size = abundance * repertoire_sample

    if sample_size >= simulated_sequences_embeddings.shape[0]:
        logger.warning("Sample size is larger than number of simulated sequences. Using all simulated sequences.")
        sample_size = simulated_sequences_embeddings.shape[0]
    logger.info("Sampling %d sequences from simulated sequences", int(sample_size))
    sample = simulated_sequences_embeddings.sample(
        n=int(sample_size),
        weights=weighted_dist_to_centroids_rescaled,
        random_state=random_seed,
        axis=0
    )

    return sample

# transform data to numpy array
logger.info("Selecting simulated sequences...")


target_arr = np.array(target_embed)

simulated_seqs = sim_embed
#sim_seqs_arr = np.array(simulated_seqs)


centr_arr = np.array(centroids.X)
n_centroids = centr_arr.shape[0]

picked_simulated_seqs = pick_simulated_sequences(
    simulated_sequences_embeddings=simulated_seqs,
    centroids_embeddings=centr_arr,
    target_embedding=target_arr,
    abundance=abundance,
    repertoire_sample=repertoire_sample
)

assert(all(simulated_seqs.index == sim_meta.index))


picked_simulated_seqs_meta = sim_meta.loc[picked_simulated_seqs.index]


# Plot distribution of distance to target in original sequences
plt.figure(figsize=(10, 6))
sns.histplot(sim_meta['dist_to_target'], bins=50, kde=False)
plt.xlabel("Distance to Target")
plt.ylabel("Frequency")
plt.title("Distribution of Distance to Target in Original Sequences")
plt.savefig(figures_dir + "dist_to_target_original.png")


# Plot picked simulated meta dist_to_target distribution
plt.figure(figsize=(10, 6))
sns.histplot(picked_simulated_seqs_meta['dist_to_target'], bins=50, kde=False)
plt.xlabel("Distance to Target")
plt.ylabel("Frequency")
plt.title("Distribution of Distances to Target in Picked Simulated Sequences")
plt.savefig(figures_dir + "dist_to_target_selected_sequences.png")


## UMAP of repertoire with all simulated sequences
## -----------------------------------------------
logger.info("Plotting results...")

target_df = pd.DataFrame(target_embed)
target_df.index = pd.Series(["target"])

centroids_df = pd.DataFrame(centroids.X)
centroids_df.index = centroids.obs["leiden_res_1.00"]

centroids_df.columns = simulated_seqs.columns
target_df.columns = simulated_seqs.columns

annd_centroids_simulated = ad.AnnData(pd.concat([rep_embed, simulated_seqs, centroids_df, target_df], axis=0))

annd_centroids_simulated.obs["type"] = ["repertoire"] * rep_embed.shape[0] + ["simulated"] * simulated_seqs.shape[0] + ["centroid"] * centroids_df.shape[0] + ["target"]
annd_centroids_simulated.obs["size"] = [2] * rep_embed.shape[0] + [2] * simulated_seqs.shape[0] + [50] * centroids_df.shape[0] + [50]
annd_centroids_simulated.obs["alpha"] = [0.5] * rep_embed.shape[0] + [0.5] * simulated_seqs.shape[0] + [1.0] * centroids_df.shape[0] + [1.0]

annd_centroids_simulated.layers['counts'] = annd_centroids_simulated.X.copy()

pca_proj_all = pca.transform(annd_centroids_simulated.X)
umap_proj_all = reduce_umap.transform(pca_proj_all)

annd_centroids_simulated.obsm["X_umap"] = umap_proj_all

sc.pl.umap(
    annd_centroids_simulated,
    color=['type'],
    alpha=annd_centroids_simulated.obs['alpha'],
    size=annd_centroids_simulated.obs['size'],
    ncols=3,
    #legend_loc="none",
    show=False,
    title = ["Clusters with centroids and all simulated sequences"],
    frameon=False,
    save="_centroids_with_all_simulated_sequences.png",
)

## UMAP of selected target sequences
annd_centroids_simulated_picked = ad.AnnData(pd.concat([rep_embed, picked_simulated_seqs, centroids_df, target_df], axis=0))

annd_centroids_simulated_picked.obs["type"] = ["repertoire"] * rep_embed.shape[0] + ["simulated"] * picked_simulated_seqs.shape[0] + ["centroid"] * centroids_df.shape[0] + ["target"]
annd_centroids_simulated_picked.obs["size"] = [2] * rep_embed.shape[0] + [1] * picked_simulated_seqs.shape[0] + [50] * centroids_df.shape[0] + [50]
annd_centroids_simulated_picked.obs["alpha"] = [0.5] * rep_embed.shape[0] + [0.5] * picked_simulated_seqs.shape[0] + [1.0] * centroids_df.shape[0] + [1.0]

annd_centroids_simulated_picked.layers['counts'] = annd_centroids_simulated_picked.X.copy()

pca_proj_selected = pca.transform(annd_centroids_simulated_picked.X)
umap_proj_selected = reduce_umap.transform(pca_proj_selected)

annd_centroids_simulated_picked.obsm["X_umap"] = umap_proj_selected

sc.pl.umap(
    annd_centroids_simulated_picked,
    color=['type'],
    alpha=annd_centroids_simulated_picked.obs['alpha'],
    size=annd_centroids_simulated_picked.obs['size'],
    ncols=3,
    #legend_loc="none",
    show=False,
    frameon=False,
    title = ["Target and picked simulated sequences"],
    save="_centroids_with_picked_simulated_sequences.png",
)

logger.info("Saving results...")

### Save results
repertoire_with_simulated_meta = pd.concat([rep_meta, picked_simulated_seqs_meta])
repertoire_with_simulated_meta.reset_index(inplace=True, names="sequence_id")
repertoire_with_simulated_meta.to_csv(dir+subj_id+"_repertoire_with_simulated_meta.tsv", sep="\t", index=False)

repertoire_with_simulated_embed = pd.concat([rep_embed, picked_simulated_seqs])
repertoire_with_simulated_embed.reset_index(inplace=True, names="sequence_id")
repertoire_with_simulated_embed.to_csv(dir+subj_id+"_repertoire_with_simulated_embedding.tsv", sep="\t", index=False)
logger.info("Done!")
