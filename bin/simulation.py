# !/usr/bin/env python
# This script simulates the convergence of antibody sequences towards a target sequence using a naive stereotyping approach.
# It generates child sequences from a starting sequence and evaluates their distance to a target sequence.
# Usage: python simulation.py --input_repertoire <input_file.tsv> --target_aa <target_sequence.tsv> --embedding_model_input <model_name> --distance_function_input <distance_function> --n_children_per_generation <number_of_children> --mutation_rate <mutation_rate> --output_prefix <output_prefix>
# The input repertoire should be a TSV file with columns for nucleotide and amino acid sequences. It requires the following columns:
# - sequence: Nucleotide sequence
# - sequence_aa: Amino acid sequence
# - v_call: V gene call
# - j_call: J gene call
# - junction_aa: Junction amino acid sequence
# - sequence_id: Unique identifier for each sequence
# - cdr1, fwr2, cdr2, fwr3, cdr3, fwr4: CDR and FWR regions of the sequence



import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import r
import matplotlib.pyplot as plt
import logomaker
import amulety
import torch
import warnings
import argparse
import os
warnings.filterwarnings("ignore", module="amulety.amulety")

# Read simulated repertoire

# Repertoire values

parser = argparse.ArgumentParser(description="Simulate convergence of antibody sequences.")
parser.add_argument("--input_repertoire", type=str, help="Path to input repertoire TSV file")
parser.add_argument("--target_aa", type=str, help="Path to target sequence TSV file")
parser.add_argument("--embedding_model_input", type=str, default="None", help="Embedding model to use")
parser.add_argument("--distance_function_input", type=str, default="hamming_distance", help="Distance function to use")
parser.add_argument("--n_children_per_generation", type=int, default=100, help="Number of children per generation")
parser.add_argument("--mutation_rate", type=float, default=0.001, help="Mutation rate per generation")
parser.add_argument("--output_prefix", type=str, default=None, help="Prefix for output files")
parser.add_argument("--sequence_nt_column", type=str, default="sequence", help="Column name for nucleotide sequence")
parser.add_argument("--sequence_aa_column", type=str, default="sequence_vdj_aa", help="Column name for amino acid sequence")
parser.add_argument("--random_seed", type=int, default=42, help="Random seed for reproducibility")

args = parser.parse_args()

input_repertoire = args.input_repertoire
sequence_nt_column = args.sequence_nt_column
sequence_aa_column = args.sequence_aa_column
target_aa = args.target_aa
embedding_model_input = args.embedding_model_input
distance_function_input = args.distance_function_input
n_children_per_generation = args.n_children_per_generation
mutation_rate = args.mutation_rate

if args.output_prefix is None:
    output_prefix = os.path.basename(input_repertoire).replace('.tsv', '') + "_simulation_with_target_" + target_aa
else:
    output_prefix = args.output_prefix

# Enable automatic pandas conversion if needed
pandas2ri.activate()

# Import the shazam package
shazam = importr('shazam')
alakazam = importr('alakazam')


# %% [markdown]
# ## Functions

# %%
class Sequence:
    def __init__(self, sequence_aa, sequence_id, sequence_nt, dist_to_target=None, embedding=None):
        self.sequence_aa = sequence_aa
        self.sequence_id = sequence_id
        self.sequence_nt = sequence_nt
        self.dist_to_target = dist_to_target
        self.embedding = embedding

def hamming_distance(sequences: pd.Series,
                    target: pd.Series,
                    normalized=False):
    """
    Calculate the Hamming distance between two strings of equal length.
    Hamming distance is the number of positions at which the corresponding symbols are different.
    Args:
        str1 (str): First string.
        str2 (str): Second string.
        normalized (bool): Whether to normalize the distance by the length of the strings.
    Returns:
        int or float: Hamming distance between the two strings, normalized if specified.
    """
    def hdist(str1, str2):
        return sum(c1 != c2 for c1, c2 in zip(str1, str2))

    # Check equal length
    seqs_len = sequences.str.len().values
    target_len = target.str.len().values
    if not all(seqs_len == target_len):
        raise ValueError("All sequences must be of equal length to the target for Hamming distance calculation.")

    dist = sequences.apply(lambda x: hdist(x, target.iloc[0])).values
    if normalized:
        return dist / target_len
    else:
        return dist

def compute_embedding(sequences, model='antiberta2'):
    if model == 'antiberta2':
        embedding = amulety.amulety.antiberta2(sequences=sequences)
    elif model == 'esm2':
        embedding = amulety.amulety.esm2(sequences=sequences)
    elif model == 'antiberty':
        embedding = amulety.amulety.antiberty(sequences=sequences)
    elif model == 'balm_paired':
        embedding = amulety.amulety.balm_paired(sequences=sequences)
    else:
        raise ValueError(f"Unsupported model: {model}")
    return embedding

class Stereotyping:
    """
    A class to perform naive stereotyping of sequences based on a target sequence.
    This class generates child sequences from a starting sequence and evaluates their distance to a target sequence.
    Args:
        target_aa (str): The target amino acid sequence to compare against.
        start_nt (str): The starting nucleotide sequence.
        start_aa (str): The starting amino acid sequence.
        targeting_model (str): The model used for generating child sequences.
        distance_function (callable): Function to calculate distance between sequences.
        distance_normalized (bool): Whether to normalize the distance.
    Attributes:
        distance_function (callable): Function to calculate distance between sequences.
        distance_normalized (bool): Whether to normalize the distance.
        start (Sequence): The starting sequence object.
        target (Sequence): The target sequence object.
        selected (Sequence): The currently selected sequence based on distance to target.
        n_generations (int): Number of generations processed.
        n_children (int): Total number of children generated.
        targeting_model (str): The model used for generating child sequences.
        generations (dict): Dictionary to store sequences generated in each generation.
    """

    def find_closest_seq_to_target(
                    self,
                    sequences):
        """
        Find the closest sequence in compatible_naives pd.Series to the target_aa sequence.
        Returns the closest row and its minimum distance.
        """
        if not isinstance(sequences, pd.Series):
            raise TypeError("compatible_naives must be a pandas Series.")

        if self.embedding_model is not None:
            embeddings = compute_embedding(
                sequences=sequences,
                model=self.embedding_model
            )
            # Compute euclidean distances
            distances = torch.cdist(embeddings, self.target.embedding)
            min_index = int(distances.argmin())
            return embeddings, distances.numpy().flatten(), min_index

        else:
            distances = self.distance_function(
                sequences, # Convert to string for distance calculation
                pd.Series([self.target.sequence_aa]),
                normalized=self.distance_normalized
            )
            min_index = int(distances.argmin())
            return "", distances, min_index

    def __init__(self,
                target_aa,
                compatible_naives,
                targeting_model,
                distance_function=hamming_distance,
                distance_normalized=True,
                sequence_nt_column="sequence",
                sequence_aa_column="sequence_aa",
                embedding_model=None):
        self.distance_function = distance_function
        self.distance_normalized = distance_normalized
        self.embedding_model = embedding_model
        self.target = Sequence(
            sequence_aa=target_aa,
            sequence_id="target_sequence",
            sequence_nt=None,  # Target sequence does not have a nucleotide representation
            dist_to_target=0,  # Distance to itself is zero
            embedding=compute_embedding(pd.Series([target_aa]), model=embedding_model) if embedding_model else None
        )

        # Select closest sequence from compatible_naives
        if not isinstance(compatible_naives, pd.DataFrame):
            raise TypeError("compatible_naives must be a pandas DataFrame.")
        if sequence_nt_column not in compatible_naives.columns or sequence_aa_column not in compatible_naives.columns:
            raise ValueError(f"DataFrame must contain columns '{sequence_nt_column}' and '{sequence_aa_column}'.")

        # Find the closest sequence to the target in the compatible_naives DataFrame
        (embeddings, distances, closest_index) = self.find_closest_seq_to_target(compatible_naives[sequence_aa_column])

        self.start = Sequence(
            sequence_aa=compatible_naives[sequence_aa_column].iloc[closest_index],
            sequence_id=compatible_naives["sequence_id"].iloc[closest_index],
            sequence_nt=compatible_naives[sequence_nt_column].iloc[closest_index].upper(),
            dist_to_target=distances[closest_index],
            embedding=embeddings[closest_index] if self.embedding_model else None
        )
        print(f"Selected start sequence: {self.start.sequence_id} with distance {self.start.dist_to_target:.4f} to target.")
        print (f"Start sequence: {self.start.sequence_aa}")

        self.selected = dict()
        self.selected[0] = self.start
        self.n_generations = 0
        self.n_children = 0
        self.targeting_model = targeting_model
        self.generations = dict()

    def get_children(self, n_children, mutation_rate=0.001, include_parent=True):
        """
        Generate child sequences based on the selected sequence.
        Args:
            n_children (int): Number of child sequences to generate.
            mutation_rate (float): Rate of mutations to apply to the selected sequence.
            include_parent (bool): Whether to include the parent sequence in the children list.
        Returns:
            list: A list of Sequence objects representing the generated children.
        """
        children_nt = []
        children_aa = []
        children_id = []
        for i in range(n_children):
            child_nt = str(shazam.shmulateSeq(
                    self.selected[self.n_generations].sequence_nt,
                    targetingModel=self.targeting_model,
                    numMutations=mutation_rate,
                    frequency=True
                    )[0])
            child_aa = str(alakazam.translateDNA(child_nt)[0])
            # Ensure the child sequence has same length as the target sequence
            if len(child_aa) != len(self.target.sequence_aa):
                print("Skipping child with stop codons.")
                print(f"Child {i} sequence: {child_aa}")
                print(f"Target sequence: {self.target.sequence_aa}")
                next

            children_nt.append(child_nt)
            children_aa.append(child_aa)
            children_id.append(f"generation_{self.n_generations + 1}_child_{i}")

        # Include the parent sequence in the children list
        if include_parent:
            children_nt.append(self.selected[self.n_generations].sequence_nt)
            children_aa.append(self.selected[self.n_generations].sequence_aa)
            children_id.append(self.selected[self.n_generations].sequence_id)
            n_children += 1  # Increment the number of children to account for the parent

        (embeddings, distances, min_dist) = self.find_closest_seq_to_target(pd.Series(children_aa))

        # Create Sequence objects for each child
        children = [
            Sequence(
                sequence_aa=children_aa[i],
                sequence_id=children_id[i],
                sequence_nt=children_nt[i],
                dist_to_target=distances[i],
                embedding=embeddings[i] if self.embedding_model else None
            )
            for i in range(n_children)
        ]

        # Update the generation count
        self.n_generations += 1

        # Store the generated children in the generations dictionary
        self.generations[self.n_generations] = children
        # Include the parent sequence in the children list

        # Update the generations and counters
        self.n_children += n_children
        # Update the selected sequence
        self.selected[self.n_generations] = children[min_dist]
        return children

    def iterate_to_target(self, n_children_per_generation, mutation_rate=0.001, include_parent=True):
        """
        Iterate until the selected sequence is within a specified distance to the target.
        This method generates children until a sequence is found that meets the distance criteria.
        """
        if self.n_generations > 0:
            warnings.warn(f"Current model has already a generation. Starting from the current generation: {self.n_generations}")
        while self.selected[self.n_generations].dist_to_target > 0:
            _children = self.get_children(n_children_per_generation, mutation_rate, include_parent)
            print(f"Generation {self.n_generations}: Selected sequence {self.selected[self.n_generations].sequence_id} with distance {self.selected[self.n_generations].dist_to_target:.4f} to target.")
            if self.selected[self.n_generations].dist_to_target == 0:
                print(f"Target sequence reached: {self.selected[self.n_generations].sequence_id}")
                break
        return self.generations, self.n_generations, self.n_children

    def save_generations(self, filename):
        """
        Write the generations of sequences to a CSV file.
        Args:
            filename (str): The name of the file to write the generations to.
        """
        if not filename.endswith('.tsv'):
            raise ValueError("Filename must end with '.tsv'.")

        all_sequences = []
        for gen, children in self.generations.items():
            for child in children:
                all_sequences.append({
                    'generation': gen,
                    'sequence_id': child.sequence_id,
                    'sequence_aa': child.sequence_aa,
                    'sequence_nt': child.sequence_nt,
                    'dist_to_target': child.dist_to_target
                })

        df = pd.DataFrame(all_sequences)
        df.to_csv(filename, index=False, sep='\t')

    def save_selected(self, filename):
        """
        Write all the selected sequencs to a CSV file.
        Args:
            filename (str): The name of the file to write the selected sequence to.
        """
        if not filename.endswith('.tsv'):
            raise ValueError("Filename must end with '.tsv'.")

        selected_sequences = []
        for gen, seq in self.selected.items():
            selected_sequences.append({
                'generation': gen,
                'sequence_id': seq.sequence_id,
                'sequence_aa': seq.sequence_aa,
                'sequence_nt': seq.sequence_nt,
                'dist_to_target': seq.dist_to_target
            })

        df = pd.DataFrame(selected_sequences)
        df.to_csv(filename, index=False, sep='\t')


    def plot_distance_generations(self, filename='generations_plot.png'):
        """
        Plot the distance to target for all children in each generation.
        This method uses matplotlib to visualize the distance of the selected sequence to the target over generations.
        """
        if not self.generations:
            raise ValueError("No generations to plot. Run iterate_to_target first.")

        distances = []
        generations = []

        for gen, children in self.generations.items():
            for child in children:
                distances.append(child.dist_to_target)
                generations.append(gen)

        # Prepare DataFrame for seaborn boxplot
        df_plot = pd.DataFrame({
            "generation": generations,
            "distance": distances
        })
        # Count number of sequences for each distance per generation
        count_df = df_plot.groupby(['generation', 'distance']).size().reset_index(name='count')
        plt.clf()

        plt.figure(figsize=(14, 6))
        pivot = count_df.pivot(index='generation', columns='distance', values='count').fillna(0)
        for distance in pivot.columns:
            plt.plot(pivot.index, pivot[distance], marker='o', markersize=1, label=f'Distance {distance}')

        plt.title('Number of Sequences per Distance per Generation')
        plt.xlabel('Generation')
        plt.ylabel('Number of Sequences')
        plt.legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')

    def plot_logo_generations(self, filename='generations_logo.png'):
        """
        Plot a sequence logo for all the generated sequences.
        """
        all_sequences = []
        for gen, children in self.generations.items():
            for child in children:
                all_sequences.append({
                    'generation': gen,
                    'sequence_id': child.sequence_id,
                    'sequence_aa': child.sequence_aa,
                    'sequence': child.sequence_nt,
                    'dist_to_target': child.dist_to_target
                })

        df = pd.DataFrame(all_sequences)
        counts_mat = logomaker.alignment_to_matrix(df['sequence_aa'].tolist())

        # Plot the sequence logo
        plt.figure(figsize=(25, 2))
        logo = logomaker.Logo(counts_mat, shade_below=.5, fade_below=.5, font_name='Arial')
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.ax.set_ylabel('Count')
        plt.savefig(filename, dpi=300, bbox_inches='tight')



# Get the HH_S5F model from R
hh_s5f = r["HH_S5F"]

# Read embedding model param
if embedding_model_input != "None":
    if embedding_model_input not in ["antiberta2", "esm2", "antiberty", "balm_paired"]:
        raise ValueError(f"Unsupported embedding model: {embedding_model_input}. Supported models are 'antiberta2', 'esm2', 'antiberty', 'balm_paired'.")
    embedding_model = embedding_model_input
else:
    embedding_model = None

# Read distance function
if distance_function_input == "hamming_distance":
    distance_function = hamming_distance
elif embedding_model_input == "None":
    raise ValueError(f"Unsupported distance function: {distance_function_input}. Only 'hamming_distance' is supported.")

# Read repertoire
repertoire = pd.read_csv(input_repertoire, sep="\t", header=0)

# Preprocess repertoire
repertoire["junction_length_aa"] = repertoire["junction_aa"].apply(lambda x: len(x))
repertoire["CDRH3_length_aa"] = repertoire["junction_length_aa"] - 2
repertoire["sequence_aa_length"] = repertoire["sequence_vdj_aa"].apply(lambda x: len(x))

# Read target sequence
len_target_aa = len(target_aa)

#TODO : use only unmutated for compatible, and consider all sequences of same length
# Get compatible naive sequences
compatible_naives_def = (repertoire["sequence_aa_length"] == len_target_aa)
repertoire_compatible = repertoire[compatible_naives_def]
print(f"Found {repertoire_compatible.shape[0]} compatible naive sequences for target {target_aa} with length {len_target_aa}.")

# Get distance
compatible_naives = repertoire_compatible.copy()
targeting_model = hh_s5f


sequences = compatible_naives[sequence_aa_column]


stereo = Stereotyping(
    target_aa=target_aa,
    compatible_naives=repertoire_compatible,
    targeting_model=hh_s5f,
    sequence_nt_column=sequence_nt_column,
    sequence_aa_column=sequence_aa_column,
    embedding_model=embedding_model,
    distance_function=hamming_distance,
    distance_normalized=False
)

stereo.start.dist_to_target


all_children = stereo.iterate_to_target(
    n_children_per_generation=n_children_per_generation,
    mutation_rate=mutation_rate,
    include_parent=True
)

stereo.save_generations(f"{output_prefix}_generations.tsv")
stereo.save_selected(f"{output_prefix}_selected_sequences.tsv")
stereo.plot_distance_generations(f"{output_prefix}_generations_plot.png")
stereo.plot_logo_generations(f"{output_prefix}_generations_logo.png")
