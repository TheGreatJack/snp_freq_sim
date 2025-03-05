from random import choices
import random
from collections import Counter
import numpy as np
import argparse
import hashlib # For generating unique hashes
import itertools # For creating combinations
import os # For file path operations

def get_observable_frequencies(n):
    """
    Generates a list of lists of lists describing observable class frequencies for n samplings from n items.

    Args:
        n: The number of samplings and total items.

    Returns:
        A list of lists of lists. Each inner list of lists corresponds to a number of classes 'c' (from 2 to n).
        Each list within that is a sorted list of observable class frequencies for that 'c'.
    """
    result_list = []

    def get_partitions(remaining_sum, parts_needed, current_partition):
        if parts_needed == 0:
            if remaining_sum == 0:
                yield current_partition[:]  # Yield a copy of the current partition
            return
        for i in range(1, remaining_sum - (parts_needed - 1) + 1):
            if i > 0:
                current_partition.append(i)
                yield from get_partitions(remaining_sum - i, parts_needed - 1, current_partition)
                current_partition.pop()  # Backtrack


    for c in range(2, n + 1):

        unique_frequency_sets = set()
        partitions = get_partitions(n, c, [])
        partitions = [p for p in partitions]

        for partition in partitions:
            frequencies = [round(p / n, 5) for p in partition] # Round to avoid floating point issues, 5 decimal places should be enough
            frequencies.sort()
            unique_frequency_sets.add(tuple(frequencies))

        unique_frequency_lists = [list(freq_tuple) for freq_tuple in unique_frequency_sets]
        result_list.append(unique_frequency_lists)

    return result_list


def snp_num_prob_gen(p, factor):
    probs = {2:1}
    for i in range(p,2,-1):
        probs[i] = 10 ** -(i+factor-2)

    total_prob = sum(probs.values())

    for i in probs.keys():
        probs[i] = probs[i] / total_prob

    return probs


def sample_truncated_poisson(lam, minimum, n):
    """
    Samples n numbers from a Poisson distribution with parameter lambda,
    with the constraint that no sampled number is lower than minimum.

    Parameters:
    lam (float or int): Lambda parameter of the Poisson distribution (rate parameter, must be non-negative).
    minimum (int): The minimum value allowed for the sampled numbers (must be a non-negative integer).
    n (int): The number of samples to generate.

    Returns:
    list: A list of n integers sampled from the truncated Poisson distribution.
    """
    if lam < 0:
        raise ValueError("Lambda parameter must be non-negative.")
    if minimum < 0 or not isinstance(minimum, int):
        raise ValueError("Minimum parameter must be a non-negative integer.")
    if n <= 0 or not isinstance(n, int):
        raise ValueError("Number of samples n must be a positive integer.")

    samples = []# Generate ploidy related class frequency distributions
    while len(samples) < n:
        sample = np.random.poisson(lam)
        if sample >= minimum:
            samples.append(sample)
    return samples

def sample_dna_bases(n):
    """
    Samples n DNA bases with diversity constraints.

    Args:
        n (int): The number of DNA bases to sample.

    Returns:
        list: A list of n DNA bases (strings 'A', 'T', 'C', 'G').

    Raises:
        ValueError: If n is not a positive integer or if n < 2.
    """
    # Input validation
    if not isinstance(n, int) or n <= 0: # Corrected to be n <= 0 as the constraint starts from n>=2
        raise ValueError("n must be a positive integer greater than or equal to 2.")
    if n < 2:
        raise ValueError("n must be greater than or equal to 2 to satisfy the diversity constraints.")


    bases = ['A', 'T', 'C', 'G']

    if 2 <= n <= 4:
        # Case 1: All bases must be different
        sampled_bases = random.sample(bases, n)
    elif n > 4:
        # Case 2: At least 4 bases must be different
        distinct_bases = random.sample(bases, 4)
        remaining_bases_count = n - 4
        additional_bases = random.choices(bases, k=remaining_bases_count) # Use choices for repetition
        sampled_bases = distinct_bases + additional_bases
        random.shuffle(sampled_bases) # Shuffle for randomness
    else: # Should not reach here due to validation, but just in case.
        raise ValueError("n value not handled by function logic (should be >= 2).") # More descriptive error

    return sampled_bases

def parse_comma_separated_int(value):
    """
    Parses a comma-separated string of integers into a list of integers.
    Raises argparse.ArgumentTypeError if parsing fails or values are not positive integers.
    """
    try:
        int_list = [int(x) for x in value.split(',')] # Splits the string by comma and converts each part to int
        for i in int_list: # Checks each integer in the list
            if i < 1: # If any integer is less than 1
                raise argparse.ArgumentTypeError(f"Values must be positive integers, but got: {value}") # Raise error if not positive
        return int_list # Returns the list of integers
    except ValueError: # Catches errors during int conversion
        raise argparse.ArgumentTypeError(f"Invalid comma-separated integers: {value}") # Raise error if conversion fails

def generate_output_filename(base_hash, n_snps, p, f, lam, rep_num = None, data_type="raw_data"):
    """
    Generates a unique output filename based on parameters and a hash.
    """
    filename = f"{base_hash}_n_snps_{n_snps}_p_{p}_f_{f}_lam_{lam}" # Base filename with hash and parameters
    if rep_num is not None: # If repetition number is provided, add it to the filename
        filename += f"_rep_{rep_num}" # Adds repetition number
    filename += f"_{data_type}.tsv" # Adds data type and extension
    return filename # Returns the generated filename

def calculate_hash(n_snps, p, f, lam):
    """
    Calculates an MD5 hash from the combination of parameters to ensure unique filenames.
    """
    param_str = f"n_snps_{n_snps}_p_{p}_f_{f}_lam_{lam}" # Creates a string from parameters
    hash_object = hashlib.md5(param_str.encode()) # Encodes the string and creates MD5 hash object
    return hash_object.hexdigest() # Returns the hexadecimal representation of the hash

def main():
    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(description="Simulate SNP data for multiple parameter combinations and repetitions.")
    parser.add_argument('--n_snps', type=parse_comma_separated_int, default="1000", help='Number of SNPs to simulate (integer >= 1), comma-separated list for multiple values')
    parser.add_argument('--p', type=parse_comma_separated_int, default="6", help='Ploidy (integer >= 1), comma-separated list for multiple values')
    parser.add_argument('--f', type=parse_comma_separated_int, default="0", help='Correction factor (positive integer), comma-separated list for multiple values')
    parser.add_argument('--d', type=str, default="poisson", choices=['poisson'], help='Sequence depth distribution (default: poisson)')
    parser.add_argument('--lam', type=parse_comma_separated_int, default="10", help='Lambda for Poisson distribution (integer >= 1), comma-separated list for multiple values')
    parser.add_argument('--reps', type=int, default=1, help='Number of repetitions per parameter combination (positive integer)')
    parser.add_argument('--k_bins', type=int, default=10, help='Number of bins for the histogram of observed proportions (positive integer)')


    args = parser.parse_args()

    # --- Input Validation ---
    if args.reps < 1:
        raise ValueError("reps must be a positive integer")
    if args.k_bins < 1:
        raise ValueError("k_bins must be a positive integer")


    n_snps_values = args.n_snps # List of n_snps values from CLI
    p_values = args.p # List of p values from CLI
    f_values = args.f # List of f values from CLI
    d = args.d # Distribution type, single value
    lam_values = args.lam # List of lam values from CLI
    reps = args.reps # Number of repetitions, single value
    k_bins = args.k_bins # Number of histogram bins, single value


    # --- Parameter Combinations and Simulation ---
    param_combinations = itertools.product(n_snps_values, p_values, f_values, lam_values) # Creates all combinations of parameters

    for n_snps, p, f, lam in param_combinations: # Loops through each parameter combination
        base_hash = calculate_hash(n_snps, p, f, lam) # Calculate hash for the parameter combination
        all_obs_props = [] # List to store all observed proportions for histogram generation

        for rep_num in range(1, reps + 1): # Loops through repetitions
            output_filename_raw = generate_output_filename(base_hash, n_snps, p, f, lam, rep_num, "raw_data") # Generates filename for raw data
            print(f"Running simulation for n_snps={n_snps}, p={p}, f={f}, lam={lam}, rep={rep_num}") # Prints current simulation parameters

            snp_num_prob = snp_num_prob_gen(p, f) # Generate SNP number probabilities
            snp_num_counts= Counter(choices(list(snp_num_prob.keys()),weights=snp_num_prob.values(), k = n_snps)) # Count of each snp number

            obs_freqs = get_observable_frequencies(p) # Get observable frequencies for ploidy p

            site_id = 0
            current_run_obs_props = [] # List to store observed proportions for the current run
            output_lines = [] # List to store lines to write to the raw data output file

            for snp_num,count in snp_num_counts.items(): # Iterate through each snp number and its count
                for depth in sample_truncated_poisson(lam, snp_num, count): # Sample depth for each snp based on poisson distribution
                    random_snp_weights = choices(obs_freqs[snp_num-2], k=1)[0] # Choose random weights for the snp
                    dna_bases_site = sample_dna_bases(snp_num) # Sample DNA bases for the site
                    dna_bases_counts = Counter(choices(dna_bases_site,weights=random_snp_weights, k = depth)) # Sample base counts based on weights and depth
                    while len(set(dna_bases_counts.keys())) < len(set(dna_bases_site)): # Ensure diversity of bases matches snp_num
                        dna_bases_counts = Counter(choices(dna_bases_site,weights=random_snp_weights, k = depth)) # Resample if diversity is not met
                    for base, count in dna_bases_counts.items(): # Iterate through each base and its count
                        obs_prop = count / depth # Calculate observed proportion
                        output_line = f"{site_id}\t{snp_num}\t{base}\t{count}\t{depth}\t{obs_prop:.5f}" # Format output line with observed proportion
                        output_lines.append(output_line) # Add line to output lines
                        current_run_obs_props.append(obs_prop) # Add observed proportion to current run list
                        #print(output_line) # Print raw data to standard output (optional, for debugging)
                        site_id += 1 # Increment site ID

            all_obs_props.append(current_run_obs_props) # Extend the list of all observed proportions with the current run's proportions

            # --- Save Raw Data to File ---
            with open(output_filename_raw, 'w') as outfile: # Open file for writing raw data
                outfile.write("site_id\tsnp_num\tbase\tcount\tdepth\tobs_prop\n") # Write header
                outfile.write("\n".join(output_lines)) # Write all output lines to the file
            print(f"Raw data saved to: {output_filename_raw}") # Print message indicating raw data file saved

        # --- Generate and Save Histogram ---
        
        output_filename_hist = generate_output_filename(base_hash, n_snps, p, f, lam, data_type="histograms") # Generate filename for histograms

        with open(output_filename_hist, 'w') as outfile: # Open file for writing histograms
            for data_list in all_obs_props: # Iterate through histogram bins
                hist, bin_edges = np.histogram(data_list, bins=k_bins, range=(0.0, 1.0), density = False) # Generate histogram of observed proportions
                hist = hist/len(data_list)
                line_to_write = ""
                for i, num in enumerate(hist):
                    formatted_num = "{:.5f}".format(float(num)) # Format to 5 decimal places
                    line_to_write += formatted_num
                    if i < len(hist) - 1: # Add tab if not the last number
                        line_to_write += "\t"
                outfile.write(line_to_write + "\n") # Add newline at the end of the line
    
        print(f"Histograms saved to: {output_filename_hist}") # Print message indicating histogram file saved

if __name__ == "__main__":
    main()