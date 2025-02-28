from random import choices
import random
from collections import Counter
import numpy as np
import argparse # Import the argparse module

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



def main():
    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(description="Simulate SNP data with varying parameters.") # Creates an ArgumentParser object which will hold the information necessary to parse the command line into Python data types.
    parser.add_argument('--n_snps', type=int, default=1000, help='Number of SNPs to simulate (integer >= 1)') # Adds argument for n_snps, with type int, default 1000 and help message.
    parser.add_argument('--p', type=int, default=6, help='Ploidy (integer >= 1)') # Adds argument for ploidy
    parser.add_argument('--f', type=int, default=0, help='Correction factor (positive integer)') # Adds argument for correction factor
    parser.add_argument('--d', type=str, default="poisson", choices=['poisson'], help='Sequence depth distribution (default: poisson)') # Adds argument for distribution, with choices limited to 'poisson' for now.
    parser.add_argument('--lam', type=int, default=10, help='Lambda for Poisson distribution (integer >= 1)') # Adds argument for lambda

    args = parser.parse_args() # Parses the arguments and returns an object containing the arguments.

    # --- Input Validation ---
    if args.n_snps < 1: # Checks if n_snps is less than 1
        raise ValueError("n_snps must be an integer greater than or equal to 1") # Raises ValueError if condition is met
    if args.p < 1: # Checks if p is less than 1
        raise ValueError("p must be an integer greater than or equal to 1") # Raises ValueError if condition is met
    if args.f < 0: # Checks if f is less than 0 (not positive)
        raise ValueError("f must be a positive integer") # Raises ValueError if condition is met
    if args.lam < 1: # Checks if lam is less than 1
        raise ValueError("lam must be an integer greater than or equal to 1") # Raises ValueError if condition is met


    n_snps = args.n_snps # Assigns the parsed n_snps value
    p = args.p # Assigns the parsed p value
    f = args.f # Assigns the parsed f value
    d = args.d # Assigns the parsed d value
    lam = args.lam # Assigns the parsed lam value


    snp_num_prob = snp_num_prob_gen(p, f) # Calls snp_num_prob_gen with p and f from arguments

    snp_num_counts= Counter(choices(list(snp_num_prob.keys()),weights=snp_num_prob.values(), k = n_snps)) # Uses n_snps from arguments

    print(snp_num_counts)

    # Generate ploidy related class frequency distributions

    obs_freqs = get_observable_frequencies(p) # Uses p from arguments


    site_id = 0

     # Generate depth distributions per each snp num
        # Iterate over each of the previous distributions
            # Select a random possible state from the class frequencies
            # Select snp_num random bases for the SNPs (if 4 >= snp_num >= 2 at least snp_num bases must be different), if snp_num >4 atleast 4 bases must be different
            # Sample bases using as weights the  possible state from class frequencies
            # With the counts from the sampling:
                # Print:
                    # The site number, for each snp one site
                    # The Base that was sampled
                    # The counts of that base
                    # The depth of the site


    for snp_num,count in snp_num_counts.items():
        for depth in sample_truncated_poisson(lam, snp_num, count): # Uses lam from arguments
            random_snp_weights = choices(obs_freqs[snp_num-2], k=1)[0]
            dna_bases_site = sample_dna_bases(snp_num)
            #print(dna_bases_site,random_snp_weights)
            dna_bases_counts = Counter(choices(dna_bases_site,weights=random_snp_weights, k = depth))
            while len(set(dna_bases_counts.keys())) < len(set(dna_bases_site)):
                dna_bases_counts = Counter(choices(dna_bases_site,weights=random_snp_weights, k = depth))
            for base,count in dna_bases_counts.items():
                print(site_id,snp_num,base,count,depth,sep="\t")
            site_id += 1


if __name__ == "__main__":
    main()