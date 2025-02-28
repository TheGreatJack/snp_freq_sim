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
        for partition in partitions:
            frequencies = [round(p / n, 5) for p in partition] # Round to avoid floating point issues, 5 decimal places should be enough
            frequencies.sort()
            unique_frequency_sets.add(tuple(frequencies))

        unique_frequency_lists = [list(freq_tuple) for freq_tuple in unique_frequency_sets]
        result_list.append(unique_frequency_lists)

    return result_list

# Example usage and verification
n_value_3 = 3
result_3 = get_observable_frequencies(n_value_3)
print(f"Result for n={n_value_3}: {result_3}")
# Expected for n=3: [[[0.33333, 0.66667]], [[0.33333, 0.33333, 0.33333]]]
# Output will be slightly rounded

n_value_4 = 4
result_4 = get_observable_frequencies(n_value_4)
print(f"Result for n={n_value_4}: {result_4}")
# Expected for n=4: [[[0.25, 0.75], [0.5, 0.5]], [[0.25, 0.25, 0.5]], [[0.25, 0.25, 0.25, 0.25]]]
# Output will be slightly rounded

# To get exactly the format in the example, we can format the output numbers.
def format_result(result):
    formatted_result = []
    for c_list in result:
        formatted_c_list = []
        for freq_list in c_list:
            formatted_freq_list = [str(f) for f in freq_list]
            formatted_c_list.append(formatted_freq_list)
        formatted_result.append(formatted_c_list)
    return formatted_result

# To get the fractions as in the original example (approximately), we can use fractions.
from fractions import Fraction

def get_observable_frequencies_fraction(n):
    result_list = []

    def get_partitions(remaining_sum, parts_needed, current_partition):
        if parts_needed == 0:
            if remaining_sum == 0:
                yield current_partition[:]
            return
        for i in range(1, remaining_sum - (parts_needed - 1) + 1):
            if i > 0:
                current_partition.append(i)
                yield from get_partitions(remaining_sum - i, parts_needed - 1, current_partition)
                current_partition.pop()

    for c in range(2, n + 1):
        unique_frequency_sets = set()
        partitions = get_partitions(n, c, [])
        for partition in partitions:
            frequencies = [Fraction(p, n) for p in partition]
            frequencies.sort()
            unique_frequency_sets.add(tuple(frequencies))

        unique_frequency_lists = [list(freq_tuple) for freq_tuple in unique_frequency_sets]
        result_list.append(unique_frequency_lists)
    return result_list


result_3_fraction = get_observable_frequencies_fraction(n_value_3)
print(f"Fraction Result for n={n_value_3}: {result_3_fraction}")

result_4_fraction = get_observable_frequencies_fraction(n_value_4)
print(f"Fraction Result for n={n_value_4}: {result_4_fraction}")