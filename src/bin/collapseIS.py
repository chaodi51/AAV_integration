import pandas as pd
import ast
from collections import Counter
import random
import sys

pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.max_colwidth', None) 

infile, outfile = sys.argv[1:]

col = ['sample', 'vector', 'chrom', 'IS_pos', 'vector_align_clip_pos', 'nReads','sonication_sites']
# df = pd.read_csv('~/aav_integration/v1.0.0/data/interim/IS_results/IS_events_anno/shortRead/240112_VH00163_77_AAF5V37M5/012_S11_S11_a1__CAG_EGFP_ITR_corrected__IS.csv', usecols = col)
# df = pd.read_csv('~/aav_integration/v1.0.0/data/interim/IS_results/IS_events_anno/shortRead/230816_M71053_29_000000000-LBGW4/MI2077x04_S4_a1__CAG_EGFP_ITR_corrected__IS.csv', usecols = col)

df = pd.read_csv(infile, usecols = col)
# Prepare an empty list to hold the grouped data
grouped_data = []

# Function to extract the last two strings separated by ";"
def extract_chrom_pos(s):
    # Evaluate the string representation of the list to a list
    s_eval = ast.literal_eval(s)
    # Extract the last two strings for each element in the list
    new_list = [";".join(item.split(";")[-2:]) for item in s_eval]

    return list(set(new_list))

if df.empty:
    grouped_df = pd.DataFrame(columns=['sample', 'vector', 'chrom', 'IS_pos', 'nReads', 'unique_sonications', 'vector_clip_pos'])
    grouped_df.to_csv(outfile, index = False)
else: # process non-empty files
    # convert string representation of list to actual list
    df['vector_align_clip_pos'] = df['vector_align_clip_pos'].apply(ast.literal_eval)

    # Apply the function to each element of the column
    df['sonication_sites'] = df['sonication_sites'].apply(extract_chrom_pos)

    # Sort the dataframe by 'chrom' and 'IS_pos'
    df = df.sort_values(by=['chrom', 'IS_pos'])
    # print(df)

    # Variables to hold the cumulative sums, the min and max position for each group, and a list of original positions and original vector_align_clip_pos
    last_chrom = None
    min_pos = None
    max_pos = None
    sum_reads = 0
    # sum_sonications = 0
    original_IS = []
    original_clips = []
    original_sonications = []

    for index, row in df.iterrows():
        sample = row['sample']
        vector = row['vector']
        if row['chrom'] != last_chrom or row['IS_pos'] > max_pos + 10:
            # If the chromosome changed or the position gap is more than 10bp, save the previous group
            if last_chrom is not None:
                midpoint = (min_pos + max_pos) // 2  # Calculate the midpoint of the group
                grouped_data.append([sample, vector, last_chrom, midpoint, sum_reads, 
                    ";".join(map(str, original_IS)), original_clips, original_sonications])
            # Reset sums, min and max position, and original positions for the new group
            last_chrom = row['chrom']
            min_pos = row['IS_pos']
            max_pos = row['IS_pos']
            sum_reads = row['nReads']
            # sum_sonications = row['unique_sonications']
            original_IS = [row['IS_pos']]
            original_clips = row['vector_align_clip_pos']
            original_sonications = row['sonication_sites']
        else:
            # If within the same group, update the sums, adjust min and max if necessary, and add the current position
            sum_reads += row['nReads']
            # sum_sonications += row['unique_sonications']
            min_pos = min(min_pos, row['IS_pos'])
            max_pos = max(max_pos, row['IS_pos'])
            original_IS.append(row['IS_pos'])
            original_clips += row['vector_align_clip_pos']
            original_sonications += row['sonication_sites']

    # Append the last group after finishing the loop
    midpoint = (min_pos + max_pos) // 2  # Calculate the midpoint for the last group
    grouped_data.append([sample, vector, last_chrom, midpoint, sum_reads, ";".join(map(str, original_IS)), original_clips, original_sonications])

    # Convert the grouped data back to a DataFrame
    grouped_df = pd.DataFrame(grouped_data, columns=['sample', 'vector', 'chrom', 'IS_pos', 'nReads', 'original_IS_pos', 'original_clips', 'original_sonications'])

    # Function to find the most abundant value
    def find_most_abundant_value(value_list):
        if not value_list:
            return None
        count = Counter(value_list)
        max_count = max(count.values())
        most_common = [val for val, cnt in count.items() if cnt == max_count]
        return random.choice(most_common)  # Randomly choose among the most common if there's a tie

    grouped_df['unique_sonications'] = grouped_df['original_sonications'].apply(lambda x: len(list(set(x))))
    grouped_df['vector_clip_pos'] = grouped_df['original_clips'].apply(find_most_abundant_value)
    
    # print(grouped_df)

    grouped_df.drop(['original_IS_pos', 'original_clips', 'original_sonications'], axis = 1, inplace = True)
    grouped_df.to_csv(outfile, index = False)
