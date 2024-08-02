import sys
import pandas as pd

pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.max_colwidth', None) 

run, mapfile, infile, outfile = sys.argv[1:]

runid = run.split('/')[1] # shortRead/240125_A01959_0044_BHVG7FDRX3
master_map = pd.read_csv(mapfile, delimiter = '\t')
run_samples = master_map.loc[master_map['run'] == runid, ['run', 'Sample', 'biosample']]

df = pd.read_csv(infile)

# add biosample names to the IS table
df = pd.merge(df, run_samples, left_on='sample', right_on='Sample', how='left')

special_runs = ['240405_M05240_0149_000000000-L9FLR', '240328_VH00163_94_AAFJHJGM5','240405_M05240_0149_000000000-L9FLR_mouse', '240328_VH00163_94_AAFJHJGM5_mouse']
if runid in special_runs:
    df.drop(['run', 'Sample'], axis = 1, inplace = True)
    df.to_csv(outfile, index = False)
    sys.exit()
else:
    df.drop(['run', 'Sample'], axis = 1, inplace = True)

samples = ['S' + str(i) for i in range(5, 13)]

# function to compare collisions (close IS in unrelated samples)
# rules from Ali:
#		a) the integration site with the higher sonication number is kept and the others are discarded. 
#		b) If the sonication number is the same the integration site with the higher read number is kept.
#		c) If both read count and sonication number are the same both are kept
#		d) Correction of point c (011524) if both read count and sonication site are the same they are all removed.

def compare_info(row):
    parts = row['original_info'].split('--')
    if len(parts) == 1:
        return row['original_info']  # Only one part, no comparison needed
    
    # Parse each part
    parsed_parts = [part.split(';') for part in parts]
    # Convert relevant fields to integers for comparison
    parsed_parts = [(part, int(part[4]), int(part[3])) for part in parsed_parts]  # (original_string, sonication_sites, nRead)
    
    # Sort based on sonication_sites, then nRead, both in descending order
    sorted_parts = sorted(parsed_parts, key=lambda x: (x[1], x[2]), reverse=True)
    
    # Keep the part with the larger sonication_sites, or nRead in case of a tie
    if sorted_parts[0][1] > sorted_parts[1][1]:
        return ';'.join(sorted_parts[0][0])  # Return the string of the kept part
    elif sorted_parts[0][2] > sorted_parts[1][2]:
        return ';'.join(sorted_parts[0][0])
    else:
        return "drop_IS"  # Indicate a tie on both fields, row should be removed

# Function to find clusters of IS_pos within 10 bp
def collision_rm(df_comp, current_sample, comparison_sample):
    grouped_data = []
    last_chrom = None
    max_pos = None
    original_IS = []

    for index, row in df_comp.iterrows():
        sample = row['sample']
        vector = row['vector']
        original_IS = row.tolist()

        if row['chrom'] != last_chrom or int(row['IS_pos']) > max_pos + 10:
            # If the chromosome changed or the position gap is more than 10bp, save the previous group
            if last_chrom is not None:
                grouped_data.append([sample, vector, last_chrom, max_pos,
                    "--".join(map(str, original_info))])
            # Reset info for the new group
            last_chrom = row['chrom']
            max_pos = int(row['IS_pos'])
            original_info = [";".join(map(str, original_IS))]
        else:
            # If within the same group add the current IS info
            original_info += [";".join(map(str, original_IS))]
            max_pos = max(max_pos, int(row['IS_pos']))

    grouped_data.append([sample, vector, last_chrom, max_pos, "--".join(map(str, original_info))])
    grouped_df = pd.DataFrame(grouped_data, columns=['sample', 'vector', 'chrom', 'IS_pos', 'original_info'])
    # print(grouped_df)
    grouped_df = grouped_df.apply(lambda x: compare_info(x), axis=1).to_frame('updated_info')
    # print(grouped_df)
    grouped_df = grouped_df[grouped_df['updated_info'] != 'drop_IS']

    # split the updated_info to get the same format as original df
    df_split = grouped_df['updated_info'].str.split(';', expand=True)
    # print(df_split)
    df_split.columns = list(df_comp.columns)
    df_split[['IS_pos','nReads','unique_sonications','vector_clip_pos']] = df_split[['IS_pos','nReads','unique_sonications','vector_clip_pos']].astype('int')

    return df_split

# for test
# current_sample = 'S5'
# comparison_sample = 'S7'
# df_comp = df[df['biosample'].isin([current_sample, comparison_sample])]
# df_comp = df_comp.sort_values(by=['chrom', 'IS_pos'])
# print(df_comp.shape)
# print(df_comp)
# df = collision_rm(df_comp, current_sample, comparison_sample)

sample_pair = {
    'S5':'S9', 'S9':'S5',
    'S6':'S10', 'S10':'S6',
    'S7':'S11', 'S11':'S7',
    'S8':'S12', 'S12':'S8'
}

for index, current_sample in enumerate(samples):
    # Create a new list that includes all samples after the current one
    following_samples = samples[index + 1:]
    
    # Exclude the sample's paired partner
    paired_partner = sample_pair[current_sample]
    filtered_samples = [s for s in following_samples if s != paired_partner]
    
    # compare current_sample vs each in the filtered_samples
    for comparison_sample in filtered_samples:
        print(current_sample, "vs", comparison_sample)
        df_comp = df[df['biosample'].isin([current_sample, comparison_sample])]
        df_rest = df[~df['biosample'].isin([current_sample, comparison_sample])]
        df_comp = df_comp.sort_values(by=['chrom', 'IS_pos'])
        # print(df_comp.shape)
        # print(df.dtypes)
        # print(df_comp)
        df_f = collision_rm(df_comp, current_sample, comparison_sample)
        df = pd.concat([df_rest, df_f], ignore_index=True) 

df.to_csv(outfile, index = False)

