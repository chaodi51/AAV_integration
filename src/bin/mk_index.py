import sys
import pandas as pd
from collections import defaultdict

template, datasets, output = sys.argv[1:4]
rmd_ls = sys.argv[4:]

rmd_dict = defaultdict(str)
for rmd in rmd_ls:
    s = rmd.split("/")[-1].split(".Rmd")[0]
    run = s.split('--')[1]
    rmd_dict[run] = rmd


data = pd.read_excel(datasets, engine = 'openpyxl')

with open(template) as f, open(output, "w") as fout:
    for line in f:
        print(line.strip(), file = fout)

    fout.write('| ' + ' | '.join(data.columns) + ' |\n')
    fout.write('| ' + ' | '.join(['---'] * len(data.columns)) + ' |\n')

    for index, row in data.iterrows():
        # Convert each cell to a string, handling missing values
        row_data = [str(cell) if pd.notnull(cell) else '' for cell in row]
        if isinstance(row['date sequencing'], pd.Timestamp):
            row['date sequencing'] = row['date sequencing'].strftime('%Y-%m-%d')
        
        runID = row['Illumina ID']
        if runID in rmd_dict:
            # Replace 'Run ID' with markdown link
            run_id_html = rmd_dict[runID].replace('Rmd', 'html').split('/')[-1]
            row_data[0] = f"[{row['Run ID']}]({run_id_html})"

        # Write the row to the .rmd file
        fout.write('| ' + ' | '.join(row_data) + ' |\n')