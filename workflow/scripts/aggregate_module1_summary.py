"""Collects one JSON summary per sample (written by qc_preprocess_sample.py)
into a single module1_summary.csv, replicating the combined table the
original single-script run printed/wrote at the end of its `main()`."""
import json
import pandas as pd

rows = [json.load(open(p)) for p in snakemake.input.summaries]
df = pd.DataFrame(rows)
df.to_csv(snakemake.output[0], index=False)
print(df.to_string(index=False))
