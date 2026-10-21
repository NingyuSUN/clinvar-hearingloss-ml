"""Descriptive only: repeat splits of the same cohort are not independent."""
import argparse
import json
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--statistics', required=True, type=Path)
a = p.parse_args()
s = pd.read_csv(a.statistics / 'seed_metrics.csv')
c = pd.read_csv(a.statistics / 'seed_contrasts.csv')
r = {'scope': 'Descriptive stability across ten reused-cohort splits; not ten independent data sets',
     'arms': {}, 'contrasts': {}}
for arm, data in s.groupby('arm'):
    r['arms'][arm] = {k: {'min': float(data[k].min()), 'max': float(data[k].max())}
                      for k in ['auc', 'FN', 'FP', 'recall']}
for contrast, data in c.groupby('contrast'):
    r['contrasts'][contrast] = {k: {'min': float(data[k].min()), 'max': float(data[k].max()),
        'positive_seeds': int((data[k] > 0).sum()), 'negative_seeds': int((data[k] < 0).sum()),
        'zero_seeds': int((data[k] == 0).sum())} for k in ['auc', 'recall', 'fpr']}
(a.statistics / 'seed_stability_summary.json').write_text(json.dumps(r, indent=2) + '\n')
