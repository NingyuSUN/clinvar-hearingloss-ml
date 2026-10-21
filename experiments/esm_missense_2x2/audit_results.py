"""Independent server-side recomputation with sklearn and the frozen base AUC code."""
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--predictions', type=Path, required=True)
    ap.add_argument('--statistics', type=Path, required=True)
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    sys.path.insert(0, str(args.repo_root / 'src'))
    from hlpath.compare import paired_group_bootstrap_auc

    result = json.loads((args.statistics / 'summary.json').read_text())
    assert digest(args.predictions) == result['input_sha256']
    frame = pd.read_csv(args.predictions, dtype={'VariationID': str, 'group': str},
                        float_precision='round_trip')
    assert len(frame) == 4 * 1328 * 10
    assert not frame.duplicated(['arm', 'VariationID', 'seed']).any()
    assert (frame.decision == (frame.score >= frame.threshold).astype(int)).all()
    assert frame.groupby('VariationID').y.nunique().max() == 1
    assert frame.groupby('VariationID')['group'].nunique().max() == 1
    aggregate = {}
    for arm in ['A', 'B', 'C', 'D']:
        rows = frame[frame.arm == arm]
        assert set(rows.seed) == set(range(101, 111))
        assert (rows.groupby('VariationID').size() == 10).all()
        v = rows.groupby('VariationID', sort=True).agg(y=('y', 'first'),
            group=('group', 'first'), score=('score', 'mean'), votes=('decision', 'sum'))
        v['pred'] = (v.votes >= 5).astype(int)
        assert len(v) == 1328 and v.y.sum() == 925 and v.group.nunique() == 135
        tn, fp, fn, tp = confusion_matrix(v.y, v.pred, labels=[0, 1]).ravel()
        expected = result['summary'][arm]
        for k, value in {'TP': tp, 'FN': fn, 'FP': fp, 'TN': tn}.items():
            assert int(value) == expected[k], (arm, k, value, expected[k])
        assert abs(roc_auc_score(v.y, v.score) - expected['auc']) < 1e-12
        aggregate[arm] = v
    c, d = aggregate['C'], aggregate['D']
    pd.testing.assert_index_equal(c.index, d.index)
    assert (c.y == d.y).all() and (c.group == d.group).all()
    official = paired_group_bootstrap_auc(c.y.to_numpy(), c.score.to_numpy(),
        d.score.to_numpy(), groups=c.group.to_numpy(), n_boot=2000, seed=0)
    primary = {x['metric']: x for x in result['contrasts'] if x['contrast'] == 'D-C'}
    delta = official['delta_b_minus_a']
    for a, b in [('estimate', 'difference'), ('lo', 'ci95_low'), ('hi', 'ci95_high')]:
        assert abs(delta[a] - primary['auc'][b]) < 1e-12, (a, delta[a], primary['auc'][b])

    # Independently expand sampled groups into repeated rows for decision deltas.
    groups = sorted(c.group.unique())
    row_indices = [np.flatnonzero(c.group.to_numpy() == group) for group in groups]
    y, cp, dp = c.y.to_numpy(), c.pred.to_numpy(), d.pred.to_numpy()
    rng = np.random.default_rng(0)
    draws = {'recall': [], 'fpr': [], 'specificity': []}
    for _ in range(2000):
        ix = np.concatenate([row_indices[g] for g in rng.integers(0, len(groups), len(groups))])
        yy, cc, dd = y[ix], cp[ix], dp[ix]
        pos, neg = yy == 1, yy == 0
        if not pos.any() or not neg.any():
            continue
        recall = dd[pos].mean() - cc[pos].mean()
        fpr = dd[neg].mean() - cc[neg].mean()
        draws['recall'].append(recall)
        draws['fpr'].append(fpr)
        draws['specificity'].append(-fpr)
    for metric, values in draws.items():
        assert len(values) >= 1800
        lo, hi = np.quantile(values, [.025, .975])
        assert abs(lo - primary[metric]['ci95_low']) < 1e-12
        assert abs(hi - primary[metric]['ci95_high']) < 1e-12
    changes = {'rescued_pathogenic': int(((y == 1) & (cp == 0) & (dp == 1)).sum()),
               'new_pathogenic_misses': int(((y == 1) & (cp == 1) & (dp == 0)).sum()),
               'removed_benign_false_positives': int(((y == 0) & (cp == 1) & (dp == 0)).sum()),
               'new_benign_false_positives': int(((y == 0) & (cp == 0) & (dp == 1)).sum())}
    dc = result['C_D_discordance']
    assert changes['rescued_pathogenic'] == dc['pathogenic']['C0_D1']
    assert changes['new_pathogenic_misses'] == dc['pathogenic']['C1_D0']
    assert changes['removed_benign_false_positives'] == dc['benign']['C1_D0']
    assert changes['new_benign_false_positives'] == dc['benign']['C0_D1']
    audit = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'passed': True,
        'predictions_sha256': digest(args.predictions), 'script_sha256': digest(__file__),
        'statistics_sha256': digest(args.statistics / 'summary.json'),
        'checks': ['53120 rows, matched identities and 10 predictions per variant/arm',
                   'full-precision threshold decisions and >=5 votes',
                   'four confusion matrices and sklearn AUCs',
                   'primary paired AUC CI against unchanged base implementation',
                   'primary decision CIs by explicit group-row replication',
                   'paired error rescues and new errors'],
        'primary_comparison': primary, 'C_D_changes': changes,
        'limitations': 'Verifies arithmetic and stored prediction integrity; training isolation requires separate split and fit audit.'}
    args.out.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'passed': True, 'C_D_changes': changes}))


if __name__ == '__main__':
    main()
