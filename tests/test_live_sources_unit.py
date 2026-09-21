"""Offline regression contracts for live inputs; no service or API key needed."""
import copy
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'predict'))
import live_sources as ls
import predict_novel as novel
from prediction_core import Bundle, MODELS, extract_features


def response(monkeypatch, payload, method='post', status=200):
    monkeypatch.setattr(ls.requests, method, lambda *a, **kw: SimpleNamespace(status_code=status, json=lambda: payload))


def vep_record():
    return {'assembly_name': 'GRCh38', 'seq_region_name': '16', 'start': 2496587,
            'end': 2496587, 'strand': 1, 'allele_string': 'G/C',
            'transcript_consequences': [{'transcript_id': 'ENST1', 'gene_symbol': 'TBC1D24',
                'consequence_terms': ['missense_variant'], 'mane_select': 'NM1',
                'biotype': 'protein_coding', 'protein_start': 10, 'domains': []}]}


@pytest.mark.parametrize('payload', [{}, {'errors': [{'message': 'unavailable'}]}, {'data': None},
    {'data': {}}, {'data': {'variant': None}, 'errors': [{'message': 'failed'}]}, {'data': {'variant': []}}])
def test_api_failure_is_not_absence(monkeypatch, payload):
    response(monkeypatch, payload)
    with pytest.raises(ls.AnnotationError):
        ls.fetch_gnomad_af('16', 2496587, 'G', 'C')


def test_verified_absence_matches_all_frozen_absent_rows(monkeypatch):
    response(monkeypatch, {'data': {'variant': None}})
    features = ls.fetch_gnomad_af('16', 2496587, 'G', 'C')
    cached = [json.loads(line) for line in (ROOT/'predict/bundle/annotation_cache.jsonl').read_text().splitlines()]
    absent = [a['features']['p4_freq_log10']['value'] for a in cached
              if a['features']['p2_not_observed_in_gnomad_r4']['value'] == 1]
    assert len(absent) == 2085
    assert all(v == features['p4_freq_log10']['value'] for v in absent)


@pytest.mark.parametrize('part', [{'af': .1, 'an': 100}, {'af': .1, 'an': 100, 'filters': None},
    {'af': float('nan'), 'an': 100, 'filters': []}, {'af': 1.1, 'an': 100, 'filters': []},
    {'af': True, 'an': 100, 'filters': []}, {'af': .1, 'an': 0, 'filters': []},
    {'af': .1, 'an': True, 'filters': []}])
def test_bad_frequency_rejected(monkeypatch, part):
    response(monkeypatch, {'data': {'variant': {'exome': part, 'genome': None}}})
    with pytest.raises(ls.AnnotationError):
        ls.fetch_gnomad_af('16', 2496587, 'G', 'C')


def test_true_zero_and_filtered_have_distinct_semantics(monkeypatch):
    response(monkeypatch, {'data': {'variant': {'exome': {'af': 0, 'an': 3000, 'filters': []}, 'genome': None}}})
    zero = ls.fetch_gnomad_af('16', 2496587, 'G', 'C')
    assert zero['p4_freq_log10']['value'] == -8
    assert zero['p2_not_observed_in_gnomad_r4']['value'] == 0
    response(monkeypatch, {'data': {'variant': {'exome': {'filters': ['AC0']}, 'genome': None}}})
    filtered = ls.fetch_gnomad_af('16', 2496587, 'G', 'C')
    values, errors = extract_features({'origin': novel.ORIGIN, 'features': filtered}, ['p4_freq_log10'], {'p4_freq_log10'})
    assert values is None and errors == ['p4_freq_log10:filtered_or_unavailable']


@pytest.mark.parametrize('error', [requests.Timeout(), requests.ConnectionError(), ValueError('bad JSON')])
def test_transport_and_decode_errors_normalized(monkeypatch, error):
    def fail(*a, **kw):
        raise error
    monkeypatch.setattr(ls.requests, 'post', fail)
    with pytest.raises(ls.AnnotationError):
        ls.fetch_gnomad_af('16', 2496587, 'G', 'C')


@pytest.mark.parametrize('field,value', [('assembly_name', 'GRCh37'), ('seq_region_name', '1'),
    ('start', 2496588), ('end', 2496588), ('strand', -1), ('allele_string', 'A/C')])
def test_vep_checks_identity_even_after_http_200(monkeypatch, field, value):
    record = vep_record(); record[field] = value
    response(monkeypatch, [record], 'get')
    with pytest.raises(ls.AnnotationError, match='identity_mismatch'):
        ls.fetch_vep('16', 2496587, 'G', 'C')


@pytest.mark.parametrize('transcripts', [[], None, [None], [{'transcript_id': 'T', 'consequence_terms': 'missense_variant'}]])
def test_vep_malformed_transcripts_fail_as_annotation(monkeypatch, transcripts):
    record = vep_record(); record['transcript_consequences'] = transcripts
    response(monkeypatch, [record], 'get')
    with pytest.raises(ls.AnnotationError):
        ls.fetch_vep('16', 2496587, 'G', 'C')


def test_vep_valid_response(monkeypatch):
    response(monkeypatch, [vep_record()], 'get')
    assert ls.fetch_vep('16', 2496587, 'G', 'C')['variant_type'] == 'missense'


def test_phyloP_never_fills_training_feature(monkeypatch):
    response(monkeypatch, {'phyloP100way': [{'value': 4.2}]}, 'get')
    features = ls.fetch_conservation('16', 2496587)
    assert features['ucsc_phyloP100way_reference']['value'] == 4.2
    assert features['ensembl_conservation'] == {'value': None, 'status': 'unverified_source'}
    assert extract_features({'features': features, 'origin': novel.ORIGIN}, ['ensembl_conservation'], {'ensembl_conservation'})[0] is None


@pytest.mark.parametrize('gene', [{}, {'gnomad_constraint': []}, {'gnomad_constraint': {'oe_lof_upper': -1, 'oe_mis_upper': 1}},
    {'gnomad_constraint': {'oe_lof_upper': float('inf'), 'oe_mis_upper': 1}}])
def test_constraint_invalid_structure_or_values(monkeypatch, gene):
    response(monkeypatch, {'data': {'gene': gene}})
    with pytest.raises(ls.AnnotationError):
        ls.fetch_gnomad_constraint('TBC1D24')


def test_optional_gpn_dependency_missing_is_row_failure(monkeypatch):
    monkeypatch.setitem(sys.modules, 'polars', None)
    result = ls.fetch_gpn_star('16', 2496587, 'G', 'C')
    assert len(result) == 3 and all(v['value'] is None for v in result.values())


@pytest.fixture(scope='module')
def bundle():
    return Bundle(ROOT/'predict/bundle')


def annotation(bundle):
    ann = copy.deepcopy(bundle.cache['16-2496587-G-C'])
    ann['error'] = None
    ann['features']['ensembl_conservation'] = {'value': None, 'status': 'unverified_source'}
    return ann


def test_actual_models_only_score_compatible_gpn(bundle):
    scores = novel.score_annotation(annotation(bundle), bundle)
    assert scores['GPN3']['status'] == 'scored_research'
    for name in set(MODELS) - {'GPN3'}:
        assert scores[name]['score'] is None
        assert 'ensembl_conservation:unverified_source' in scores[name]['reasons']


def test_unsupported_type_cannot_reach_thresholds(bundle):
    ann = annotation(bundle); ann['variant_type'] = 'truncating'
    assert all(m['status'] == 'unsupported_variant_type' for m in novel.score_annotation(ann, bundle).values())


def test_model_that_saw_group_fails_integrity(bundle, monkeypatch):
    model = bundle.model
    def bad(*args):
        booster, meta, thresholds = model(*args)
        return booster, {**meta, 'excluded_groups': []}, thresholds
    monkeypatch.setattr(bundle, 'model', bad)
    with pytest.raises(ValueError, match='saw held-out'):
        novel.score_annotation(annotation(bundle), bundle)


def test_cli_retains_failed_rows_and_honest_membership(tmp_path, bundle, monkeypatch):
    input_path = tmp_path/'input.csv'; out = tmp_path/'out'
    input_path.write_text('variant_id,assembly,chrom,pos,ref,alt\nknown,GRCh38,16,2496587,G,C\n'
        'failed,GRCh38,1,10,A,C\nunknown,GRCh38,1,11,A,C\nbad,GRCh37,1,12,A,C\n', encoding='utf-8')
    def annotate(row, key):
        ann = annotation(bundle)
        ann['variant_key'] = f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}"
        if row['variant_id'] == 'failed':
            ann['error'] = 'vep_failed:vep_transport_error'
        return ann
    monkeypatch.setattr(novel, 'annotate_variant', annotate)
    monkeypatch.setattr(novel, 'Bundle', lambda p: bundle)
    monkeypatch.setattr(sys, 'argv', ['predict_novel', '--input', str(input_path), '--bundle', str(ROOT/'predict/bundle'), '--output', str(out)])
    monkeypatch.delenv('ALPHAGENOME_API_KEY', raising=False)
    novel.main()
    rows = [json.loads(line) for line in (out/'predictions.jsonl').read_text(encoding='utf-8').splitlines()]
    assert [r['input_index'] for r in rows] == [0, 1, 2, 3]
    assert rows[0]['in_development_7125'] is True and rows[0]['in_common_4050'] is True
    assert rows[1]['annotation_status'] == 'vep_failed:vep_transport_error'
    assert rows[2]['in_development_7125'] is None and rows[2]['in_common_4050'] is None
    assert rows[3]['annotation_status'] == 'unsupported_assembly'


def test_vep_domain_database_must_be_scalar(monkeypatch):
    record = vep_record()
    record['transcript_consequences'][0]['domains'] = [{'db': ['Pfam']}]
    response(monkeypatch, [record], 'get')
    with pytest.raises(ls.AnnotationError, match='invalid_transcript'):
        ls.fetch_vep('16', 2496587, 'G', 'C')


def test_gpn_schema_error_is_contained_per_feature(monkeypatch):
    class ColumnError(Exception):
        pass
    class Frame:
        height = 1
        def filter(self, expr): return self
        def collect(self): return self
        def is_empty(self): return False
        def __getitem__(self, key): raise ColumnError('missing llr_calibrated')
    class Expr:
        def __eq__(self, value): return self
        def __and__(self, other): return self
    monkeypatch.setitem(sys.modules, 'polars', SimpleNamespace(scan_parquet=lambda p: Frame(), col=lambda c: Expr()))
    result = ls.fetch_gpn_star('16', 2496587, 'G', 'C')
    assert all(r['value'] is None and 'llr_calibrated' in r['error'] for r in result.values())
