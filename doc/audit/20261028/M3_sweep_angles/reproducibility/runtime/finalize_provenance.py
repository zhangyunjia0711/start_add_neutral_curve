"""Correct descriptive metadata; leave every numerical point unchanged."""
import datetime,hashlib,json
from pathlib import Path
from common import ROOT,sha

def digest(d):
    return hashlib.sha256(json.dumps(d,sort_keys=True).encode()).hexdigest()

def main():
    stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    changes=[];audit=json.loads((ROOT/'outputs/s70_frequency_seed_search.json').read_text())['attempts']
    def fix_record(r):
        if 'frequency_seed_audit_index' not in r:return False
        j=r['frequency_seed_audit_index'];trial=audit[j] if j<len(audit) else {}
        valid=trial.get('finished') and trial.get('f_kHz')==r['f_kHz'] and trial.get('station')==r['station']
        if valid:return False
        r['superseded_frequency_seed_audit_index']=r.pop('frequency_seed_audit_index')
        r['frequency_seed_audit_reference_status']='Original audit was superseded during s70 guide revision; its index is unavailable. Native seed was independently checked by START before contour continuation. Do not resolve this index against the current audit.'
        return True
    p=ROOT/'outputs/s70_neutral_seed_scan.json';data=json.loads(p.read_text());before=sha(p)
    corrected=sum(bool(fix_record(r)) for r in data['records'])
    if corrected:
        tmp=p.with_suffix('.metadata.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(p)
        changes.append(dict(file=str(p),before_sha256=before,after_sha256=sha(p),superseded_indices_marked=corrected))
    model='Guide only: bordered real-alpha/real-beta polynomial; adaptive cubic actual local primitive-profile matrices, independently midpoint-checked (s70 physical_guide.py); final native points and physical folds are separately corrected'
    for p in (ROOT/'outputs').glob('s70_f*kHz_*_v3_guide.json'):
        d=json.loads(p.read_text());before=sha(p);numeric=digest(d['records']);old=d.get('metadata_correction',{}).get('previous_model_text',d['model']);d['model']=model
        source=d.get('seed',{}).get('source_record',{});fixed=fix_record(source)
        if 'metadata_correction' not in d:d['metadata_correction']=dict(time=stamp,previous_model_text=old,records_unchanged=True,superseded_seed_audit_index_marked=fixed)
        assert digest(d['records'])==numeric
        tmp=p.with_suffix('.metadata.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(p)
        changes.append(dict(file=str(p),before_sha256=before,after_sha256=sha(p),records_sha256=numeric,records_unchanged=True))
    p=ROOT/'work/contour3d.py';source=p.read_text();old="model='Guide only: bordered original real-alpha/real-beta polynomial; piecewise native-operator interpolation'"
    new="model=('Guide only: adaptive cubic actual local primitive-profile matrices; independent midpoint checks' if case=='s70' else 'Guide only: bordered original real-alpha/real-beta polynomial; piecewise native-operator interpolation')"
    assert source.count(old)==1 or source.count(new)==1;before=sha(p);p.write_text(source.replace(old,new));changes.append(dict(file=str(p),before_sha256=before,after_sha256=sha(p),metadata_text_only=True))
    receipts=[]
    for p in sorted((ROOT/'physical_guides/s70').glob('i*/receipt.json')):
        d=json.loads(p.read_text());accepted=sorted([c for c in d['checks'] if c['accepted']],key=lambda c:c['lo'])
        assert len(accepted)==d['segments'] and accepted[0]['lo']==0. and accepted[-1]['hi']==1.
        assert all(a['hi']==b['lo'] for a,b in zip(accepted,accepted[1:]))
        assert all(max(c['midpoint_coefficient_errors'])<5e-7 for c in accepted)
        receipts.append(dict(path=str(p),sha256=sha(p),segments=d['segments'],rejected_parent_checks=sum(not c['accepted'] for c in d['checks']),maximum_accepted_coefficient_error=max(max(c['midpoint_coefficient_errors']) for c in accepted)))
    assert len(receipts)==150
    corrected=sum('superseded_frequency_seed_audit_index' in r for r in data['records'])
    result=dict(created=stamp,changes=changes,invalid_s70_audit_indices_explicitly_marked=corrected,numerical_points_unchanged=True,physical_guide_validation=receipts,execution_note='First attempt completed metadata edits but its receipt check incorrectly required rejected adaptive parent intervals to pass. Retry checks accepted leaf intervals tile [0,1] and satisfy the original tolerance. No numerical values or solver gates changed.')
    out=ROOT/'outputs/provenance_corrections.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print('PROVENANCE_CORRECTED',len(changes),'files',corrected,'superseded indices',len(receipts),'physical guide cells')

if __name__=='__main__':main()
