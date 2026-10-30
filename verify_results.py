"""Compare immutable recorded results with an isolated rerun (stdlib only)."""
import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def close(a, b):
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-12)
    except (ValueError, TypeError):
        return a == b

def main():
    recorded = Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'results/recorded'
    rerun = Path(sys.argv[2]) if len(sys.argv)>2 else ROOT/'results/verified-rerun'
    specs = {
        'cases.csv': ['n','rep','regime','kind','level','beta','lambda_old','lambda_new','truth','stable_to_unstable'] + [m+s for m in ['warm','cold'] for s in ['_verdict','_steps','_lo','_hi','_resolved_1','_resolved_2','_resolved_4','_resolved_8','_resolved_16']],
        'duration_sensitivity.csv': ['n','rep','h','period','target','lambda_old','calibration_error','kind','level','lambda_new'] + [m+s for m in ['warm','cold'] for s in ['_verdict','_steps','_resolved_1','_resolved_2','_resolved_4','_resolved_8','_resolved_16']],
        'heldout_cases.csv': ['n','rep','regime','kind','level','lambda_old','lambda_new','truth','route','adaptive_verdict','adaptive_steps','warm_verdict','warm_steps','cold_verdict','cold_steps'],
        'cache_cases.csv': ['n','rep','regime','kind','level','cache_verdict','cache_lo','cache_hi'],
        'resolution_summary.csv': None,
        'duration_summary.csv': None,
        'cache_summary.csv': None,
    }
    checks=[]
    for name, columns in specs.items():
        a,b=read(recorded/name),read(rerun/name)
        columns=columns or list(a[0])
        differences=[]
        if len(a)!=len(b): differences.append({'row_counts':[len(a),len(b)]})
        for i,(x,y) in enumerate(zip(a,b)):
            for k in columns:
                if not close(x[k],y[k]): differences.append({'row':i,'field':k,'recorded':x[k],'rerun':y[k]})
        checks.append({'file':name,'rows':len(b),'fields':columns,'differences':differences,'passed':not differences})
    for name in ['configs.json','heldout_configs.json','edge_checks.json','interval_checks.json']:
        same=json.loads((recorded/name).read_text())==json.loads((rerun/name).read_text())
        checks.append({'file':name,'exact_json_match':same,'passed':same})
    verdict_errors={}
    for name,methods in [('cases.csv',['warm','cold']),('heldout_cases.csv',['warm','cold','adaptive'])]:
        rs=read(rerun/name)
        verdict_errors[name]={m:sum(int(r[m+'_verdict'])!=0 and r[m+'_verdict']!=r['truth'] for r in rs) for m in methods}
    duration=read(rerun/'duration_sensitivity.csv')
    verdict_errors['duration_sensitivity.csv']={m:sum(int(r[m+'_verdict'])!=0 and int(r[m+'_verdict'])!=(-1 if float(r['lambda_new'])<0 else 1) for r in duration) for m in ['warm','cold']}
    timings=[]
    for name,group,fields in [('timing_summary.csv','n',['dense_ms','warm_ms','cold_ms','arpack_ms']),('heldout_summary.csv','regime',['adaptive_us','warm_us','cold_us'])]:
        for a,b in zip(read(recorded/name),read(rerun/name)):
            for field in fields:
                timings.append({'table':name,'group':a[group],'metric':field,'recorded':float(a[field]),'rerun':float(b[field])})
    manifest=json.loads((recorded/'manifest.json').read_text())
    original_mismatches=[name for name,digest in manifest.items() if hashlib.sha256((recorded/name).read_bytes()).hexdigest()!=digest]
    result={'all_comparisons_passed':all(c['passed'] for c in checks), 'comparison_tolerance':{'relative':1e-10,'absolute':1e-12},'checks':checks,'decisive_verdict_errors':verdict_errors,'timings':timings,'original_manifest_mismatches':original_mismatches,'timing_note':'Timings are measured, not expected to reproduce exactly. Recorded results remain unchanged.'}
    target=ROOT/'results/comparison';target.mkdir(exist_ok=True)
    (target/'comparison.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'all_comparisons_passed':result['all_comparisons_passed'],'checks':[(x['file'],x['passed']) for x in checks],'verdict_errors':verdict_errors,'original_manifest_mismatches':original_mismatches},indent=2))
    if not result['all_comparisons_passed'] or any(v for d in verdict_errors.values() for v in d.values()):
        raise SystemExit(1)

if __name__=='__main__':main()
