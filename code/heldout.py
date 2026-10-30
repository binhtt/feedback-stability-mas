"""Independent-seed assessment of a fixed adaptive gate after exploratory cache results.
Gate chosen before this run: |old lambda| >= 0.10, unchanged phase directions.
No candidate/reference spectral result enters the adaptive decision.
"""
import numpy as np,time,json
from pathlib import Path
import study as s
import cache_screen as c
P=Path(__file__).resolve().parent
SEED=202610300
THRESHOLD=.10

def adaptive(old,new,cache,v,old_lambda,kind):
    if abs(old_lambda)>=THRESHOLD and kind!='order':
        delta=c.delta_prepare(old,new);ver,_,_=c.delta_bounds(old,new,cache,delta)
        if ver:return ver,'cache',0
    ver,k,_,_,_=s.screen(new,v)
    return ver,'warm',k

def main():
    rng=np.random.default_rng(SEED);rows=[];configs=[];start=time.perf_counter()
    for n in [12,48,192]:
        for rep in range(10):
            F,G,d,w=s.base_model(rng,n);corebeta=float(rng.uniform(1,6))
            configs.append(dict(n=n,rep=rep,F=F.tolist(),G=G.tolist(),decay=d.tolist(),word=w.tolist(),core_beta=corebeta))
            for regime,target in [('core',None),('far',-.2),('near',-.02),('critical',-.002)]:
                beta=corebeta if target is None else s.calibrate(F,G,d,w,target)
                old=s.prepare(F,G,d,beta,w);t=time.perf_counter();qold,v=s.reference(old,True);cache=c.cache_trajectory(old,v);offline=time.perf_counter()-t;lamold=np.log(qold)/old.period
                for kind in s.KINDS:
                    for level in s.LEVELS:
                        ff,gg,dd,ww=s.mutate(F,G,d,w,kind,level);new=s.prepare(ff,gg,dd,beta,ww)
                        t=time.perf_counter();q=s.reference(new);dense=time.perf_counter()-t;truth=-1 if q<1-s.GUARD else 1 if q>1+s.GUARD else 0
                        times={k:[] for k in ['adaptive','warm','cold']};out={}
                        for repeat in range(3):
                            methods=list(times)
                            if repeat%2:methods.reverse()
                            for method in methods:
                                t=time.perf_counter()
                                if method=='adaptive':ver,route,steps=adaptive(old,new,cache,v,lamold,kind);out['route']=route
                                else:ver,steps,_,_,_=s.screen(new,v if method=='warm' else np.ones(n))
                                times[method].append(time.perf_counter()-t);out[method+'_verdict']=ver;out[method+'_steps']=steps
                                assert not ver or ver==truth,(n,rep,regime,kind,level,method,ver,truth)
                        rows.append(dict(n=n,rep=rep,regime=regime,kind=kind,level=level,lambda_old=lamold,lambda_new=np.log(q)/new.period,truth=truth,offline_s=offline,dense_s=dense,**out,**{k+'_median_s':float(np.median(t)) for k,t in times.items()}))
        s.dump_csv('heldout_cases.csv',rows);print(n,len(rows),'elapsed',time.perf_counter()-start,flush=True)
    (P/'heldout_configs.json').write_text(json.dumps(configs))
    (P/'heldout_metadata.json').write_text(json.dumps(dict(seed=SEED,n_cases=len(rows),threshold=THRESHOLD,sizes=[12,48,192],reps=10,timing_repeats=3,gate='abs(previous lambda)>=0.10 and update kind != order; cache then warm fallback',design='Fixed after exploratory cache results, before independent-seed assessment; no candidate spectral result enters gate',elapsed_s=time.perf_counter()-start),indent=2))
if __name__=='__main__':main()
