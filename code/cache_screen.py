"""Exploratory phase-cache delta screen, followed by warm CW when inconclusive."""
from pathlib import Path
import numpy as np,csv,json,time
import study as s
P=Path(__file__).resolve().parent

def cache_trajectory(op,v):
    z=v.copy();states=[z.copy()];p=op.n//2
    for k in op.word:
        if k==0:z=np.r_[op.d[:p]*z[:p]+op.X@z[p:],op.d[p:]*z[p:]]
        else:z=np.r_[op.d[:p]*z[:p],op.d[p:]*z[p:]+op.Y@z[:p]]
        states.append(z.copy())
    ratios=z/v
    return np.array(states),float(min(ratios)),float(max(ratios))

def delta_prepare(old,new):
    dd=new.d-old.d;dx=new.X-old.X;dy=new.Y-old.Y
    return dd,[(r,dx[r].copy()) for r in np.flatnonzero(np.any(dx!=0,axis=1))],[(r,dy[r].copy()) for r in np.flatnonzero(np.any(dy!=0,axis=1))]

def delta_bounds_scalar(old,new,cache,delta):
    states,qlo,qhi=cache;dd,dx,dy=delta;p=old.n//2;loglo=np.log(qlo);loghi=np.log(qhi)
    for t,(ko,kn) in enumerate(zip(old.word,new.word)):
        v=states[t];z=states[t+1]
        if ko==kn:
            change=dd*v
            if kn==0:
                for r,row in dx:change[r]+=row@v[p:]
            else:
                for r,row in dy:change[p+r]+=row@v[:p]
        else:
            # Order mutation changes a whole phase direction; no sparse shortcut.
            change=dd*v
            if ko==0:change[:p]-=old.X@v[p:];change[p:]+=new.Y@v[:p]
            else:change[p:]-=old.Y@v[:p];change[:p]+=new.X@v[p:]
        ratios=1+change/z
        if min(ratios)<=0:raise ValueError('Loss of positivity / cancellation in delta ratio')
        loglo+=np.log(min(ratios));loghi+=np.log(max(ratios))
    lo=float(np.exp(loglo));hi=float(np.exp(loghi))
    verdict=-1 if hi<1-s.GUARD else 1 if lo>1+s.GUARD else 0
    return verdict,lo,hi

def delta_bounds(old,new,cache,delta):
    states,qlo,qhi=cache;dd,dx,dy=delta;p=old.n//2
    v=states[:-1];change=v*dd
    maskx=(old.word==0)&(new.word==0);masky=(old.word==1)&(new.word==1)
    for row,values in dx:change[maskx,row]+=v[maskx,p:]@values
    for row,values in dy:change[masky,p+row]+=v[masky,:p]@values
    xy=(old.word==0)&(new.word==1);yx=(old.word==1)&(new.word==0)
    if np.any(xy):
        change[xy,:p]-=v[xy,p:]@old.X.T;change[xy,p:]+=v[xy,:p]@new.Y.T
    if np.any(yx):
        change[yx,p:]-=v[yx,:p]@old.Y.T;change[yx,:p]+=v[yx,p:]@new.X.T
    ratios=1+change/states[1:]
    if np.min(ratios)<=0:raise ValueError('Loss of positivity / cancellation')
    lo=float(np.exp(np.log(qlo)+np.log(np.min(ratios,axis=1)).sum()))
    hi=float(np.exp(np.log(qhi)+np.log(np.max(ratios,axis=1)).sum()))
    verdict=-1 if hi<1-s.GUARD else 1 if lo>1+s.GUARD else 0
    return verdict,lo,hi

def main():
    configs=json.loads((P/'configs.json').read_text());oldrows=list(csv.DictReader((P/'baselines.csv').open()));mainrows=list(csv.DictReader((P/'cases.csv').open()));lookup={(int(r['n']),int(r['rep']),r['regime'],r['kind'],float(r['level'])):r for r in mainrows}
    rows=[];bench=[];cachecost=[];start=time.perf_counter()
    for cfg in configs:
        n=cfg['n'];rep=cfg['rep'];F=np.array(cfg['F']);G=np.array(cfg['G']);d=np.array(cfg['decay']);w=np.array(cfg['word'])
        for baseline in [r for r in oldrows if int(r['n'])==n and int(r['rep'])==rep]:
            regime=baseline['regime'];beta=float(baseline['beta']);old=s.prepare(F,G,d,beta,w);_,v=s.reference(old,True)
            t=time.perf_counter();cache=cache_trajectory(old,v);offline=time.perf_counter()-t;cachecost.append(dict(n=n,rep=rep,regime=regime,cache_s=offline,cache_payload_bytes=cache[0].nbytes))
            for kind in s.KINDS:
                for lev in s.LEVELS:
                    ff,gg,dd,ww=s.mutate(F,G,d,w,kind,lev);new=s.prepare(ff,gg,dd,beta,ww)
                    t=time.perf_counter();delta=delta_prepare(old,new);ds=time.perf_counter()-t
                    t=time.perf_counter();ver,lo,hi=delta_bounds(old,new,cache,delta);sec=time.perf_counter()-t
                    sv,sl,sh=delta_bounds_scalar(old,new,cache,delta)
                    assert sv==ver and np.allclose([sl,sh],[lo,hi],rtol=1e-12,atol=1e-12)
                    r=lookup[(n,rep,regime,kind,lev)];q=np.exp(float(r['lambda_new'])*new.period);truth=int(r['truth'])
                    assert lo<=q+1e-8 and q<=hi+1e-8,(n,rep,regime,kind,lo,q,hi)
                    assert not ver or ver==truth
                    rows.append(dict(n=n,rep=rep,regime=regime,kind=kind,level=lev,cache_verdict=ver,cache_lo=lo,cache_hi=hi,delta_prepare_s=ds,cache_eval_s=sec,hybrid_s=ds+sec+(float(r['warm_s']) if not ver else 0),warm_s=float(r['warm_s'])))
                    if rep<3 and regime=='near' and lev==.1 and kind in ['weight','joint','order']:
                        times={m:[] for m in ['cache','hybrid','warm','cold','dense']}
                        for repeat in range(7):
                            methods=list(times)
                            if repeat%2:methods.reverse()
                            for method in methods:
                                t=time.perf_counter()
                                if method in ['cache','hybrid']:
                                    dt=delta_prepare(old,new);vr,_,_=delta_bounds(old,new,cache,dt)
                                    if method=='hybrid' and not vr:s.screen(new,v)
                                elif method=='warm':s.screen(new,v)
                                elif method=='cold':s.screen(new,np.ones(n))
                                else:s.reference(new)
                                times[method].append(time.perf_counter()-t)
                        bench.append(dict(n=n,rep=rep,kind=kind,cache_resolved=int(ver!=0),**{m+'_median_s':float(np.median(t)) for m,t in times.items()}))
    s.dump_csv('cache_cases.csv',rows);s.dump_csv('cache_timings.csv',bench);s.dump_csv('cache_offline.csv',cachecost)
    (P/'cache_metadata.json').write_text(json.dumps(dict(design='Exploratory algorithm extension after main CW and duration studies; all main candidates retained',n_cases=len(rows),n_timing_cases=len(bench),elapsed_s=time.perf_counter()-start,scope='float64 comparison bounds with guard; no interval validation of this layer yet'),indent=2))
    print('Cache cases',len(rows),'elapsed',time.perf_counter()-start,flush=True)
    for reg in ['core','far','near','critical']:
        rr=[r for r in rows if r['regime']==reg];print(reg,'resolved',sum(r['cache_verdict']!=0 for r in rr)/len(rr),flush=True)
if __name__=='__main__':main()
