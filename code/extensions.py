"""Post-main-study sensitivity: explicitly exploratory, no replacement of main results."""
import study as s
import numpy as np,json,time
from scipy.linalg import expm
from pathlib import Path
import mpmath as mp
P=Path(__file__).resolve().parent

def interval_action(F,G,d,beta,word,v,h):
    # Inputs converted from exact binary64 values: point intervals enclose those values.
    iv=mp.iv;iv.dps=60;n=len(d);p=n//2
    ds=[iv.mpf(float(x)) for x in d];a=[iv.exp(-x*h) for x in ds];b=iv.mpf(float(beta))
    X=[[None]*p for _ in range(p)];Y=[[None]*p for _ in range(p)]
    def c(i,j,w):
        # Integral response, continuous extension if losses coincide.
        q=(h*a[i] if d[i]==d[j] else (a[j]-a[i])/(ds[i]-ds[j]))
        return b*iv.mpf(float(w))*q*a[i]
    for i in range(p):
        for j in range(p):X[i][j]=c(i,p+j,F[i,j]);Y[i][j]=c(p+i,j,G[i,j])
    z=[iv.mpf(float(x)) for x in v];di=[x*x for x in a]
    for k in word:
        if k==0:z=[di[i]*z[i]+sum(X[i][j]*z[p+j] for j in range(p)) for i in range(p)]+[di[p+i]*z[p+i] for i in range(p)]
        else:z=[di[i]*z[i] for i in range(p)]+[di[p+i]*z[p+i]+sum(Y[i][j]*z[j] for j in range(p)) for i in range(p)]
    rat=[z[i]/iv.mpf(float(v[i])) for i in range(n)]
    # Convert outward interval endpoints to strings for an auditable enclosure.
    lo=min(r.a for r in rat);hi=max(r.b for r in rat)
    ver=-1 if hi<1 else 1 if lo>1 else 0
    return ver,str(lo),str(hi)

def main():
    cfgs=json.loads((P/'configs.json').read_text());rows=[];checks=[];start=time.perf_counter()
    # Validate analytic phase formula against general dense expm, independent implementation.
    s.H=.125
    for cfg in [c for c in cfgs if c['n'] in [6,12,24] and c['rep']<3]:
        n=cfg['n'];p=n//2;F=np.array(cfg['F']);G=np.array(cfg['G']);d=np.array(cfg['decay']);w=np.array(cfg['word']);beta=3.8
        op=s.prepare(F,G,d,beta,w);M=np.eye(n);idle=expm(-np.diag(d)*s.H)
        for k in w:
            B=-np.diag(d)
            if k==0:B[:p,p:]+=beta*F
            else:B[p:,:p]+=beta*G
            M=idle@expm(B*s.H)@M
        err=float(np.max(abs(M-op.dense())));assert err<1e-10
        checks.append(dict(test='analytic_phase_vs_expm',n=n,rep=cfg['rep'],error=err))
        rng=np.random.default_rng(721+n+cfg['rep']);v=rng.uniform(.2,2,n)
        err=float(np.max(abs(op.apply(v)-M@v)));assert err<1e-10
        checks.append(dict(test='structured_action_vs_dense_expm',n=n,rep=cfg['rep'],error=err))
    for h in [.005,.025,.125,.5]:
        s.H=h
        for cfg in [c for c in cfgs if c['n'] in [24,96] and c['rep']<5]:
            F=np.array(cfg['F']);G=np.array(cfg['G']);d=np.array(cfg['decay']);w=np.array(cfg['word']);n=cfg['n']
            for target in [-.2,-.02,-.002]:
                beta=s.brentq(lambda b: np.log(s.reference(s.prepare(F,G,d,b,w)))/(2*h*len(w))-target,0,32,xtol=1e-11);old=s.prepare(F,G,d,beta,w);qold,v=s.reference(old,True)
                # Dense root reference avoids a slow-mixing calibration bias.
                actual=float(np.log(qold)/old.period)
                assert abs(actual-target)<1e-9
                for kind in s.KINDS:
                    for lev in s.LEVELS:
                        ff,gg,dd,ww=s.mutate(F,G,d,w,kind,lev);op=s.prepare(ff,gg,dd,beta,ww);q=s.reference(op);truth=-1 if q<1-s.GUARD else 1 if q>1+s.GUARD else 0
                        fields={}
                        for method,vv in [('warm',v),('cold',np.ones(n))]:
                            t=time.perf_counter();ver,k,lo,hi,hist=s.screen(op,vv);sec=time.perf_counter()-t
                            assert not ver or ver==truth;assert lo<=q+1e-8 and q<=hi+1e-8
                            fields.update({method+'_steps':k,method+'_verdict':ver,method+'_s':sec})
                            for budget in [1,2,4,8,16]:fields[f'{method}_resolved_{budget}']=int(any(x[2] for x in hist[:budget]))
                        rows.append(dict(n=n,rep=cfg['rep'],h=h,period=op.period,target=target,lambda_old=actual,calibration_error=actual-target,kind=kind,level=lev,lambda_new=float(np.log(q)/op.period),**fields))
        print('duration',h,'rows',len(rows),flush=True)
        s.dump_csv('duration_sensitivity.csv',rows)
    # Interval action checks on fixed first three configs at n=6 and 12, five mutations each.
    ivrows=[];s.H=.125
    for cfg in [c for c in cfgs if c['n'] in [6,12] and c['rep']<3]:
        F=np.array(cfg['F']);G=np.array(cfg['G']);d=np.array(cfg['decay']);w=np.array(cfg['word']);n=cfg['n']
        beta=s.calibrate(F,G,d,w,-.002);_,v=s.reference(s.prepare(F,G,d,beta,w),True)
        for kind in s.KINDS:
            ff,gg,dd,ww=s.mutate(F,G,d,w,kind,.1);op=s.prepare(ff,gg,dd,beta,ww);q=s.reference(op);vv=v.copy()
            # Refine approximate vector, then independently validate action using interval formulas.
            for _ in range(16):vv=op.apply(vv);vv/=max(vv)
            ver,lo,hi=interval_action(ff,gg,dd,beta,ww,vv,s.H)
            assert ver==(-1 if q<1 else 1)
            ivrows.append(dict(n=n,rep=cfg['rep'],kind=kind,verdict=ver,lo=lo,hi=hi,rho_reference=q))
    (P/'interval_checks.json').write_text(json.dumps(ivrows,indent=2));(P/'independent_checks.json').write_text(json.dumps(checks,indent=2))
    (P/'extensions_metadata.json').write_text(json.dumps(dict(design='Exploratory extension chosen after main results; durations fixed before extension execution',n_duration_cases=len(rows),n_interval_checks=len(ivrows),n_independent_checks=len(checks),elapsed_s=time.perf_counter()-start,mpmath=mp.__version__),indent=2))
    print('DONE',len(rows),len(ivrows),len(checks),flush=True)
if __name__=='__main__':main()
