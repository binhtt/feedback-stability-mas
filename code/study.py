"""Reproducible synthetic positive-system evolution study; numerical screening only.
Run with single-thread BLAS. No LLM calls, external data, or formal-certification claim.
"""
from pathlib import Path
from dataclasses import dataclass
import time,json,csv,platform,os
import numpy as np
from scipy.linalg import eig
from scipy.special import exprel
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator,eigs
import scipy

ROOT=Path(__file__).resolve().parent
SEED=202610290
SIZES=[6,12,24,48,96,192]
REPS=10
LEVELS=[.05,.1,.2]
KINDS=['weight','decay','order','joint','redistribute']
GUARD=1e-10
MAX_STEPS=16
H=.125

@dataclass
class Operator:
    d:np.ndarray
    X:np.ndarray
    Y:np.ndarray
    word:np.ndarray
    @property
    def n(self):return len(self.d)
    @property
    def period(self):return 2*H*len(self.word)
    def apply(self,v):
        z=v.copy();p=self.n//2
        for k in self.word:
            # off-diagonal block uses the input before diagonal scaling.
            if k==0:
                top=self.d[:p]*z[:p]+self.X@z[p:]
                z=np.concatenate((top,self.d[p:]*z[p:]))
            else:
                bottom=self.d[p:]*z[p:]+self.Y@z[:p]
                z=np.concatenate((self.d[:p]*z[:p],bottom))
        return z
    def dense(self):
        M=np.eye(self.n);p=self.n//2
        for k in self.word:
            if k==0:
                top=self.d[:p,None]*M[:p]+self.X@M[p:]
                M=np.vstack((top,self.d[p:,None]*M[p:]))
            else:
                bot=self.d[p:,None]*M[p:]+self.Y@M[:p]
                M=np.vstack((self.d[:p,None]*M[:p],bot))
        return M

def prepare(F,G,decay,beta,word):
    # Exact analytic phase formula in real arithmetic, evaluated in float64.
    # E_ij = beta F_ij h exp(-d_i h) exprel((d_i-d_j)h).
    p=len(decay)//2
    a=np.exp(-decay*H)
    cross=H*a[:,None]*exprel(H*(decay[:,None]-decay[None,:]))
    X=beta*F*cross[:p,p:]*a[:p,None]  # subsequent idle interval
    Y=beta*G*cross[p:,:p]*a[p:,None]
    return Operator(a*a,X,Y,np.array(word,copy=True))

def reference(op,vector=False):
    M=op.dense()
    if vector:
        vals,vecs=eig(M);j=np.argmax(np.abs(vals));q=float(abs(vals[j]));v=np.real(vecs[:,j]);v*=np.sign(v.sum());assert np.all(v>0);v/=v.max()
        return q,v
    return float(max(abs(np.linalg.eigvals(M))))

def screen(op,v,steps=MAX_STEPS):
    v=v.copy()/max(v);lo=0.;hi=np.inf;history=[]
    for k in range(1,steps+1):
        y=op.apply(v);rat=y/v
        lo=max(lo,float(min(rat)));hi=min(hi,float(max(rat)))
        assert lo<=hi+1e-11
        verdict=-1 if hi<1-GUARD else 1 if lo>1+GUARD else 0
        history.append((lo,hi,verdict))
        if verdict:return verdict,k,lo,hi,history
        v=y/max(y)
    return 0,steps,lo,hi,history

def base_model(rng,n):
    p=n//2;F=rng.uniform(.1,1,(p,p));G=rng.uniform(.1,1,(p,p))
    F/=np.mean(F.sum(axis=1));G/=np.mean(G.sum(axis=1));d=rng.uniform(.8,1.2,n)
    w=np.array([0]*8+[1]*8);rng.shuffle(w)
    return F,G,d,w

def mutate(F,G,d,w,kind,level):
    f=F.copy();g=G.copy();dd=d.copy();ww=w.copy()
    if kind in ['weight','joint']:f[0]*=1+level
    if kind=='joint':g[0]*=1+level
    if kind=='decay':dd[0]*=1-level
    if kind=='redistribute':f[0]=(1-level)*f[0]+level*np.roll(f[0],1)
    if kind=='order':
        for _ in range({.05:1,.1:2,.2:4}[level]):
            ys=np.where(ww==1)[0];xs=np.where(ww==0)[0]
            if ys[0]>=xs[-1]:break
            ww[ys[0]],ww[xs[-1]]=0,1
    return f,g,dd,ww

def calibrate(F,G,d,w,target):
    # Matrix-free power root estimation. Dense reference validates resulting margin.
    def fun(beta):
        op=prepare(F,G,d,beta,w);v=np.ones(len(d));q=1.
        for _ in range(100):
            y=op.apply(v);q=float(max(y));nv=y/q
            if max(abs(nv-v))<1e-13:break
            v=nv
        return np.log(q)/op.period-target
    return brentq(fun,0,32,xtol=1e-11)

def dump_csv(name,rows):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)

def main():
    rng=np.random.default_rng(SEED);rows=[];bases=[];configs=[];bench=[];failures=[];start=time.perf_counter()
    for n in SIZES:
        nstart=time.perf_counter()
        for rep in range(REPS):
            F,G,d,w=base_model(rng,n)
            corebeta=float(rng.uniform(1,6))
            configs.append(dict(n=n,rep=rep,F=F.tolist(),G=G.tolist(),decay=d.tolist(),word=w.tolist(),core_beta=corebeta))
            for regime,target in [('core',None),('far',-.2),('near',-.02),('critical',-.002)]:
                beta=corebeta if target is None else calibrate(F,G,d,w,target)
                old=prepare(F,G,d,beta,w);t=time.perf_counter();qold,vold=reference(old,True);offline=time.perf_counter()-t
                lamold=float(np.log(qold)/old.period)
                if target is not None:assert abs(lamold-target)<1e-9,(n,rep,regime,lamold)
                bases.append(dict(n=n,rep=rep,regime=regime,beta=beta,lambda_old=lamold,offline_reference_s=offline))
                for kind in KINDS:
                    for level in LEVELS:
                        ff,gg,dd,ww=mutate(F,G,d,w,kind,level)
                        t=time.perf_counter();op=prepare(ff,gg,dd,beta,ww);prep=time.perf_counter()-t
                        t=time.perf_counter();q=reference(op);dense_s=time.perf_counter()-t
                        truth=-1 if q<1-GUARD else 1 if q>1+GUARD else 0
                        fields={}
                        for method,v in [('warm',vold),('cold',np.ones(n))]:
                            t=time.perf_counter();ver,k,lo,hi,hist=screen(op,v);sec=time.perf_counter()-t
                            if ver and ver!=truth:failures.append(dict(n=n,rep=rep,regime=regime,kind=kind,level=level,method=method,q=q,lo=lo,hi=hi))
                            assert lo<=q+1e-8 and q<=hi+1e-8,(q,lo,hi)
                            fields.update({method+'_verdict':ver,method+'_steps':k,method+'_s':sec,method+'_lo':lo,method+'_hi':hi})
                            for budget in [1,2,4,8,16]:
                                fields[f'{method}_resolved_{budget}']=int(any(h[2] for h in hist[:budget]))
                        row=dict(n=n,rep=rep,regime=regime,kind=kind,level=level,beta=beta,lambda_old=lamold,lambda_new=float(np.log(q)/op.period),truth=truth,stable_to_unstable=int(qold<1-GUARD and q>1+GUARD),prep_s=prep,dense_s=dense_s,**fields)
                        rows.append(row)
                        # Dedicated timing suite: repeat timings, alternate order to limit drift.
                        if rep<3 and regime=='near' and level==.1 and kind in ['weight','joint','order']:
                            times={m:[] for m in ['dense','warm','cold','arpack']};arerrs=[]
                            for repeat in range(7):
                                methods=list(times)
                                if repeat%2:methods.reverse()
                                for method in methods:
                                    t=time.perf_counter()
                                    if method=='dense':reference(op)
                                    elif method in ['warm','cold']:screen(op,vold if method=='warm' else np.ones(n))
                                    else:
                                        arq=eigs(LinearOperator((n,n),matvec=op.apply,dtype=float),k=1,which='LM',v0=np.ones(n),tol=1e-10,maxiter=500,return_eigenvectors=False)[0]
                                        arerrs.append(float(abs(abs(arq)-q)))
                                    times[method].append(time.perf_counter()-t)
                            bench.append(dict(n=n,rep=rep,kind=kind,warm_verdict=fields['warm_verdict'],cold_verdict=fields['cold_verdict'],**{m+'_median_s':float(np.median(t)) for m,t in times.items()},arpack_max_rho_error=max(arerrs)))
        dump_csv('cases.csv',rows);dump_csv('baselines.csv',bases);dump_csv('timings.csv',bench)
        print(f'n={n}: cases={len(rows)}, elapsed={time.perf_counter()-start:.1f}s (size {time.perf_counter()-nstart:.1f}s)',flush=True)
    (ROOT/'configs.json').write_text(json.dumps(configs))
    (ROOT/'failures.json').write_text(json.dumps(failures,indent=2))
    assert not failures,failures
    meta=dict(seed=SEED,sizes=SIZES,seeds_per_size=REPS,n_cases=len(rows),n_baselines=len(bases),levels=LEVELS,kinds=KINDS,max_steps=MAX_STEPS,guard=GUARD,h=H,period=4.,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,elapsed_s=time.perf_counter()-start,threads={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']},numerical_scope='Float64 CW inequalities with guard, NOT outward-rounded certification',design='Fixed before new study run; pilot informed choices; core uncalibrated, stress groups calibrated; no seed filtering')
    (ROOT/'metadata.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__':main()
