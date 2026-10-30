"""Focused edge cases for the screener, not implementation-mirroring unit tests."""
import json
import numpy as np
import study as s
from pathlib import Path

def main():
    results=[]
    for q,expected in [(.99,-1),(1.,0),(1.01,1)]:
        # Reducible diagonal system: bound validity does not require irreducibility.
        op=s.Operator(np.ones(6)*q**(1/16),np.zeros((3,3)),np.zeros((3,3)),np.array([0,1]*8))
        ver,k,lo,hi,hist=s.screen(op,np.arange(1,7,dtype=float))
        assert ver==expected,(q,ver,lo,hi)
        results.append(dict(test='diagonal_or_exact_boundary',rho=q,expected=expected,actual=ver))
    cfg=json.loads((Path(__file__).parent/'configs.json').read_text())[0]
    F=np.array(cfg['F']);G=np.array(cfg['G']);d=np.array(cfg['decay']);word=np.array(cfg['word']);beta=3.9
    op=s.prepare(F,G,d,beta,word);v=np.linspace(.5,2,6)
    _,_,_,_,a=s.screen(op,v);_,_,_,_,b=s.screen(op,100*v)
    assert np.allclose(a,b,rtol=1e-12,atol=1e-12)
    results.append(dict(test='positive_vector_scaling',passed=True))
    for scale in [.5,1.,1.5]:
        qo=s.reference(s.prepare(F,G,d,beta*scale,word))
        if scale==.5:prev=qo
        else:assert qo>prev;prev=qo
    results.append(dict(test='coupling_monotonicity',passed=True))
    (Path(__file__).parent/'edge_checks.json').write_text(json.dumps(results,indent=2));print('Passed',len(results),'edge checks')
if __name__=='__main__':main()
