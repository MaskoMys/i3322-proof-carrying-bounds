from __future__ import annotations
import sys,pickle
sys.path.append('/mnt/data/i3322_project')
from solve_sdp import words, bell_coeffs, star_mul, E
import numpy as np, cvxpy as cp

def sigma_basis(W):
    idx={w:i for i,w in enumerate(W)};seen=set();plus=[];minus=[]
    for i,w in enumerate(W):
        if i in seen: continue
        sw=(w[1],w[0]); j=idx[sw]; s=(-1)**(len(w[0])+len(w[1]))
        if i==j:
            assert s==1; v=np.zeros(len(W));v[i]=1;plus.append(v);seen.add(i)
        else:
            assert i!=j
            v=np.zeros(len(W));v[i]=1;v[j]=s;plus.append(v)
            v=np.zeros(len(W));v[i]=1;v[j]=-s;minus.append(v)
            seen.add(i);seen.add(j)
    return np.array(plus).T,np.array(minus).T

def solve_sym(k,beta_target=None,margin=False,eps=1e-8,solver='CLARABEL',verbose=False):
    W=words(k);n=len(W);Tp,Tm=sigma_basis(W);npdim=Tp.shape[1];nmdim=Tm.shape[1]
    Zp=cp.Variable((npdim,npdim),symmetric=True);Zm=cp.Variable((nmdim,nmdim),symmetric=True)
    Y=Tp@Zp@Tp.T+Tm@Zm@Tm.T
    beta=cp.Variable() if beta_target is None else beta_target
    cons=[Zp>>0,Zm>>0]
    t=None
    if margin:
        t=cp.Variable();cons=[Zp-t*np.eye(npdim)>>0,Zm-t*np.eye(nmdim)>>0]
    pairs={}
    for i,u in enumerate(W):
      for j,v in enumerate(W): pairs.setdefault(star_mul(u,v),[]).append((i,j))
    B=bell_coeffs();allg=set(pairs)|set(B)|{E}
    for g in allg:
      expr=cp.sum([Y[i,j] for i,j in pairs.get(g,[])])
      target=(beta if g==E else 0)-float(B.get(g,0))
      cons.append(expr==target)
    obj=cp.Maximize(t) if margin else cp.Minimize(beta)
    prob=cp.Problem(obj,cons)
    
    if solver=='CLARABEL': kwargs={'tol_gap_abs':eps,'tol_feas':eps,'tol_gap_rel':eps,'max_iter':3000,'verbose':verbose}
    elif solver=='SCS': kwargs={'eps':eps,'max_iters':1000000,'verbose':verbose}
    elif solver=='CVXOPT': kwargs={'abstol':eps,'reltol':eps,'feastol':eps,'max_iters':500,'kktsolver':'robust','verbose':verbose}
    else: kwargs={'verbose':verbose}
    prob.solve(solver=solver,**kwargs)
    Yv=Tp@Zp.value@Tp.T+Tm@Zm.value@Tm.T
    bval=beta if isinstance(beta,(float,int)) else beta.value
    print('k',k,'dims',n,npdim,nmdim,'status',prob.status,'beta',bval,'I',(bval-4)/4,'t',None if t is None else t.value,'eigmin',np.linalg.eigvalsh((Yv+Yv.T)/2)[0])
    # residual
    rs=[]
    for g in allg:
      val=sum(Yv[i,j] for i,j in pairs[g]);target=(bval if g==E else 0)-float(B.get(g,0));rs.append(abs(val-target))
    print('res max/sum',max(rs),sum(rs),len(rs))
    return W,Yv,bval,Tp,Tm,Zp.value,Zm.value,pairs

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--k',type=int,default=3);ap.add_argument('--target',type=float);ap.add_argument('--margin',action='store_true');ap.add_argument('--out')
 a=ap.parse_args();res=solve_sym(a.k,a.target,a.margin,eps=1e-9)
 if a.out:
  with open(a.out,'wb') as f:pickle.dump(res,f)
