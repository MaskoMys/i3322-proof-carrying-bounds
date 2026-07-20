from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Tuple, List, Dict
import cvxpy as cp
import numpy as np

# Word represented by pair (aword,bword), each tuple of ints 1..3 reduced adjacent unequal
Word = Tuple[Tuple[int,...], Tuple[int,...]]
E: Word = ((),())

def reduce_single(seq):
    st=[]
    for x in seq:
        if st and st[-1]==x:
            st.pop()
        else:
            st.append(x)
    return tuple(st)

def mul(u:Word,v:Word)->Word:
    return (reduce_single(u[0]+v[0]), reduce_single(u[1]+v[1]))

def inv(u:Word)->Word:
    return (tuple(reversed(u[0])), tuple(reversed(u[1])))

def star_mul(u:Word,v:Word)->Word:
    return mul(inv(u),v)

def single_words(maxlen):
    out=[()]
    def rec(pref,n):
        if n==0:
            out.append(tuple(pref)); return
        for i in (1,2,3):
            if not pref or pref[-1]!=i:
                rec(pref+[i],n-1)
    for l in range(1,maxlen+1): rec([],l)
    return out

def words(k:int)->List[Word]:
    sw=single_words(k)
    out=[]
    for a in sw:
        for b in sw:
            if len(a)+len(b)<=k:
                out.append((a,b))
    out.sort(key=lambda w:(len(w[0])+len(w[1]),len(w[0]),w[0],w[1]))
    return out

def bell_coeffs()->Dict[Word,Fraction]:
    c={}
    def add(w,val): c[w]=c.get(w,Fraction(0))+Fraction(val)
    # A1+A2-B1-B2
    add(((1,),()),1); add(((2,),()),1)
    add(((),(1,)),-1); add(((),(2,)),-1)
    # correlations
    for j in (1,2,3): add(((1,),(j,)),1)
    add(((2,),(1,)),1);add(((2,),(2,)),1);add(((2,),(3,)),-1)
    add(((3,),(1,)),1);add(((3,),(2,)),-1)
    return c

def solve_dual(k:int, solver='CLARABEL', eps=1e-9, verbose=False):
    W=words(k); n=len(W)
    pairs={}
    for i,u in enumerate(W):
        for j,v in enumerate(W):
            g=star_mul(u,v)
            pairs.setdefault(g,[]).append((i,j))
    allg=set(pairs)|set(bell_coeffs())|{E}
    Y=cp.Variable((n,n), symmetric=True)
    beta=cp.Variable()
    constraints=[Y >> 0]
    B=bell_coeffs()
    for g in allg:
        expr=cp.sum([Y[i,j] for i,j in pairs.get(g,[])]) if g in pairs else 0
        target=(beta if g==E else 0) - float(B.get(g,0))
        constraints.append(expr==target)
    prob=cp.Problem(cp.Minimize(beta),constraints)
    kwargs={}
    if solver=='SCS': kwargs={'eps':eps,'max_iters':1000000,'verbose':verbose}
    elif solver=='CLARABEL': kwargs={'tol_gap_abs':eps,'tol_feas':eps,'tol_gap_rel':eps,'max_iter':1000,'verbose':verbose}
    prob.solve(solver=solver,**kwargs)
    print('k',k,'n',n,'status',prob.status,'beta',beta.value,'I',(beta.value-4)/4)
    vals=np.linalg.eigvalsh((Y.value+Y.value.T)/2)
    print('eig min/max',vals[0],vals[-1], 'rank >1e-7',np.sum(vals>1e-7))
    return W,Y.value,beta.value,pairs

if __name__=='__main__':
    for k in [1,2,3]:
        solve_dual(k, solver='CLARABEL',eps=1e-10)

def solve_margin(k:int,beta_target:float, solver='CLARABEL', eps=1e-10, verbose=False):
    W=words(k); n=len(W)
    pairs={}
    for i,u in enumerate(W):
        for j,v in enumerate(W):
            g=star_mul(u,v); pairs.setdefault(g,[]).append((i,j))
    allg=set(pairs)|set(bell_coeffs())|{E}
    Y=cp.Variable((n,n), symmetric=True); t=cp.Variable()
    constraints=[Y-t*np.eye(n) >> 0]
    B=bell_coeffs()
    for g in allg:
        expr=cp.sum([Y[i,j] for i,j in pairs.get(g,[])]) if g in pairs else 0
        target=(beta_target if g==E else 0) - float(B.get(g,0))
        constraints.append(expr==target)
    prob=cp.Problem(cp.Maximize(t),constraints)
    kwargs={'tol_gap_abs':eps,'tol_feas':eps,'tol_gap_rel':eps,'max_iter':2000,'verbose':verbose} if solver=='CLARABEL' else {'eps':eps,'max_iters':1000000,'verbose':verbose}
    prob.solve(solver=solver,**kwargs)
    vals=np.linalg.eigvalsh((Y.value+Y.value.T)/2)
    print('margin k',k,'target',beta_target,'status',prob.status,'t',t.value,'eigmin',vals[0])
    # constraint residual max and l1
    res=[]
    for g in allg:
        val=sum(Y.value[i,j] for i,j in pairs.get(g,[]))
        target=(beta_target if g==E else 0)-float(B.get(g,0))
        res.append(abs(val-target))
    print('constraint max',max(res),'sum',sum(res),'num',len(res))
    return W,Y.value,float(t.value),pairs
