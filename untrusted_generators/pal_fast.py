from __future__ import annotations
import numpy as np
from scipy.linalg import eigh_tridiagonal

def iterate_fast(n:int, cn:float=-1.0, C:float=0.9, maxit=10000,tol=1e-15):
    c=np.empty(n+1,dtype=float); c[0]=1.0;c[n]=cn
    c[1:n]=np.where(np.arange(1,n)<n/2,C,-C if cn<0 else C)
    lam=np.ones(n)/np.sqrt(n)
    last=-1e9
    for it in range(maxit):
        diag=c[:-1]*c[1:]+(c[:-1]-c[1:])/2-1
        diag[-1]+=(c[n]+1)/2
        off=np.sqrt(np.maximum(0,1-c[1:n]**2))/2
        vals,vecs=eigh_tridiagonal(diag,off,select='i',select_range=(n-1,n-1),check_finite=False)
        lam=np.abs(vecs[:,0]); lam/=np.linalg.norm(lam)
        tau=(1+2*c[2:n+1])*lam[1:]**2-(1-2*c[0:n-1])*lam[:-1]**2
        den=np.sqrt(tau*tau+4*lam[:-1]**2*lam[1:]**2)
        cnew=tau/den
        delta=np.max(np.abs(c[1:n]-cnew))
        c[1:n]=cnew
        val=float(vals[0])
        if abs(val-last)<tol and delta<tol: break
        last=val
    diag=c[:-1]*c[1:]+(c[:-1]-c[1:])/2-1;diag[-1]+=(c[n]+1)/2
    off=np.sqrt(np.maximum(0,1-c[1:n]**2))/2
    vals,vecs=eigh_tridiagonal(diag,off,select='i',select_range=(n-1,n-1),check_finite=False)
    lam=np.abs(vecs[:,0]);lam/=np.linalg.norm(lam)
    return float(vals[0]),c,lam,it+1

if __name__=='__main__':
 import time
 for n in [99,199,499,999,1999,4999]:
  t=time.time();v,c,l,it=iterate_fast(n,-1.0,tol=1e-14,maxit=5000);print(n,v,0.250875384514-v,it,time.time()-t)
