from __future__ import annotations
from fractions import Fraction
import json,hashlib,math,sys
sys.set_int_max_str_digits(1000000)
sys.path.append('/mnt/data/i3322_project')
from pal_fast import iterate_fast

def dyad(x,b):
    den=1<<b;return Fraction(int(round(float(x)*den)),den)

def rat_cs(c,b):
    s=math.sqrt(max(0.0,1-c*c)); t=dyad(s/(1+c),b) if c>-1+1e-15 else None
    if t is None:return Fraction(-1),Fraction(0)
    return (1-t*t)/(1+t*t),2*t/(1+t*t)

def fs(x):return f'{x.numerator}/{x.denominator}'

def exact_value(c,s,lam,cn):
    n=len(lam); norm=sum(x*x for x in lam)
    # Eq 23 tridiagonal M, exact
    num=Fraction(0)
    for i0 in range(n):
        i=i0+1
        diag=c[i-1]*c[i]+(c[i-1]-c[i])/2-Fraction(1)
        if i==n:diag+=(c[n]+1)/2
        num+=diag*lam[i0]*lam[i0]
        if i<n:num+=s[i]*lam[i0]*lam[i0+1]  # 2*(s_i/2)
    return num/norm,norm,num

def create(n=499,cn=-1.0,cbits=32,lbits=44,out='lower_compact.json'):
    val,cnum,lamnum,it=iterate_fast(n,cn,tol=1e-15,maxit=10000)
    c={0:Fraction(1),n:Fraction(int(cn))};s={}
    for i in range(1,n):c[i],s[i]=rat_cs(cnum[i],cbits)
    lam=[dyad(x,lbits) for x in lamnum]
    I,norm,num=exact_value(c,s,lam,Fraction(int(cn)))
    data={'format':'i3322-pal-rational-lower-compact-v1','dimension':n,'cn':int(cn),'c_bits':cbits,'state_bits':lbits,
          'numerical_I':val,'iterations':it,
          'parameters':[{'i':i,'c':fs(c[i]),'s':fs(s[i])} for i in range(1,n)],
          'state':[fs(x) for x in lam],
          'target_I':'2508753844/10000000000',
          'target_B':'12508753844/2500000000',
          'exact_value_above_target': bool(I > Fraction(2508753844,10000000000)),
          'construction':'Pal-Vertesi block projectors; manuscript mapping A=(Pal A2,A1,A3), B=(Pal B2,B1,B3); observables=2P-I'}
    payload=json.dumps(data,sort_keys=True,separators=(',',':')).encode();data['sha256_without_hash']=hashlib.sha256(payload).hexdigest()
    with open(out,'w') as f:json.dump(data,f,indent=2,sort_keys=True)
    print('n',n,'numerical',val,'exact',float(I),'loss',val-float(I),'B',float(4+4*I),'size approx',len(payload))
    return data
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=499);ap.add_argument('--cbits',type=int,default=32);ap.add_argument('--lbits',type=int,default=44);ap.add_argument('--out',default='lower_compact.json')
 a=ap.parse_args();create(a.n,-1.0,a.cbits,a.lbits,a.out)
