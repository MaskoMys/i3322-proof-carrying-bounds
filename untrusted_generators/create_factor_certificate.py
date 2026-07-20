from __future__ import annotations
from fractions import Fraction
import pickle, numpy as np, json, hashlib, sys
sys.path.append('/mnt/data/i3322_project')
from solve_sdp import bell_coeffs, star_mul, E
from exact_cert import residual, word_to_json, frac_to_str

def create(pkl,out,bits=40):
    with open(pkl,'rb') as f: res=pickle.load(f)
    W,Y,b,*_=res;Y=(Y+Y.T)/2
    vals,V=np.linalg.eigh(Y)
    keep=vals>0
    vals=vals[keep];V=V[:,keep]
    Rf=V*np.sqrt(vals)
    den=1<<bits
    R=[[Fraction(int(round(float(Rf[i,j])*den)),den) for j in range(Rf.shape[1])] for i in range(Rf.shape[0])]
    n=len(W);r=len(R[0])
    Yq=[[sum(R[i][k]*R[j][k] for k in range(r)) for j in range(n)] for i in range(n)]
    beta=Fraction(str(b))
    Res,l1=residual(W,Yq,beta)
    certified=beta+l1
    data={
      'format':'i3322-level3-factor-absorption-v1','word_length':3,
      'words':[word_to_json(w) for w in W], 'rank':r,'factor_bits':bits,
      'beta':frac_to_str(beta),
      'R':[[frac_to_str(x) for x in row] for row in R],
      'residual_l1':frac_to_str(l1),
      'certified_beta':frac_to_str(certified),
      'certified_I':frac_to_str((certified-Fraction(4))/4),
      'residual_terms':len(Res),
      'metadata':{'source_pickle':pkl,'numerical_beta':b,'numerical_min_eigenvalue':float(np.min(np.linalg.eigvalsh(Y)))}
    }
    payload=json.dumps(data,sort_keys=True,separators=(',',':')).encode();data['sha256_without_hash']=hashlib.sha256(payload).hexdigest()
    with open(out,'w') as f:json.dump(data,f,indent=2,sort_keys=True)
    print('rank',r,'l1',float(l1),'certI',float((certified-4)/4),'size',len(payload))
    return data
if __name__=='__main__':create('/mnt/data/i3322_project/k3_margin_5_003505.pkl','/mnt/data/i3322_project/certificate_level3_factor.json',40)
