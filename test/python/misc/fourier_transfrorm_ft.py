import numpy as np

ft = np.loadtxt("ftau_up.dat")
n_t = len(ft)
n_w=128
beta=20.0
dtau = beta/(n_t-1)

print('beta=', beta)

f=open('fw_from_t.dat','w')

for wn in range(n_w):
 iw = complex(0.,(2*wn+1)*np.pi/beta)
 s = complex(0.,0.)
 for tn in range(n_t):
   weight = 0.5 if(tn==0 or tn==n_t-1) else 1.0
   tau = tn*dtau
   s+= weight * ft[tn] * np.exp(iw*tau) * dtau
 f.write('%f %f\n'%(s.real,s.imag))

f.close()



