import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import (figure, axes, plot, xlabel, ylabel, title,
     grid, savefig, show)
import math
from scipy import signal
from scipy.signal import lfilter,dimpulse,impulse
from pandas import Series
from matplotlib import pyplot
import itertools
# ----------------------------
from boxsim import boxsim
from boxres import boxres
from auto_v import auto_v
from cal_teta import cal_teta
from cal_parameters import cal_parameters

def Levenberg_marquadrt(u,y,teta,nb,nc,nd,nf,NP,delta,mu,gamma,D,N1):
    VV_new = np.zeros((1, 20))
    VV_new = VV_new[0, :]
    for i in range(0,20): # This is the place where the sum square errors goes down

        mu = mu/gamma
        e= boxres(teta,nb,nc,nd,nf,u,y,D)

        e1 = e.reshape(-1, e.shape[0])
        e2 = e1.T
        V = np.matmul(e1 , e2)
        VV_new[i] = V
        teta_dum=teta.copy()
        J_t_J, J_t_e, V = cal_teta(u, y, teta,teta_dum, nb, nc, nd, nf, delta, NP, D,N1)
        #output=(J_t_J, J_t_e, V)

        while 1:  #loop with going up the sum square error
            delta_teta, teta_new, V_new, e = cal_parameters(J_t_J, teta_dum, nb, nc, nd, nf, J_t_e, mu, u, y, D)
            #out=[delta_teta, teta_new, V_new, e]

            if V_new < V or (mu > 1e10):

               break
            elif V_new  > V or np.isnan(V_new ) or np.isinf(V_new ):
                mu = mu * gamma  #Increase mu by 10  (% End of bad Loop(Error goes up))
        teta = teta_new

    return[teta,VV_new,J_t_J]