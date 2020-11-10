#!/usr/bin/python3

import numpy as np

from lmfit import Parameters, Minimizer, fit_report, conf_interval, report_ci
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from glob import glob
import numpy as np
from uncertainties import correlated_values
import uncertainties.unumpy as unp
from uncertainties import ufloat

###GAUSSIANA
def volume_fcn(x, a0, a1, a2, a3, a4, a5, a6, a7, a8):
    return a0 + np.exp(-(a1/x))*(a2/x**2.0 + a3/x + a4 + a5*x + a6*np.exp(-a7*(-a8+x)**2.0))

def CET_fcn(x, a0, a1, a2, a3, a4, a5, a6, a7, a8):
    dVdT = np.exp(-a1/x)*(a5-(2.0*a2)/x**3.0 - a3/x**2.0) + (a1*np.exp(-a1/x)*(a4+a2/x**2.0+a3/x+a5*x))/x**2.0
    CTE = dVdT/volume_fcn(x,a0, a1, a2, a3, a4, a5, a6, a7, a8)
    return CTE

def volume_fcn_unc(x, a0, a1, a2, a3, a4, a5, a6, a7, a8):
    return a0 + unp.exp(-(a1/x))*(a2/x**2.0 + a3/x + a4 + a5*x + a6*unp.exp(-a7*(-a8+x)**2.0))

def CET_fcn_unc(x, a0, a1, a2, a3, a4, a5, a6, a7, a8):
    dVdT = unp.exp(-a1/x)*(a5-(2.0*a2)/x**3.0 - a3/x**2.0) + (a1*unp.exp(-a1/x)*(a4+a2/x**2.0+a3/x+a5*x))/x**2.0
    CTE = dVdT/volume_fcn_unc(x,a0, a1, a2, a3, a4, a5, a6, a7, a8)
    return CTE

def func2min(params, T, V):
    parametros = [params[k] for k in params.keys()]
    model_vol = volume_fcn(T, *parametros)
    residual_vol = model_vol - V
    return residual_vol



def fit_volume_vs_temperature(T, V):
    
    params = Parameters()
    params.add_many(('a0', 400,  True,       0,   8000,   None, None),
                    ('a1',   100,  True,     0,   2000,   None, None),
                    ('a2',   0,  False,   -150000,   150000,   None, None),
                    ('a3',   0,  True,   -3000,   3000,   None, None),
                    ('a4',   0,  True,   -5,   5,   None, None),
                    ('a5',   0.001,  True,   0.00001,   0.15,   None, None),
                    ('a6',   0.0,  False,   -5,   0.0,   None, None),
                    ('a7',   0.001,  False,   0.00001,   0.001,   None, None),
                    ('a8',   0.001,  False,   0,   300,   None, None)
                    )
   
    
    minimizer = Minimizer(func2min, params, fcn_args=(T, V))
    
    kws={'popsize':60,  'mutation':(0.8,1.2), 'recombination':0.8, 'updating':'deferred', 'workers':-1}
    out = minimizer.minimize(method='differential_evolution', max_nfev=50000000,**kws)
    
    fit = func2min(out.params, T, V)
    
    print(fit_report(out), flush=True)
    print('Cost:',np.sum(fit))
    
    
    print('---------- DE bounds ----------')
    for k in out.params.keys():
        print('{0:3}:  {1: >8}\t{2: >12.6g}\t{3: >8}'.format(k,out.params[k].min,out.params[k].value,out.params[k].max))
        
    return out

def second_fit(T,V,out):
    minimizer = Minimizer(func2min, out.params, fcn_args=(T, V))
    result2 = minimizer.minimize(method='nelder')
    print(fit_report(result2), flush=True)
    return result2

################################################################################################


file = './volume_vs_temperature.dat'


print('{0:#^80}'.format(' '+str(file)+' '))
        
M = np.genfromtxt(file)
temperature = M[:,0]
volume = M[:,1]


result0 = fit_volume_vs_temperature(temperature, volume)
result = second_fit(temperature, volume, result0)

has_uncertainties = False

t_fit = np.array(temperature)

if not np.any(np.diag(result.covar) < 0):
    has_uncertainties = True        
else:
    print('!!!!!!!!!!!!!! Negative variances. Running Nelder-Mead.')
    result = second_fit(temperature, volume, result)
    if np.any(np.diag(result.covar) >= 0):
        has_uncertainties = True


(a0, a1, a2, a3, a4, a5, a6, a7, a8) = [result.params[k].value for k in result.params.keys()]

if has_uncertainties:
    try:
        (a0, a1, a3, a4, a5) = correlated_values([a0, a1, a3, a4, a5], result.covar)
    except Exception as e:
        print("The execption", e.__class__, "occurred when trying to construct correlated variables a0, a1, a4, a4, a5.")
        print("Ignoring correlation and treating variables as independent")
        a0 = ufloat(a0,np.sqrt(np.diag(result.covar)[0]))
        a1 = ufloat(a1,np.sqrt(np.diag(result.covar)[1]))
        a3 = ufloat(a3,np.sqrt(np.diag(result.covar)[2]))
        a4 = ufloat(a4,np.sqrt(np.diag(result.covar)[3]))
        a5 = ufloat(a5,np.sqrt(np.diag(result.covar)[4]))
    
    v_fit = volume_fcn_unc(t_fit, a0, a1, a2, a3, a4, a5, a6, a7, a8)
    v_nominal = np.array([x.nominal_value for x in v_fit])
    v_unc = np.array([x.std_dev for x in v_fit])
    
    cte_fit = CET_fcn_unc(t_fit, a0, a1, a2, a3, a4, a5, a6, a7, a8)
    cte_nominal = np.array([x.nominal_value for x in cte_fit])
    cte_unc = np.array([x.std_dev for x in cte_fit])
    
    
else:
    v_nominal = volume_fcn(t_fit,a0, a1, a2, a3, a4, a5, a6, a7, a8)
    v_unc = np.zeros(len(v_nominal))
    
    cte_nominal = CET_fcn_unc(t_fit,a0, a1, a2, a3, a4, a5, a6, a7, a8)
    cte_unc = np.zeros(len(v_nominal))
    

    

######SAVES DATA TO FILE
#Experimental linear thermal expansion coefficient
# Physical Review B 83, 104102 (2011)
b=3.6e-14
c=1.21e-11
T0=212
deltaT0=47
W = 1.0/(1+np.exp((t_fit-T0)/deltaT0))
y_exp=3*(b*W*t_fit**3.0 + c*t_fit**2.0*(1-W))
#Times 3 to get the volumetric TEC

with open('fitting_results_diamond.dat','w') as f:
    f.write('#{T:5}\t{V:10}\t{Vfit:10}\t{Vunc:10}\t{CET:10}\t{CETfit:10}\t{CETunc:10}\n'.format(T='T(K)',V=' V(A3)',Vfit=' Vfit(A3)',Vunc=' Vunc(A3)',CET=' CET(1/K)',CETfit=' CETfit(1/K)',CETunc=' CETunc(1/K)'))
    for i in range(len(temperature)):
        f.write('{T:5.1f}\t{V:10f}\t{Vfit:10f}\t{Vunc:10f}\t{CET: 10e}\t{CETfit: 10e}\t{CETunc: 10e}\n'.format(T=temperature[i],V=volume[i],Vfit=v_nominal[i],Vunc=v_unc[i],CET=y_exp[i],CETfit=cte_nominal[i],CETunc=cte_unc[i]))



######## PLOTS
fig,ax1 = plt.subplots(nrows=1, ncols=1,dpi=120,figsize=(6,4))

ax1.plot(temperature, volume, '.k',label='Data')
ax1.plot(t_fit, v_nominal, label='Fit', color='red')
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Volume ($\AA^3$)')
ax1.legend()
#plt.suptitle(file)
plt.savefig('volume_vs_temperature.pdf')


#########################################################################

fig,ax2 = plt.subplots(nrows=1, ncols=1,dpi=120,figsize=(6,4))

N = 5
cte_up = cte_nominal+N*cte_unc
cte_down = cte_nominal-N*cte_unc
#ax2.fill_between(t_fit, cte_up, cte_down, facecolor='#888888', alpha=0.5, label = f'{N}$\sigma$', edgecolor=None)
ax2.plot(t_fit, 1e6*y_exp, '.k',label='PRB 83, 104102 (2011)')
ax2.plot(t_fit, 1e6*cte_nominal, label='Fit', color='red')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('CET (10$^6$K$^{-1}$)')
ax2.legend()
#plt.suptitle(file)
plt.savefig('thermal_expansion_vs_temperature.pdf')
