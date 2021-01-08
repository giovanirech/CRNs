import numpy as np

from lmfit import Parameters, Minimizer, fit_report, conf_interval, report_ci
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp

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


def CET_numerical(T,V):
    #linear extrapolation of volume at T=0
    v0 = V[0]-((V[1]-V[0])/(T[1]-T[0]))*T[0]
    T = np.insert(T, 0, 0.0)
    V = np.insert(V, 0, v0)
    CET = [0.0]
    for i in range(1,len(T)):
        CET.append((1.0/V[i])*(V[i]-V[i-1])/(T[i]-T[i-1]))
    return np.array(CET[1:])

#def func2min(params, T, V, CET):
#    parametros = [params[k] for k in params.keys()]
#    model_vol = volume_fcn(T, *parametros)
#    weights_volume = np.ones(len(T))
#    #weights_volume = np.log(temperature)
#    residual_vol = ((model_vol - V)/(V[-1]-V[0]))*weights_volume
#    
#    #model_cte = CET_fcn(T,*parametros)
#    #weights_cte = np.ones(len(T))
#    #weights_cte=np.log(temperature)
#    #residual_cte = ((model_cte - CET)/(CET[-1]-CET[0]))*weights_cte
#    
#    #residuals = np.concatenate((residual_vol, residual_cte))
#    return residual_vol

def func2min(params, T, V, CET):
    parametros = [params[k] for k in params.keys()]
    model_vol = volume_fcn(T, *parametros)
    residual_vol = model_vol - V  
    return residual_vol

def fit_volume_vs_temperature(T, V, CET):
    
    params = Parameters()
    #GAUSSIANA
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
    #EXPANDIDA
    #params.add_many(('a0', 400,  True,       0,   8000,   None, None),
    #                ('a1',   416.5,  True,     0,   2000,   None, None),
    #                ('a2',   0,  True,   -150000,   150000,   None, None),
    #                ('a3',   139,  True,   -10000,   10000,   None, None),
    #                ('a4',   -1.24,  True,   -50,   50,   None, None),
    #                ('a5',   0.001,  True,   0.0001,   0.15,   None, None),
    #                ('a6',   0.0,  True,   -1e-5,   1e-5,   None, None),
    #                )
    
    
    
    minimizer = Minimizer(func2min, params, fcn_args=(T, V, CET))
    
    #kws={'popsize':60,  'mutation':(0.8,1.2), 'recombination':0.8, 'updating':'deferred', 'workers':-1}
    #out = minimizer.minimize(method='differential_evolution', max_nfev=50000000,**kws)
    
    #kws = {'local':'L-BFGS-B', 'totaliter': 50, 'maxiter':20}
    kws = {'local':'L-BFGS-B', 'totaliter': 20, 'maxiter':5}
    out = minimizer.minimize(method='ampgo', **kws)
    
    fit = func2min(out.params, T, V,CET)
    
    print(fit_report(out), flush=True)
    print('Cost:',np.sum(fit))
    
    
    print('---------- DE bounds ----------')
    for k in out.params.keys():
        print('{0:3}:  {1: >8}\t{2: >12.6g}\t{3: >8}'.format(k,out.params[k].min,out.params[k].value,out.params[k].max))
        
    return out

def second_fit(T,V,CET,out):
    minimizer = Minimizer(func2min, out.params, fcn_args=(T, V, CET))
    result2 = minimizer.minimize(method='nelder')
    ci = conf_interval(minimizer, result2)
    report_ci(ci)
    print(fit_report(result2), flush=True)
    return result2

def plot(x_data, y_data, x_fit, y_fit, z_data, z_fit, file):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,dpi=120,figsize=(10,3))
    ax1.plot(x_data, y_data, 'o', label='data')
    ax1.plot(x_fit, y_fit, '-',label='fit')
    ax1.set_ylabel('Volume (A3)')
    ax1.set_xlabel('Temperature (K)')
    ax1.legend()
    
    ax2.plot(x_data, 1e6*z_data, 'x',label='Data')
    ax2.plot(x_fit, 1e6*z_fit, label='Fit')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('CET (K$^{-6}$)')
    ax2.legend()
    plt.suptitle(file)
    plt.savefig(file.rstrip('dat')+'png')
    #plt.show()
    
def plot_unc(x_data, y_data, x_fit, y_fit, y_unc, z_data, z_fit,z_unc, file):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,dpi=120,figsize=(10,3))
    #uncertainty band
    y_up = y_fit+y_unc
    y_down = y_fit-y_unc
    ax1.fill_between(x_fit, y_up, y_down, facecolor='#888888', alpha=0.25, label = '1$\sigma$', edgecolor=None)
    ax1.plot(x_data, y_data, '.', label='data')
    ax1.plot(x_fit, y_fit, '-',label='fit')
    ax1.set_ylabel('Volume (A3)')
    ax1.set_xlabel('Temperature (K)')
    ax1.legend()
    
    z_up = z_fit+z_unc
    z_down = z_fit-z_unc
    ax2.fill_between(x_fit, z_up, z_down, facecolor='#888888', alpha=0.25, label = '1$\sigma$', edgecolor=None)
    ax2.plot(x_data, 1e6*z_data, 'x',label='Data')
    ax2.plot(x_fit, 1e6*z_fit, label='Fit')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('CET (K$^{-6}$)')
    ax2.legend()
    plt.suptitle(file)
    plt.savefig(file.rstrip('dat')+'png')
    #plt.show()
    
def plot2(x_data, y_data, x_fit, y_fit, y_fit2, z_data, z_fit, z_fit2, file):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,dpi=120,figsize=(10,3))
    ax1.plot(x_data, y_data, 'o', label='data')
    ax1.plot(x_fit, y_fit, ':',label='fit ampgo')
    ax1.plot(x_fit, y_fit2, '-',label='fit nm')
    ax1.set_ylabel('Volume (A3)')
    ax1.set_xlabel('Temperature (K)')
    ax1.legend()
    
    ax2.plot(x_data, 1e6*z_data, 'x',label='Data')
    ax2.plot(x_fit, 1e6*z_fit, ':', label='Fit ampgo')
    ax2.plot(x_fit, 1e6*z_fit2,'-', label='fit nm')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('CET (K$^{-6}$)')
    ax2.legend()
    plt.suptitle(file)
    plt.savefig(file.rstrip('dat')+'png')
    #plt.show()
    
##########################################################################

from glob import glob
import numpy as np
from uncertainties import correlated_values
import uncertainties.unumpy as unp
from uncertainties import ufloat

X_fit = []
Y_fit = []
Z_fit = []
X_data = []
Y_data = []
Z_data = []

shrinks = [4, 6, 8, 10]
labels = ['4x4x4', '6x6x6', '8x8x8', '10x10x10']


print('{0:#^100}'.format(' Shrink '))

for shrink in shrinks:
    
    file = './volume_vs_temperature{0}.dat'.format(shrink)
    
    print('{0:#^80}'.format(' '+str(file)+' '))
            
    M = np.genfromtxt(file)
    temperature = M[:,0]
    volume = M[:,1]
    
    cte_data = CET_numerical(temperature, volume)
    
    result0 = fit_volume_vs_temperature(temperature, volume, cte_data)
    result = second_fit(temperature, volume, cte_data, result0)
    
    has_uncertainties = False
    
    t_fit = np.array(temperature)
    
    if not np.any(np.diag(result.covar) < 0):
        has_uncertainties = True        
    else:
        print('!!!!!!!!!!!!!! Negative variances. Running Nelder-Mead.')
        result = second_fit(temperature, volume, cte_data, result)
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
        
        plot_unc(temperature, volume, t_fit, v_nominal, v_unc, cte_data, cte_nominal, cte_unc, file)
        
    else:
        v_nominal = volume_fcn(t_fit,a0, a1, a2, a3, a4, a5, a6, a7, a8)
        v_unc = np.zeros(len(v_nominal))
        
        cte_nominal = CET_fcn_unc(t_fit,a0, a1, a2, a3, a4, a5, a6, a7, a8)
        cte_unc = np.zeros(len(v_nominal))
        
        plot(temperature, volume, t_fit, v_nominal, cte_data, cte_nominal, file)
            
    with open(file.rstrip('dat')+'fit','w') as f:
        f.write('#{T:5}\t{V:10}\t{Vfit:10}\t{Vunc:10}\t{CET:10}\t{CETfit:10}\t{CETunc:10}\n'.format(T='T(K)',V='V(A3)',Vfit='Vfit(A3)',Vunc='Vunc(A3)',CET='CET(1/K)',CETfit='CETfit(1/K)',CETunc='CETunc(1/K)'))
        for i in range(len(temperature)):
            f.write('{T:5.1f}\t{V:10f}\t{Vfit:10f}\t{Vunc:10f}\t{CET: 10e}\t{CETfit: 10e}\t{CETunc: 10e}\n'.format(T=temperature[i],V=volume[i],Vfit=v_nominal[i],Vunc=v_unc[i],CET=cte_data[i],CETfit=cte_nominal[i],CETunc=cte_unc[i]))
    
    X_data.append(temperature)
    Y_data.append(volume)
    Z_data.append(cte_data)
    
    X_fit.append(t_fit)
    Z_fit.append(cte_nominal)
    Y_fit.append(v_nominal)
                  

fig,ax1 = plt.subplots(nrows=1, ncols=1,dpi=180,figsize=(6,4))

for i in range(len(shrinks)):
    temperature = X_data[i]
    volume = Y_data[i]
    t_fit = X_fit[i]
    v_nominal = Y_fit[i]
    l, = ax1.plot(temperature, volume, '.')
    ax1.plot(t_fit, v_nominal, label=labels[i], color=l.get_color())
    
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Volume ($\AA^3$)')
ax1.legend()
plt.tight_layout()
plt.savefig('volume_vs_temp_all_shrinks.png')

#########################################################################

fig,ax2 = plt.subplots(nrows=1, ncols=1,dpi=120,figsize=(6,4))

for i in range(len(shrinks)):
    y_exp = Z_data[i]
    temp = X_data[i]
    t_fit = X_fit[i]
    cte_nominal = Z_fit[i]
    
    l, = ax2.plot(temp, 1e6*y_exp, '.')
    ax2.plot(t_fit, 1e6*cte_nominal, label=labels[i], color=l.get_color())
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('CET (10$^6$K$^{-1}$)')
ax2.legend()
plt.tight_layout()  
plt.savefig('CET_vs_temp_all_shrinks.png')

