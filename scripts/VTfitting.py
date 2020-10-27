#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameters, Minimizer, fit_report, Model
from uncertainties import correlated_values,nominal_value
from uncertainties import unumpy as unp
from time import time
import datetime

def func(x, params):
    a0, a1, a2,a3, a4, a5, a6, a7, a8 = params
    return a0 + np.exp(-(a1/x))*(a2/x**2.0 + a3/x + a4 + a5*x + a6*np.exp(-a7*(-a8+x)**2.0))


def func_model(x, a0, a1, a2,a3, a4, a5, a6, a7, a8):
    return a0 + np.exp(-(a1/x))*(a2/x**2.0 + a3/x + a4 + a5*x + a6*np.exp(-a7*(-a8+x)**2.0))

def func_unc(x, a0, a1, a2,a3, a4, a5, a6, a7, a8):
    return a0 + unp.exp(-(a1/x))*(a2/x**2.0 + a3/x + a4 + a5*x + a6*unp.exp(-a7*(-a8+x)**2.0))

def CTE(x, params):
    a0, a1, a2, a3,a4, a5, a6, a7, a8 = params
    dVdT = (a1*np.exp(-a1/x)*(a2/x**2+a3/x+a4+a5*x+a6*np.exp(-a7*(x-a8)**2)))/x**2+np.exp(-a1/x)*(-(2*a2)/x**3-a3/x**2+a5-2*a6*a7*(x-a8)*np.exp(-a7*(x-a8)**2))
    CTE = dVdT/func(x,params)
    return CTE

def CTE_model(x, a0, a1, a2, a3, a4, a5, a6, a7, a8):
    dVdT = (a1*np.exp(-a1/x)*(a2/x**2+a3/x+a4+a5*x+a6*np.exp(-a7*(x-a8)**2)))/x**2+np.exp(-a1/x)*(-(2*a2)/x**3-a3/x**2+a5-2*a6*a7*(x-a8)*np.exp(-a7*(x-a8)**2))
    CTE = dVdT/func_model(x, a0, a1, a2,a3, a4, a5, a6, a7, a8)
    return CTE

def CTE_unc(x, a0, a1, a2, a3, a4, a5, a6, a7, a8):
    V = a0 + unp.exp(-(a1/x))*(a2/x**2.0 + a3/x + a4 + a5*x + a6*unp.exp(-a7*(-a8+x)**2.0))
    dVdT = (a1*unp.exp(-a1/x)*(a2/x**2+a3/x+a4+a5*x+a6*unp.exp(-a7*(x-a8)**2)))/x**2+unp.exp(-a1/x)*(-(2*a2)/x**3-a3/x**2+a5-2*a6*a7*(x-a8)*unp.exp(-a7*(x-a8)**2))
    CTE = dVdT/V
    return CTE

def func2min(params, T, V):
    parametros = [params[k] for k in params.keys()]
    model = func(T, parametros)
    return model - V

def fit_volume_vs_temperature(T, V):

    params = Parameters()
    params.add_many(('a0', 400,  True,       0,   8000,   None, None),
                    ('a1',   0,  True,     0,   1000,   None, None),
                    ('a2',   0,  True,   -15000,   15000,   None, None),
                    ('a3',   0,  True,   -150,   150,   None, None),
                    ('a4',   0,  True,   -2,   2,   None, None),
                    ('a5',   0.001,  True,   -0.01,   0.01,   None, None),
                    ('a6',   0.0,  True,   -5,   0.0,   None, None),
                    ('a7',   0.001,  True,   0.00001,   0.001,   None, None),
                    ('a8',   0.001,  True,   0,   300,   None, None)
                    )

    minimizer = Minimizer(func2min, params, fcn_args=(T, V))
    kws={'popsize':500,  'mutation':(0.8,1.2), 'recombination':0.8, 'updating':'deferred', 'workers':-1}
    #out = minimizer.minimize(method='differential_evolution', max_nfev=50000000,**kws)
    out = minimizer.minimize(method='ampgo')
    fit = func2min(out.params, T, V)

    print(fit_report(out), flush=True)
    print('Cost:',np.sum(fit))
    print('---------- DE bounds ----------')
    for k in out.params.keys():
        print(k,out.params[k].min,out.params[k].value,out.params[k].max)

    return out

def second_fit(temperature,volume,out):
    modelo = Model(func_model)
    result = modelo.fit(volume,out.params,method='nelder',x=temperature)
    print(fit_report(result), flush=True)
    return result


def plot(x_data, y_data, x_fit, y_fit,y_cte):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,dpi=120,figsize=(10,4))
    ax1.plot(x_data, y_data, 'o', label='data')
    ax1.plot(x_fit, y_fit, '-',label='fit')
    ax1.set_ylabel('Volume (A3)')
    ax1.set_xlabel('Temperature (K)')
    ax1.legend()
    ax2.plot(x_fit, y_cte)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('CET (K$^{-6}$)')
    plt.tight_layout()
    plt.savefig('plotVvsT.png')
    #plt.show()

def plot_with_uncertainties(x_data, y_data, x, y, y_unc, z, z_unc):
    #volume versus temeperature
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,dpi=120,figsize=(10,4))
    ax1.plot(x_data, y_data, 'o', label='data')
    #uncertainty band
    y_up = y+y_unc
    y_down = y-y_unc
    ax1.fill_between(x, y_up, y_down, facecolor='#888888', alpha=0.25, label = '1$\sigma$', edgecolor=None)
    ax1.plot(x, y, '-',label='fit')
    ax1.set_ylabel('Volume (A3)')
    ax1.set_xlabel('Temperature (K)')
    ax1.legend()

    #coefficitent of thermal expansion versus temperature
    z_up = z+z_unc
    z_down = z-z_unc
    ax2.fill_between(x, z_up, z_down, facecolor='#888888', alpha=0.25, label = '1$\sigma$', edgecolor=None)
    ax2.plot(x, 1e6*z)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('CET (K$^{-6}$)')
    ax2.legend()

    plt.tight_layout()
    plt.savefig('plotVvsT.png')
    #plt.show()

##################################################################################

def perform_fitting_routine(temperature,volume):
    time_i = time()

    t_fit = np.linspace(temperature[0], temperature[-1],2*len(temperature))

    ResultDE = fit_volume_vs_temperature(temperature,volume)

    if not np.isnan(np.sqrt(np.diag(ResultDE.covar))).any():
        #If uncertainties could be correctly calculated using DE, we use these results
        parametros = ResultDE.params
        incertezas = np.sqrt(np.diag(ResultDE.covar))

        (a0,a1,a2,a3,a4,a5,a6,a7,a8) = correlated_values([ResultDE.params[i].value for i in ResultDE.params.keys()], ResultDE.covar)

        v_complete = func_unc(t_fit,a0,a1,a2,a3,a4,a5,a6,a7,a8)
        v_fit = np.array([x.nominal_value for x in v_complete])
        v_unc = np.array([x.std_dev for x in v_complete])

        cte_complete = CTE_unc(t_fit, a0, a1, a2, a3, a4, a5,a6,a7,a8)
        cte_fit = np.array([x.nominal_value for x in cte_complete])
        cte_unc = np.array([x.std_dev for x in cte_complete])

        plot_with_uncertainties(temperature, volume, t_fit, v_fit, v_unc, cte_fit, cte_unc)

        with open('fitted_values.dat','w') as f:
            f.write('#{temp:<14}\t{vol:<14}\t{v_err:<14}\t{cte:<14}\t{cte_err:<14}\n'.format(
                temp='Temperature(K)',vol='Volume(A3)',v_err='Vol.Std(A3)',cte='CTE(x1E6/K)',cte_err='CTEStd(x1E6/K)'))
            for i in range(len(t_fit)):
                f.write('{temp:<14}\t{vol:<14}\t{v_err:<14}\t{cte:<14}\t{cte_err:<14}\n'.format(
                    temp=t_fit[i],vol=v_fit[i],v_err=v_unc[i],cte=1e6*cte_fit[i],cte_err=1e6*cte_unc[i]))


    else:
        print('Uncertainties could not be correctly determined with DE. Trying Nelder-Mead.')
        model_result = second_fit(temperature,volume,ResultDE)
        #Check if the uncertainties were correctly obtained by DE
        has_uncertainties = not np.isnan(np.sqrt(np.diag(model_result.covar))).any()
        print('Has Uncertainties:',has_uncertainties)
        fitted_params = [model_result.params[k].value for k in model_result.params.keys()]

        if has_uncertainties:
            #If uncertainties could not be calculated using DE, we perform a second fit using Nelder-Mead
            print('Uncertainties could not be correctly determined with DE. Trying Nelder-Mead.')
            model_result = second_fit(temperature,volume,ResultDE)

            parametros = model_result.values

            v_fit = model_result.eval(x=t_fit)
            v_unc = model_result.eval_uncertainty(x=t_fit)

            #calculating CET using correlations
            (a0,a1,a2,a3,a4,a5, a6, a7, a8) = correlated_values([parametros[i] for i in parametros], model_result.covar)
            cte_complete = CTE_unc(t_fit, a0, a1, a2, a3, a4, a5, a6, a7, a8)
            cte_fit = np.array([x.nominal_value for x in cte_complete])
            cte_unc = np.array([x.std_dev for x in cte_complete])

            plot_with_uncertainties(temperature, volume, t_fit, v_fit, v_unc, cte_fit, cte_unc)

            with open('fitted_values.dat','w') as f:
                f.write('#{temp:<14}\t{vol:<14}\t{v_err:<14}\t{cte:<14}\t{cte_err:<14}\n'.format(
                    temp='Temperature(K)',vol='Volume(A3)',v_err='Vol.Std(A3)',cte='CTE(x1E6/K)',cte_err='CTEStd(x1E6/K)'))
                for i in range(len(t_fit)):
                    f.write('{temp:<14}\t{vol:<14}\t{v_err:<14}\t{cte:<14}\t{cte_err:<14}\n'.format(
                        temp=t_fit[i],vol=v_fit[i],v_err=v_unc[i],cte=1e6*cte_fit[i],cte_err=1e6*cte_unc[i]))

        else:
            print('### WARNING: Uncertainties could not be correctly obtianed.')
            v_fit = func(t_fit, fitted_params)
            cte_fit = CTE(t_fit,fitted_params)

            plot(temperature,volume,t_fit, v_fit, cte_fit)

            with open('fitted_values.dat','w') as f:
                f.write('#{temp:<14}\t{vol:<14}\t{v_err:<14}\t{cte:<14}\t{cte_err:<14}\n'.format(
                    temp='Temperature(K)',vol='Volume(A3)',v_err='Vol.Std(A3)',cte='CTE(x1E6/K)',cte_err='CTEStd(x1E6/K)'))
                for i in range(len(t_fit)):
                    f.write('{temp:<14}\t{vol:<14}\t{v_err:<14}\t{cte:<14}\t{cte_err:<14}\n'.format(
                        temp=t_fit[i],vol=v_fit[i],v_err=0.0,cte=1e6*cte_fit[i],cte_err=0.0))

    this_time = time()-time_i
    print('Time of fitting: ',datetime.timedelta(seconds=this_time))
