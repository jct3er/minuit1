import numpy as np
from lmfit import minimize, Parameters
import ROOT as r
import matplotlib.pyplot as plt
from scipy import stats


file = r.TFile("distros.root")
hist = file.Get("dist1")

bincount = hist.GetNbinsX()
bin_centers = np.array([hist.GetBinCenter(i) for i in range(1, bincount+1)])
counts = np.array([hist.GetBinContent(i) for i in range(1, bincount+1)])
errs = np.sqrt(counts)
#errs = np.where(errs<1e-10, 1e-10, errs)

def gaus(x, A1, A2, mu1, mu2, sig1, sig2):
    return A1*np.exp(-0.5*(x-mu1)**2/sig1**2) + A2*np.exp(-0.5*(x-mu2)**2/sig2**2)

def gumbel(x, A, mu, beta, c):
    return A*np.exp((x-mu)/beta-np.exp((x-mu)/beta))+c


def objectivegaus(params, x, data, errs):
    A1 = params['A1']
    A2 = params['A2']
    mu1 = params['mu1']
    mu2 = params['mu2']
    sig1 = params['sig1']
    sig2 = params['sig2']
    
    model = gaus(x, A1, A2, mu1, mu2, sig1, sig2)
    res = []
    for i in range(len(errs)):
        if errs[i] == 0:
            continue
        res.append((model[i] - data[i]) / errs[i])
    
    return res

def objectivegum(params, x, data, errs):
    A = params['A']
    mu = params['mu']
    beta = params['beta']
    c = params['c']
    
    model = gumbel(x, A, mu, beta, c)

    res = []
    for i in range(len(errs)):
        if errs[i] == 0:
            continue
        res.append((model[i] - data[i]) / errs[i])
    
    return res


params = Parameters()
params.add("A1", value=800)
params.add("A2", value=800)
params.add("mu1", value=80)
params.add("mu2", value=90)
params.add("sig1", value=5)
params.add("sig2", value=10)
params.add("A", value=1000)
params.add("mu", value=78)
params.add("beta", value=-1)
params.add("c", value=0)


resultgaus = minimize(objectivegaus, params, args=(bin_centers, counts, errs))
resultgum = minimize(objectivegum, params, args=(bin_centers, counts, errs))

gaus_fit = gaus(bin_centers, resultgaus.params['A1'].value, resultgaus.params['A2'].value, resultgaus.params['mu1'].value, resultgaus.params['mu2'].value, resultgaus.params['sig1'].value, resultgaus.params['sig2'].value)

gum_fit = gumbel(bin_centers, resultgum.params['A'].value, resultgum.params['mu'].value, resultgum.params['beta'].value, resultgum.params['c'].value)

gauschi = resultgaus.chisqr
gumchi = resultgum.chisqr

gausndof = resultgaus.nfree
gumndof = resultgum.nfree

print(f"Chi^2 for double gaussian is {resultgaus.chisqr}")
print(f"Chi^2 for gumbel is {resultgum.chisqr}")

gaus_p = 1 - stats.chi2.cdf(gauschi, gausndof)
gum_p = 1 - stats.chi2.cdf(gumchi, gumndof)

print(f"The P-Value for the double gaussian is {gaus_p}")
print(f"The P-Value for the gumbel is {gum_p}")

fig, ax = plt.subplots(2,1, figsize=[9,6])

ax[0].errorbar(bin_centers, counts, yerr=errs, fmt='o', label='Data', markersize=4)
ax[0].plot(bin_centers, gaus_fit, 'r-', label='Fit Gaus')
ax[0].set_xlabel("Value")
ax[0].set_ylabel("Counts")
ax[0].legend()

ax[1].errorbar(bin_centers, counts, yerr=errs, fmt='o', label='Data', markersize=4)
ax[1].plot(bin_centers, gum_fit, 'r-', label='Fit Gumbel')
ax[1].set_xlabel("Value")
ax[1].set_ylabel("Counts")
ax[1].legend()

plt.show()

