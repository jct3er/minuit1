import numpy as np
from lmfit import minimize, Parameters
import ROOT as r
import matplotlib.pyplot as plt


file = r.TFile("fitInputs.root")
hist = file.Get("hdata")
hbkg = file.Get("hbkg")

xbins = hist.GetNbinsX()
ybins = hist.GetNbinsY()
x_centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, xbins+1)])
y_centers = np.array([hist.GetYaxis().GetBinCenter(i) for i in range(1, ybins+1)])
counts = []
for i in x_centers:
    temp = []
    for j in y_centers:
        temp.append(hist.GetBinContent(hist.FindBin(i,j)))
    counts.append(temp)

counts = np.array(counts)
errs = np.sqrt(counts)

def fit(x, y, A, mu1, mu2, sigma1, sigma2, norm):
    val = []
    for i in x:
        temp = []
        for j in y:
            bin_num = hbkg.FindBin(i, j)
            temp.append(A*np.exp(-(i-mu1)**2/sigma1**2)*np.exp(-(j-mu2)**2/sigma2**2)+hbkg.GetBinContent(bin_num)*norm)
        val.append(temp)

    return val

def objective(params, x, y, data, errs):
    A = params['A']
    mu1 = params['mu1']
    mu2 = params['mu2']
    sig1 = params['sigma1']
    sig2 = params['sigma2']
    norm = params['norm']

    model = fit(x, y, A, mu1, mu2, sig1, sig2, norm)

    res = []
    for i in range(len(errs)):
        for j in range(len(errs[0])):
            if errs[i][j] == 0:
                continue
            res.append((model[i][j] - data[i][j]) / errs[i][j])
    
    return res

params = Parameters()
params.add("A", value=90)
params.add("mu1", value=3)
params.add("mu2", value=2)
params.add("sigma1", value=1)
params.add("sigma2", value=1)
params.add("norm", value=1)

results = minimize(objective, params, args=(x_centers, y_centers, counts, errs))

best_fit = fit(x_centers, y_centers, results.params['A'].value, results.params['mu1'].value, results.params['mu2'].value, results.params['sigma1'].value, results.params['sigma2'].value, results.params['norm'].value)

fit_hist = r.TH2F("fit", "fit", 60,0,6,60,0,6)
res_hist = r.TH2F("res", "res", 60,0,6,60,0,6)
signal_hist = r.TH2F("sig", "sig", 60,0,6,60,0,6)

for i in range(xbins):
    for j in range(ybins):
        val = best_fit[i][j]
        fit_hist.SetBinContent(i+1, j+1, val)

        true_val = counts[i][j]
        res = val - true_val
        res_hist.SetBinContent(i+1, j+1, res)

        bkg_val = hbkg.GetBinContent(hist.FindBin(x_centers[i], y_centers[j]))
        sig = val - bkg_val*results.params['norm'].value
        signal_hist.SetBinContent(i+1, j+1, sig)
        

tc = r.TCanvas()
tc.Divide(2,2)
tc.cd(1)
hist.Draw("Lego")
tc.cd(2)
fit_hist.Draw("Lego")
tc.cd(3)
res_hist.Draw("Lego")
tc.cd(4)
signal_hist.Draw("Lego")

tc.Update()
tc.Draw()
r.gApplication.Run()



    
