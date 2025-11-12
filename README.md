# minuit1

- expFit.cpp: example using the more general fitting interface in minuit
- expFit.ipynb: equivalent example using lmfit
- SimultaneousExps(lm).ipynb: generation of histograms with correlated signals for simultaneous fit exercise
- rootExample.cpp: just another example of using ROOT classes in a C++ program
- *.root files: various input histograms for fitting exercises

-----

Team member names and computing IDs: AAA(aaa) BBB(bbb)

-----

Exercise one is in the file distrofit.py

The double gaussian fit has a chi^2 of 258.05 and a P value of 0.0

The gumbel fit has a chi^2 of 19.85 and a P value of 0.897

So the gumbel is a better fit.

Exercise 2 comments:
--

The mean of my gaussian is 74.6+- 1.08
The sigma if my gaussian is 4.596 +- 1.04

The chi^2 of my fit is 76.168, this is a good fit as the fit histograms have 50 bins each so the combined degrees of freedom is 92 so my reduced chi^2 is 0.828. So I should have a good fit without over fitting.

-----

Exercise 3 comments:
--

Exercise 3 is in the file 2dfit.py

The total number of signal events is 34026 +- 201. I calculated this by adding for each bin the hdata count in that bin subtracted by hbkg of that bin times the background normalization. I got the error by doing error propagation to get the error on the signal count in each bin and then added those in quadriture.


-----
