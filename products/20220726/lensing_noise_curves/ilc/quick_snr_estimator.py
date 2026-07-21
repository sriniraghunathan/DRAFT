import pickle, gzip, numpy as np, glob, sys, os
from pylab import *
fsval = 14

f1 = 's4wide_lmin100_lmax5000.npy'
f2 = 's4wide_lmin100_lmax5000_lmaxtt3000.npy'

d1 = np.load(f1, allow_pickle = 1, encoding='latin1').item()
d2 = np.load(f2, allow_pickle = 1, encoding='latin1').item()

els, cl_kk, nl_tt_1 = d1['els'], d1['cl_kk'], d1['Nl_TT'].real
nl_tt_2 = d2['Nl_TT'].real

snr_1 = cl_kk / nl_tt_1
snr_2 = cl_kk / nl_tt_2

snr_1[np.isnan(snr_1)] = 0.
snr_1[np.isinf(snr_1)] = 0.
snr_2[np.isnan(snr_2)] = 0.
snr_2[np.isinf(snr_2)] = 0.

clf()
ax = subplot(111, yscale = 'log')
plot(els, cl_kk, lw = 2., color = 'gray', alpha = 0.5)
plot(els, nl_tt_1, lw = 1., color = 'darkgreen', label = r'$\ell_{\rm max}^{T} = 5000$')
plot(els, nl_tt_2, lw = 1., color = 'red', label = r'$\ell_{\rm max}^{T} = 3000$')
legend(loc = 1)

xmin, xmax = 0., 5000.
ymin, ymax = 1e-9, 1e-5
ylim(ymin, ymax); xlim(xmin, xmax)
ylabel(r'$L(L+1)]^{2}$C$_{L}^{\phi \phi}/2\pi$', fontsize = fsval)
xlabel(r'Multipole $L$', fontsize = fsval)
title(r'S4-Wide', fontsize = fsval + 4)
show()
