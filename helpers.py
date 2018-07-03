#!/usr/bin/env python3
from matplotlib import pyplot as pl
import seaborn

def GiveInfo(data):
    print("data type: ", data.dtype, ", shape: ", data.shape,
          ", min: ", data.min(), ", max: ", data.max(), ", mean: %f" % data.mean())

def PercOK(data1, data2, rtol=1e-5, atol=1e-8):
    a = np.isclose(data1, data2, rtol=rtol, atol=atol)
    pa = a.sum()/a.size * 100
    print("percOK: %.2f%%" % pa)
    
def PlotDists(data1, data2, lbl1='idl',lbl2='py',dims=(6,4), 
              rng1=None, rng2 = None, title=None, col2='k'):
    f = pl.figure(figsize=dims)
    ax = f.add_subplot(111)
    if not rng1:
        rng1 = (data1.min(), data1.max())
    if not rng2:
        rng2 = (data2.min(), data2.max())
    ax.hist(data1.flatten(), bins=1000, range=rng1, label=lbl1,
            histtype='stepfilled', normed=True)
    ax.hist(data2.flatten(), bins=1000, range=rng2, label=lbl2,
            histtype='stepfilled', color=col2, alpha=0.5, normed=True)
    ax.legend()
    if title:
        ax.set_title(title)
    