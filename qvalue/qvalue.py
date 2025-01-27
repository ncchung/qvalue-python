import scipy as sp
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt

def qvalue(pv, pi0=None, m=None, lowmem=False, verbose=False, plot=False):
    """
    Estimates q-values from p-values

    Args
    =====

    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    lowmem: use memory-efficient in-place algorithm
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1

    """
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        pi0 = pi0est(pv=pv, m=m, verbose=verbose)

    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = np.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -np.inf
        prev_qv = last_pv
        for i in range(int(len(pv))-2, -1, -1):
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -np.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]
    else:
        p_ordered = np.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m/len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)

        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])

        # reorder qvalues
        qv_temp = qv.copy()
        qv = np.zeros_like(qv)
        qv[p_ordered] = qv_temp

    if plot:
        plot_qvalue(qv, pv)
        
    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv, pi0

def pi0est(pv, m=None, verbose=False):
    """
    Estimates pi0 from p-values

    Args
    =====

    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    """
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0

    # evaluate pi0 for different lambdas
    pi0 = []
    lam = np.arange(0, 0.90, 0.01)
    counts = np.array([(pv > i).sum() for i in np.arange(0, 0.9, 0.01)])
    for l in range(len(lam)):
        pi0.append(counts[l]/(m*(1-lam[l])))

    pi0 = np.array(pi0)

    # fit natural cubic spline
    tck = interpolate.splrep(lam, pi0, k=3)
    pi0 = interpolate.splev(lam[-1], tck)
    if verbose:
        print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)

    if pi0 > 1:
        if verbose:
            print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
        pi0 = 1.0

    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0

    return pi0

def plot_qvalue(qv, pv):
    """
    Visualize the q-value results

    Args
    =====
    qv: q-values
    pv: p-values
    """
    qv_thresholds = np.linspace(0, 1, 100)
    num_significant_tests = [np.sum(qv <= threshold) for threshold in qv_thresholds]
    
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    id = np.argsort(pv)
    pv = pv[id]
    qv = qv[id]

    # pv vs. qv
    axes[0].plot(pv, qv)
    axes[0].set_title('p-values vs. q-values')
    axes[0].set_xlabel('p-values')
    axes[0].set_ylabel('q-values')

    # qv vs. significant tests
    axes[1].plot(qv_thresholds, num_significant_tests)
    axes[1].set_title('q-value vs. # Significant Tests')
    axes[1].set_xlabel('q-value')
    axes[1].set_ylabel('Number of Significant Tests')

    plt.tight_layout()
    plt.show()