# Converts p-values to q-values

This Python package implements qvalue and pi0 estimation from Storey and Tibshirani (2003).

The only function that needs to be called is 'qvalue()'. 
It will accept a numpy array of pvalues and will return a numpy array of qvalues.
For instance, given some uniform p-values:

```Python
import numpy as np 
from qvalue import qvalue, pi0est, plot_qvalue

pv = np.random.uniform(0.0, 1.0, size = (1000,))
```

it's possible to convert them in q-values by calling

```Python
qv, pi0 = qvalue(pv)
```

qvalue() also includes a memory-efficient (but not CPU-efficient) procedure for the case
when it's not possible to store both the p-values and the q-values in memory at the same time.
In addition it's possible to change the number of hypotheses tested (m) or the proportion of
null hypotheses (pi0).

pi0() estimates the proportion of null hypotheses:

```Python
pi0 = pi0est(pv)
```

Visualize the q-value and p-values with:

```Python
plot_qvalue(qv,pv)
```

## Installation

Install this repo:
```
pip install --upgrade https://github.com/ncchung/qvalue-python/tarball/master
```

## References

Storey, J. D., & Tibshirani, R. (2003). Statistical significance for genomewide studies. Proceedings of the National Academy of Sciences, 100(16), 9440-9445.