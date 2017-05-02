import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


"""
First, we take in the data that was compiled using the limma, minfi, and Illumina packages from bioconductor in R. The
data is transposed to reflect that CpG island methylation statuses are the variables and the cell lines are the
observations of those variables.
"""

samples = pd.read_csv('/Users/ericprince/Desktop/EPN_850K_MethylationArray/msetsq.csv', header=0, index_col=0)


samples_t = pd.DataFrame.transpose(samples)

def inertias_fx(samples, k_range):

    inertias = []

    for k in k_range:
        model = KMeans(k)   # Create a KMeans with k clusters
        model.fit(samples)     # Fit the model to samples.
        inertias.append(model.inertia_)

    return inertias

plt.figure(1)

ks1 = range(1, 5)

plt.subplot(211)
plt.plot(ks1, inertias_fx(samples=samples, k_range=ks1), '-o')
plt.title('Col = CpG Measurement')
plt.xlabel('number of clusters, k')
plt.ylabel('inertia')
plt.xticks(ks1)

ks2 = range(1, 16)  # range can not end higher than 16, because there are only 16 different samples.

plt.subplot(212)
plt.plot(ks2, inertias_fx(samples=samples, k_range=ks2), '-o')
plt.title('Col = Cell Lines')
plt.xlabel('number of clusters, k')
plt.ylabel('inertia')
plt.xticks(ks2)

plt.show()