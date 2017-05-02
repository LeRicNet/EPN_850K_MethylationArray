import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from scipy.cluster.hierarchy import linkage, dendrogram


samples = pd.read_csv('/Users/ericprince/Desktop/EPN_850K_MethylationArray/msetsq.csv', header=0, index_col=0)
samples = pd.DataFrame.transpose(samples) # we want rows to be cell lines
sample_ids = list(samples.index)

normalized_samples = normalize(samples)

mergings = linkage(normalized_samples, method='complete')

dendrogram(mergings, labels=sample_ids, leaf_rotation=90, leaf_font_size=6)
plt.show()
