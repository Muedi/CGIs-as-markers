# PCA for CGIs
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_islands_Rnbeads(index_file, meth_site_file):
    """
    opens the betas from a dataset which betas were extracted
    with the RnBeads package for R
    the files are in the same order, so that the islands can be taken as index for the meth values
    """
    island_anno = pd.read_csv(index_file)
    beta_means = pd.read_csv(meth_site_file)
    beta_means["UCSC_CpG_Islands_Name"] = (island_anno["Chromosome"].astype("str") + ":" + island_anno["Start"].astype("str") + "-" + island_anno["End"].astype("str")).values
    beta_means = beta_means[beta_means != 'NA'].dropna()
    beta_means.set_index(['UCSC_CpG_Islands_Name'], inplace = True)
    #    return betas
    return beta_means.astype("float")


# Select if Sites or islands should be used.
sites = False
if sites:
# loads the normalized data from RNBeads sites
    meth_site_file = "Output_RnBeads/complete_HCC_sites_meth.csv"
    data_sites = pd.read_csv(meth_site_file)
    data_sites.dropna(inplace=True)
    data_sites = data_sites.T
    data_sites["Class"] = ["HCC" if ("HCC" in x or "carcinoma" in x) else ("blood") for x in data_sites.index]
    data_sites["Class_code"] = [0 if x == "blood" else 1 for x in data_sites["Class"]]

    data = data_sites

else:
    # loads the normalized data from RNBeads  CpG-islands
    index_file = "Output_RnBeads/complete_HCC_islands_anno.csv"
    meth_file = "Output_RnBeads/complete_HCC_islands_meth.csv"
    data = get_islands_Rnbeads(index_file, meth_file)
    data = data.T
    data["Class"] = ["HCC" if ("HCC" in x or "carcinoma" in x) else ("blood") for x in data.index]
    data["Class_code"] = [0 if x == "blood" else 1 for x in data["Class"]]


# data.drop(["PBMC_25", "PBMC_26"], inplace=True)

"""
                        Principal Component Analysis complete
"""
data["Class"] = ["blood" if  ("GSM" in x or "blood" in x) else "cf_DNA" if ("cf_" in x) else "HCC" if ("HCC_" in x or "carcinoma" in x) else "PBMC" if ("PBMC" in x) else "cfHCC" for x in data.index ]
data["Class_code"] = [0 if x == "blood" else 1 if x == "PBMC" else 2 if x == "cf_DNA" else 3 if x == "HCC" else 4 for x in data["Class"]]
# data = data[[True if "PBMC" in x or "GSM" in x or "Blood" in x else False for x in data.index]]
pca = PCA(n_components=3)
pca_result = pca.fit_transform(data.drop(['Class', 'Class_code'], axis=1).values)

data['pca-one'] = pca_result[:,0]
data['pca-two'] = pca_result[:,1]
data['pca-three'] = pca_result[:,2]
# prints Variation per Principal component
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

# "D-Plot of the first and second PC"
color_number = len(set(data["Class"]))
sns.color_palette("dark")
pca_plot = plt.figure(figsize=(16,10))
pca_plot = sns.scatterplot(
    x="pca-one",
    y="pca-two",
    hue="Class",
    palette=sns.color_palette("dark",color_number),
    data=data,
    legend="full",
    alpha=0.8
)
fig = pca_plot.get_figure()
# fig.savefig("X.png")
for line in range(0, data.shape[0]):
    pca_plot.text(data["pca-one"][line]+0.2, data["pca-two"][line], data.index[line], horizontalalignment="left", size="medium", color="black")
fig = pca_plot.get_figure()
# fig.savefig("X_labeled.png")
# 3D-Plot of the PCA
ax = plt.figure(figsize=(16,10)).gca(projection='3d')
ax.scatter(
    xs=data["pca-one"],
    ys=data["pca-two"],
    zs=data["pca-three"],
    c=data["Class_code"],
    cmap='PiYG'
)
ax.set_xlabel('pca-two')
ax.set_ylabel('pca-one')
ax.set_zlabel('pca-three')
fig = ax.get_figure()
#fig.savefig("X_3D.png")
plt.show()
