# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:06:02 2020

@author: ajn16
"""

# Identify C45+ and ASMA - cells

data = pd.DataFrame(adata.X, columns = adata.var.index, index= adata.obs.index)
d= data
data = np.log1p(d)

geneA = 'ASMA'
geneB = 'CD45'

# CD45
CD45 = data[geneA]
bottom = np.percentile(CD45, 75)
top = np.percentile(CD45, 99)

# ASMA
ASMA = data[geneB]
middle = np.percentile(ASMA, 90)

# Subset cells
positive = data[data[geneA].between(bottom, top)]
len(positive)

# Expression of gene B in postive cells
geneA= 'CD20'
cutoff = np.mean(positive[geneA])
cutoff = np.median(positive[geneA])
cutoff = np.percentile(positive[geneA], 25)

negative = data[data[geneB] > middle]
len(negative)

real = [x for x in positive.index if x not in negative.index]
real_p = d[d.index.isin(real)]
real_n = d[~d.index.isin(real)]
#real_n = d[d.index.isin(negative.index)]
len(real)

r_p = real_p.sample(n = 50000) 
r_n = real_n.sample(n = 50000) 

len(r_p)

r_p['label'] = 1
r_n['label'] = 0

data_model = r_p.append(r_n)
data_model = data_model[[geneA,geneB,'label']]
f = full_data[[geneA,geneB]]





for i in adata.var.index:
    try:
        del adata.obs[i]      
    except:
        pass

adata.obs[geneA] = prediction
adata.obs[geneA] = adata.obs[geneA].astype('category')

expression_prediction_plot (adata, x='X_position', y='Y_position', color=geneA, 
                            save_figure= False, plot_all=False, out_dir=out_dir, 
                            alpha= 0.3, figure_type="png", s=1, cmap='gist_heat', 
                            figsize=(12, 5))

geneA = 'CD20'
import random
exp = np.log1p(adata[:,geneA].X)
exp = random.sample(list(exp), 500)
sns.distplot(exp, rug=False, hist=True)

cutoff = 7

a = np.log1p(adata[:,geneA].X)
a = np.where(a > cutoff, 1, 0)
adata.obs[geneA] = a
adata.obs[geneA] = adata.obs[geneA].astype('category')
expression_prediction_plot (adata, x='X_position', y='Y_position', color=geneA, 
                            save_figure= False, plot_all=False, out_dir=None, 
                            alpha= 1, figure_type="png", s=0.1, cmap='gist_heat', 
                            figsize=(12, 5))



data_gm = np.array(exp).reshape(-1, 1)
gmm = GaussianMixture(n_components=2)
gmm.fit(data_gm)

gmm.means_

high = list()
low = list()

for i in range(100):
    exp = random.sample(list(exp), 500)
    data_gm = np.array(exp).reshape(-1, 1)
    gmm = GaussianMixture(n_components=2)
    gmm.fit(data_gm) 
    high.append(np.max(gmm.means_))
    low.append(np.min(gmm.means_))
    
np.mean(high)

full = np.log1p(adata[:,geneA].X)
sns.distplot(full, rug=False, hist=True)

data_gm = np.array(full).reshape(-1, 1)
gmm = GaussianMixture(n_components=3)
gmm.fit(data_gm)
gmm.means_   

np.exp(7)

# Scale data between 0 and 10
fu_d = adata[:,geneA].X
full_scale = np.interp(full, (full.min(), full.max()), (0, 10))

sns.distplot(data['CD45'], rug=False, hist=True)

data_gm = pd.DataFrame(data['CD45'])
gmm = GaussianMixture(n_components=2)
gmm.fit(data_gm)
gmm.means_ 


from sklearn.mixture import GaussianMixture 

cutoff = 0
# Get the standard deviation
a = np.log1p(adata[:,geneA].X)
a = np.where(a > cutoff, 1, 0)
adata.obs[geneA] = a
adata.obs[geneA] = adata.obs[geneA].astype('category')
expression_prediction_plot (adata, x='X_position', y='Y_position', color=geneA, 
                            save_figure= False, plot_all=False, out_dir=None, 
                            alpha= 1, figure_type="png", s=0.1, cmap='gist_heat', 
                            figsize=(12, 5))



data = d
sum_data = data.sum(axis=1) # Calculate total count for each cell
n_count = data.div(sum_data, axis=0) # Divide genes by total count for every cell
med = np.median(list(itertools.chain(*data.values.tolist()))) # Calculate median count of the entire dataset
n_count = n_count*med # Multiply by scaling fator (median count of entire dataset)
n_log = np.log1p(n_count)

scaled= n_log.apply(np.sqrt)

geneA = 'CD20'
scaled = n_log[geneA]
scaled = np.interp(scaled, (scaled.min(), scaled.max()), (-1, 1))
sns.distplot(scaled, rug=False, hist=True)


scaler = MinMaxScaler(feature_range=(-1, 1))
s = scaler.fit_transform(n_log)
s = pd.DataFrame(s)
s.columns = n_log.columns

np.array_equal(scaled, np.array(s[geneA]))


full = np.log1p(adata[:,geneA].X)
f = full - cutoff
f[f < 0] = 0
f = full/cutoff
sns.distplot(data['CD45'], rug=False, hist=True, label = 'A')
sns.distplot(data['CD30'], rug=False, hist=True, label = 'A')


full = adata[:,geneA].X
c = np.exp(cutoff)
f = full - c
f = f.tolist()
for i, x in enumerate(f):
    if x <= 0: f[i] = 0
sns.distplot(np.log1p(f), rug=False, hist=True)

data_gm = scaled.reshape(-1, 1)
gmm = GaussianMixture(n_components=2)
gmm.fit(data_gm)
gmm.means_ 

geneA = 'CD45'
geneB = 'ASMA'
geneC = 'CD3D'

sns.distplot(s[geneA], rug=False, hist=False, label = geneA)
sns.distplot(s[geneB], rug=False, hist=False, label = geneB)
sns.distplot(s[geneC], rug=False, hist=False, label = geneC)

geneA = 'CD3D'
cutoff = 0.5
# Get the standard deviation
a = data[geneA]
a = np.where(a > cutoff, 1, 0)
adata.obs[geneA] = a
adata.obs[geneA] = adata.obs[geneA].astype('category')
expression_prediction_plot (adata, x='X_position', y='Y_position', color=geneA, 
                            save_figure= False, plot_all=False, out_dir=None, 
                            alpha= 1, figure_type="png", s=0.1, cmap='gist_heat', 
                            figsize=(12, 5))


scaler = MinMaxScaler(feature_range=(-1, 1))
a = scaler.fit_transform(score.reshape(-1, 1))


geneA = 'CD20'
cutoff = 0.5
a = all_scaled_data[geneA]
a = np.where(a > cutoff, 1, 0)
c_col = ['red' if x==1 else 'black' for x in a]
fig, ax = plt.subplots()
plt.scatter(x=adata.obs['X_position'], 
            y=adata.obs['Y_position'], 
            c=c_col, edgecolors='none', s=1, alpha=0.8)
ax.invert_yaxis()


adata.obs['label'] = labels
adata.obs['label'] = adata.obs['label'].astype('category')

cluster_plot(adata, color = 'label', x='X_position', y='Y_position', 
             cluster_name = 'mature b cells', s=0.001, cmaps='gist_heat')

adata.obs['label'][adata.obs['label'] == 'macrophages']


# Sequential gating
cutoff = 0.5
B_cells = data[ (data['CD45'] > cutoff) & (data['CD20'] > cutoff)]
T_cells = data[ (data['CD45'] > cutoff) & (data['CD3D'] > cutoff)]
B_T_overlap = data[ (data['CD45'] > cutoff) & (data['CD20'] > cutoff) & (data['CD3D'] > cutoff)]

len(B_T_overlap)
# Plot
a = np.where(data.index.isin(B_T_overlap.index), 1, 0)
c_col = ['red' if x==1 else 'black' for x in a]
sns.set_style("white")
fig, ax = plt.subplots()
plt.scatter(x=adata.obs['X_position'], 
            y=adata.obs['Y_position'], 
            c=c_col, edgecolors='none', s=1, alpha=0.8)
ax.invert_yaxis()


sns.set(rc={'figure.figsize':(11.7,8.27)})
# cluster plot not working
a = np.where(adata.obs['label'] == 'nk cells', 1, 0)
c_col = ['red' if x==1 else 'black' for x in a]
sns.set_style("white")
fig, ax = plt.subplots()
plt.scatter(x=adata.obs['X_position'], 
            y=adata.obs['Y_position'], 
            c=c_col, edgecolors='none', s=2, alpha=0.8)
ax.invert_yaxis()




# Playing with gasussina mixture
from sklearn.mixture import GaussianMixture 
geneA = 'FOXP3'
full = data[geneA]
sns.distplot(full, rug=False, hist=True)

data_gm = np.array(full).reshape(-1, 1)
gmm = GaussianMixture(n_components=2, means_init=[[0.25],[0.75]], covariance_type='tied')
gmm.fit(data_gm)
gmm.means_


n_log = np.log1p(data) # Log transform data
scaler = MinMaxScaler(feature_range=(0, 1))
s = scaler.fit_transform(n_log)
s= pd.DataFrame(s, columns = adata.var.index, index= adata.obs.index)

geneA = 'CD45'
sns.distplot(log_data[geneA], rug=False, hist=True)
sns.distplot(x[geneA], rug=False, hist=True)


a = data['FOXP3'].values
b = data['CD11B'].values


scaler_high = MinMaxScaler(feature_range=(0.5, 1))
scaler_low = MinMaxScaler(feature_range=(0, 0.5))

c = scaler_high.fit_transform(a[a >= 0.3].reshape(1, -1))
c = scaler_low.fit_transform(a[a <= 0.3].reshape(1, -1))

x = scaler.fit_transform(n_log)
x = pd.DataFrame(x, columns = adata.var.index, index= adata.obs.index)


geneA = 'CD20'
cutoff = 0.5
a = adata[adata.obs['ImageId']==629148]
b = a[:,geneA].X
b = np.where(b > cutoff, 1, 0)
c_col = ['red' if x==1 else 'black' for x in b]
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_style("white")
fig, ax = plt.subplots()
plt.scatter(x=a.obs['X_position'], 
            y=a.obs['Y_position'], 
            c=c_col, edgecolors='none', s=1, alpha=0.8)
ax.invert_yaxis()

len(b[b==1]) #629148  326958


geneA = 'CD45'
sns.distplot(np.log1p(adata.raw[:,geneA].X), rug=False, hist=True)

cutoff = 6.8
a = np.log1p(adata.raw[:,geneA].X)
a = np.where(a > cutoff, 1, 0)
c_col = ['red' if x==1 else 'black' for x in a]
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_style("white")
fig, ax = plt.subplots()
plt.scatter(x=adata.obs['X_position'], 
            y=adata.obs['Y_position'], 
            c=c_col, edgecolors='none', s=1, alpha=0.8)
ax.invert_yaxis()


import matplotlib.colors as colors
import matplotlib.cm as cmx

c = data['phenotype']
uniq = list(set(c))
    
# Set the color map to match the number of species
z = range(1,len(uniq))
hot = plt.get_cmap('viridis')
cNorm  = colors.Normalize(vmin=0, vmax=len(uniq))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=hot)

# Plot each species
fig, ax = plt.subplots()
for i in range(len(uniq)):
    indx = c == uniq[i]
    pts = ax.scatter(x[indx], y[indx], s=0.01, color=scalarMap.to_rgba(i), label=uniq[i])
    ax.invert_yaxis()
    
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.show()


# covarianve
cov = gmm.covariances_
sd = [ np.sqrt(  np.trace(cov[i])/1) for i in range(0,1) ]


sns.distplot(np.log1p(adata.raw[:,geneA].X), rug=False, hist=True)
sns.distplot(adata[:,geneA].X, rug=False, hist=True)


data = pd.DataFrame(s, columns = adata.var.index, index= adata.obs.index)
sns.distplot(data[geneA], rug=False, hist=True)


bdata.obs['phenotype'] = bdata.obs['phenotype'].astype('category')


bdata = sc.pp.subsample(adata, fraction=0.9, copy=True)
sc.tl.dendrogram(adata, groupby='parc')
sc.pl.dendrogram(bdata, groupby='phenotype')
markers = list(adata.var.index)
sc.pl.matrixplot(adata, markers, groupby='kmeans_renamed', dendrogram=True, use_raw=True, cmap=c, standard_scale='var', log=True, swap_axes=True)
#sc.pl.matrixplot(adata, markers, groupby='kmeans', dendrogram=True, use_raw=False, cmap=c)
sc.pl.matrixplot(adata, markers, groupby='kmeans_renamed', dendrogram=True, use_raw=False, cmap=c, swap_axes=True)




sc.pl.heatmap(adata, markers, dendrogram=True, groupby='phenotype_renamed', use_raw=True, standard_scale='var', cmap=c, log=True)
sc.pl.heatmap(adata, markers, dendrogram=True, groupby='kmeans_renamed', use_raw=False, cmap=c, swap_axes=True)

markers = ['CD45','CD3D','CD2','CD5','CD25','CD21','ASMA','FOXP3','CD68','CD11B','CD11C','CD163',
           'CD206','CD20','HLADR','KI67','CD8A','PD1','CD4','PSTAT3','PDL1','PS6']

markers = list(adata.var.index)
markers = ["CD2", "CD5", "CD45", "PD1", "KI67", "CD7", "CD8A", "CD4", "PS6", "PSTAT3"]

sc.pl.matrixplot(adata, markers, groupby='kmeans', dendrogram=False, 
                 use_raw=True, standard_scale='var', cmap=c, log=True, swap_axes=True)



sc.tl.rank_genes_groups(adata, 'phenotype_renamed', method='t-test')

sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, fontsize = 40)

sc.pl.stacked_violin(adata, markers, groupby='phenotype_renamed', rotation=90, log=True)

# Check for batch effect
# If not, check if combining two image analysis is possible
# 

tumor.raw.X


adata = sc.read('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL2/PTCL2.h5ad')

bdata = sc.read('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/4-celltype_calling/PTCL1/PTCL1.h5ad')


geneA = 'CD45'
sns.distplot(np.log1p(adata.raw[:,geneA].X), rug=False, hist=True)
sns.distplot(np.log1p(bdata.raw[:,geneA].X), rug=False, hist=True)

pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/celltype_calling/PTCL1.csv')


image_path = ['/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL1_450.csv',
          '/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL2_552.csv']


data = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/raw_data/whole_sections/PTCL1_450.csv')



geneA = 'FOXP3'
sns.distplot(adata[:,geneA].X, rug=False, hist=True)




