#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_excel('Supplemental_Table_S2.xls')
df


# In[2]:


grouped = df.groupby('Antibody')
grouped.mean()


# In[3]:


# ２つ（以上）の属性でgroup化する場合
grouped = df.groupby(['Antibody', 'Cell'])
grouped.mean()


# In[4]:


grouped = df.groupby('Antibody')
grouped["SSP-NSC"].describe()


# ## 作図の例 (SSP-NSC)
# ### 棒グラフ (bar plot)

# In[5]:


grouped = df.groupby('Antibody')
group_nsc = grouped["SSP-NSC"]
group_nsc.mean().plot.bar()


# In[6]:


# 横棒グラフ
group_nsc.mean().plot(kind='barh')


# In[7]:


# y軸を逆順に
ax = group_nsc.mean().plot(kind='barh')
ax.invert_yaxis()


# In[8]:


# エラーバーを表示
sem = group_nsc.std()/np.sqrt(group_nsc.count())
ax = group_nsc.mean().plot(kind='barh', xerr=sem)
ax.invert_yaxis()


# ### 箱ひげ図 (box plot)
# 
# 四分位点、最大値、最小値、外れ値を表示する。

# In[9]:


# (seaborn利用)
import seaborn as sns
sns.boxplot(data=df, x='SSP-NSC', y='Antibody')


# In[10]:


# scaleをlogに
f, ax = plt.subplots(figsize=(6, 4))
ax.set_xscale("log")
sns.boxplot(data=df, x='SSP-NSC', y='Antibody')


# ### stripplot 
# 
# データ点を直接表示。N数が少ない場合に推奨される。

# In[11]:


import seaborn as sns
ax = plt.subplots(figsize=(6, 4))
ax = sns.stripplot(data=df, x='SSP-NSC', y='Antibody')
ax.set_xscale("log")


# # 相関係数ヒートマップ

# In[12]:


grouped = df.groupby('Antibody')
grouped.corr()


# In[13]:


grouped["SSP-NSC", "DROMPA peaks", "FRiP (MACS2)"].corr()


# In[14]:


import seaborn as sns

sns.heatmap(grouped["SSP-NSC", "DROMPA peaks", "FRiP (MACS2)"].corr())


# ### groupbyを使うとわかりづらいので、H3K27acだけを使って描いてみる

# In[15]:


df_27ac = df[df['Antibody'] == 'H3K27ac']                                                                                                                                                                                                
df_27ac


# In[16]:


sns.heatmap(df_27ac.corr())


# In[17]:


# color paletteを変更
cmap = sns.color_palette("coolwarm", 200)
sns.heatmap(df_27ac.corr(), cmap=cmap)


# In[18]:


sns.clustermap(df_27ac.corr(), method='ward', cmap=cmap)


# In[19]:


# お題１：下の図を作成するコードを考える。


# In[20]:


# お題２；お題１で作った図について、エラーバーを標準誤差（SEM）から信頼区間（CI）に変更する。


# In[21]:


# お題３：logスケールboxplotで外れ値を非表示にするオプションを調べる。
# 問：外れ値を表示した方がよい場合と、非表示にした方がよい場合はそれぞれどのようなケースか。


# In[22]:


# お題４：logスケールboxplotをviolin plotに変更する。
# 問：boxplotが望ましい場合とviolin plotが望ましい場合の例を考える。


# In[23]:


# お題５：SSP_Fig2Bを自分で作成してみる。


# In[24]:


# お題６：H3K36me3, Input についてそれぞれ相関係数ヒートマップを描き、Antibodyごとの傾向の違いを観察する。

