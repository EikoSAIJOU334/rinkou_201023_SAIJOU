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


# In[106]:


grouped["SSP-NSC"].describe("count")


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


# ## ここから、お題について

# In[19]:


# お題１：下の図を作成するコードを考える。


# - こちらを参考に行った
# 
# https://matplotlib.org/3.3.2/gallery/units/bar_unit_demo.html#sphx-glr-gallery-units-bar-unit-demo-py
# 
# - こちらも参考になりました
# 
# https://stackoverrun.com/ja/q/12788276

# In[105]:


# プロットするグループの数を、N =に入れておく
N = 7
nsc_means = group_nsc.mean() # SSP_NSCのmeanを代入、groupbyで作った値を使用
sem_SSP = group_nsc.std()/np.sqrt(group_nsc.count()) # SSP_NSCのsemを計算、groupbyで作ったstdをrout(n)で割る

fig, ax = plt.subplots() # figureの枠を作る

ind = np.arange(N)    # the x locations for the groups 
# このindの意味？ indのなかみは、初項が0で公差が1の等差数列を要素とするndarray[0, 1, 2, 3, 4, 5, 6]
width = 0.25         # the width of the bars
ax.barh(ind, nsc_means, width, xerr=sem_SSP, label='SSP-NSC') # barhで横棒グラフを描画
# indで指定した場所に、nsc_meansを、widthの幅で、エラーバーはsem_SSPで書く。ラベルは"SSP-NSC"

ppqt_means = group_ppqt.mean() # PPQT_RSCのmeanを代入
sem_PPQT = group_ppqt.std()/np.sqrt(group_ppqt.count()) # PPQT_RSCのsemを計算

ax.barh(ind + width, ppqt_means, width,  xerr=sem_PPQT, 
       label='PPQT-RSC') #ind+widthで、widthの分ずらした場所にプロットする、なぜこれでいけるのか。
# ind + widthで指定した場所に、ppqt_meansを、widthの幅で、エラーバーはsem_PPQTで書く。ラベルは"PPQT_RSC"

ax.set_yticks(ind + width / 2) # y軸の印を付ける場所は、ind + width /2の場所。つまり、ふたつの真ん中。
ax.set_yticklabels(pd.unique(df["Antibody"])) # y軸のtickごとのラベルは、dfのAntibodyのユニークなもの。
ax.set_ylabel('Antibody') # y軸のラベル指定。

ax.legend() # legendを表示
ax.autoscale_view() # 軸のスケールを自動的に設定

plt.show()


# In[113]:


# お題２；お題１で作った図について、エラーバーを標準誤差（SEM）から信頼区間（CI）に変更する。


# - こちらを参考にしました。
# https://obgyn.jp/ci/

# In[111]:


from scipy import stats

# t分布を仮定した場合の信頼区間
ci_SSP = stats.t.interval(alpha = 0.95,          # 信頼区間
                      df = group_nsc.count()-1, # 自由度
                      loc = group_nsc.mean(),            # 標本平均
                      scale = group_nsc.sem())           # 標準誤差

ci_PPQT = stats.t.interval(alpha = 0.95,          # 信頼区間
                      df = group_ppqt.count()-1, # 自由度
                      loc = group_ppqt.mean(),            # 標本平均
                      scale = group_ppqt.sem())           # 標準誤差


# In[112]:


N = 7
nsc_means = group_nsc.mean() # SSP_NSCのmeanを代入、groupbyで作った値を使用

fig, ax = plt.subplots() # figureの枠を作る

ind = np.arange(N)    
width = 0.25        
ax.barh(ind, nsc_means, width, xerr=ci_SSP, label='SSP-NSC') 
# xerrをciに変更

ppqt_means = group_ppqt.mean() # PPQT_RSCのmeanを代入

ax.barh(ind + width, ppqt_means, width,  xerr=ci_PPQT, 
       label='PPQT-RSC') 

ax.set_yticks(ind + width / 2) 
ax.set_yticklabels(pd.unique(df["Antibody"])) 
ax.set_ylabel('Antibody') 

ax.legend() # legendを表示
ax.autoscale_view() # 軸のスケールを自動的に設定

plt.show()


# In[ ]:


# お題３：logスケールboxplotで外れ値を非表示にするオプションを調べる。


# - こちらを参考に行った。
# https://www.haya-programming.com/entry/2019/10/07/202416

# In[114]:


f, ax = plt.subplots(figsize=(6, 4))
ax.set_xscale("log")
sns.boxplot(data=df, x='SSP-NSC', y='Antibody', sym="") # 非表示にするだけで、検出は行っている


# In[115]:


f, ax = plt.subplots(figsize=(6, 4))
ax.set_xscale("log")
sns.boxplot(data=df, x='SSP-NSC', y='Antibody', whis="range") # 検出自体を行わない、ひげが長くなる


# In[ ]:


# 問：外れ値を表示した方がよい場合と、非表示にした方がよい場合はそれぞれどのようなケースか。


# - サンプル数が多く、見た目的に気になるほどの外れ値が出てしまう場合に、見た目を整える意味で非表示にする。
# - サンプル数が少なく、外れ値にも意味がある場合には外れ値を表示する。

# In[ ]:


# お題４：logスケールboxplotをviolin plotに変更する。


# In[116]:


f, ax = plt.subplots(figsize=(6, 4))
ax.set_xscale("log")
sns.violinplot(data=df, x='SSP-NSC', y='Antibody')


# In[ ]:


# 問：boxplotが望ましい場合とviolin plotが望ましい場合の例を考える。


# - 今回のように分布の幅が広く、密度にあまり差がない場合、boxplotが見やすいので望ましい。

# In[ ]:


# お題５：SSP_Fig2Bを自分で作成してみる。


# - こちらを参考に
# https://qiita.com/txt_only/items/b954e26be739bee5621e

# In[145]:


# figureを生成する

fig = plt.figure(figsize=(10, 5), dpi=100, facecolor='w', linewidth=0, edgecolor='w')
 
# 2x3の1番目
ax1 = fig.add_subplot(2, 3, 1)
ax1.set_title('SSP-NSC(log)', fontsize=18)  # グラフタイトル
ax1.set_xscale("log")
# ax1.set_xlabel("")
sns.boxplot(data=df, x='SSP-NSC', y='Antibody', sym="", ax = ax1)

# 2x3の2番目
ax2 = fig.add_subplot(2, 3, 2)
ax2.set_title('PPQT-NSC(log)', fontsize=18)  # グラフタイトル
ax2.set_xscale("log")
sns.boxplot(data=df, x='PPQT-NSC', y='Antibody', sym="", ax = ax2)

# 2x3の3番目
ax3 = fig.add_subplot(2, 3, 3)
ax3.set_title('Q-RSC', fontsize=18)  # グラフタイトル
sns.boxplot(data=df, x='Q-RSC', y='Antibody', sym="", ax = ax3)

# 2x3の4番目
ax4 = fig.add_subplot(2, 3, 4)
ax4.set_title('SSP-RSC', fontsize=18)  # グラフタイトル
sns.boxplot(data=df, x='SSP-RSC', y='Antibody', sym="", ax = ax4)

# 2x3の5番目
ax5 = fig.add_subplot(2, 3, 5)
ax5.set_title('PPQT-RSC', fontsize=18)  # グラフタイトル
sns.boxplot(data=df, x='PPQT-RSC', y='Antibody', sym="", ax = ax5)


# 2x3の6番目
ax6 = fig.add_subplot(2, 3, 6)
ax6.set_title('JSD', fontsize=18)  # グラフタイトル
sns.boxplot(data=df, x='JSD', y='Antibody', sym="", ax = ax6)
 
# 表示する

plt.tight_layout() # 文字の重なりを解消する呪文
plt.show()


# ### 小さな問題点
# - x_labelが消せない
# - xの値を上にできない
# - 色も変えないといけない
# - y_labelも２個だけでいいのに、全部ついてしまう
# - PPQT-NSC(log)のx軸の値がつぶれている

# In[ ]:


# お題６：H3K36me3, Input についてそれぞれ相関係数ヒートマップを描き、Antibodyごとの傾向の違いを観察する。


# In[147]:


df_36me3 = df[df['Antibody'] == 'H3K36me3']                                                                                                                                                                                                
df_36me3


# In[148]:


sns.clustermap(df_36me3.corr(), method='ward', cmap=cmap)


# In[149]:


df_input = df[df['Antibody'] == 'Input']                                                                                                                                                                                                
df_input


# In[150]:


sns.clustermap(df_input.corr(), method='ward', cmap=cmap)


# - 各項目の示す意味がわかっていないので、考察できない。。。

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ## 以下は作業時に試行錯誤したもの。
# ## 不要なので無視してください。

# In[64]:


# SSP_NSCのgroupbyオブジェクトをデータフレームにする。
SSP_NSC_grouped = grouped["SSP-NSC"].describe().reset_index()

# 確認
SSP_NSC_grouped


# In[109]:


SSP_NSC_grouped["count"]


# In[66]:


# 使う部分だけ取り出す
SSP_use = SSP_NSC_grouped.loc[:,['Antibody','mean']]


# In[65]:


# from groupby object to dataframe
PPQT_RSC_grouped = grouped["PPQT-RSC"].describe().reset_index()

# confirmation
PPQT_RSC_grouped


# In[67]:


# 使う部分だけ取り出す
PPQT_use = PPQT_RSC_grouped.loc[:,['mean']]


# In[68]:


# 使う部分だけ横に結合する
concat_df_v = pd.concat([SSP_use, PPQT_use], axis=1)
# column nameを付け直す
concat_df_v.columns = ['Antibody', 'SSP_NSC', 'PPQT_RSC']
# Antibodyをindexにする
concat_df_v = concat_df_v .set_index('Antibody')
# 確認
concat_df_v 


# In[89]:


# 描画
concat_df_v.plot(kind='barh', legend=True, xerr = sem_SSP)
# semをどないしよ
# CIもつけなあかん


# In[71]:


sem_SSP = group_nsc.std()/np.sqrt(group_nsc.count())
ax1 = group_nsc.mean().plot(kind='barh', xerr=sem_SSP)
ax1.invert_yaxis()


# In[73]:


group_ppqt = grouped["PPQT-RSC"]
sem_PPQT = group_ppqt.std()/np.sqrt(group_ppqt.count())
ax2 = group_ppqt.mean().plot(kind='barh', xerr=sem_PPQT)
ax2.invert_yaxis()


# In[93]:


fig = plt.figure()
ax = group_nsc.mean().plot(kind='barh', xerr=sem_SSP, align = 'edge', width = 0.2)
ax = group_ppqt.mean().plot(kind='barh', xerr=sem_PPQT, align = 'center',  width = 0.2, color = "orange")
ax.invert_yaxis()
# alignの設定は'edge'か'center'、というわけで重なってしまう。

