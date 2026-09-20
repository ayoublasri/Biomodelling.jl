import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
HERE="/home/user/Biomodelling.jl/paper"; OUT=os.path.join(HERE,"output"); PUB="/home/user/Biomodelling.jl/video/public"
for w in ("400","600","800"):
    fp=os.path.join(HERE,"figures","fonts",f"Inter-{w}.ttf")
    if os.path.exists(fp): font_manager.fontManager.addfont(fp)
FAM="Inter" if any(f.name=="Inter" for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
C=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4"]; INK,INK2,GRID="#0b0b0b","#575652","#dedcd6"
plt.rcParams.update({"font.family":FAM,"font.size":25,"axes.labelsize":26,"xtick.labelsize":23,"ytick.labelsize":23,
    "legend.fontsize":23,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":1.8,
    "axes.edgecolor":INK2,"xtick.color":INK2,"ytick.color":INK2,"axes.labelcolor":INK,"text.color":INK,
    "legend.frameon":False,"axes.grid":True,"grid.color":GRID,"grid.linewidth":1.2,"axes.axisbelow":True,
    "figure.facecolor":"white","savefig.facecolor":"white","savefig.dpi":150,"axes.labelpad":8})
FS=(9.2,4.0)
def save(fig,name):
    fig.savefig(os.path.join(PUB,name),bbox_inches="tight",pad_inches=0.3); plt.close(fig)
    from PIL import Image; print(name, Image.open(os.path.join(PUB,name)).size)

d=pd.read_csv(os.path.join(OUT,"fig4b_killcurves.csv"))
fig,ax=plt.subplots(figsize=FS)
for i,(dose,g) in enumerate(d[d.dose>0].groupby("dose")):
    ax.plot(g.t_since_drug,g.surviving_fraction,lw=4.0,color=C[i],label=f"dose {dose:g}")
ax.set_yscale("log"); ax.set_xlabel("time since drug"); ax.set_ylabel("surviving fraction")
ax.legend(ncol=4,loc="upper center",bbox_to_anchor=(0.5,-0.26),columnspacing=1.6,handlelength=1.7)
save(fig,"chart_kill.png")

d=pd.read_csv(os.path.join(OUT,"fig3c_heritability.csv")).sort_values("k_switch")
fig,ax=plt.subplots(figsize=FS)
for i,(col,lab) in enumerate((("mother_daughter","mother–daughter"),("sisters","sisters"),("cousins","cousins"))):
    ax.plot(d.k_switch,d[col],marker="o",ms=9,lw=4.0,color=C[i],label=lab)
ax.axvline(np.log(2)/20,color=INK2,ls=":",lw=2.2)
ax.text(np.log(2)/20*1.18,0.9,"1 / cycle",fontsize=21,color=INK2)
ax.set_xscale("log"); ax.set_xlabel("promoter switching rate"); ax.set_ylabel("correlation"); ax.set_ylim(-0.08,1.05)
ax.legend(ncol=3,loc="upper center",bbox_to_anchor=(0.5,-0.26),columnspacing=1.8,handlelength=1.7)
save(fig,"chart_memory.png")

d=pd.read_csv(os.path.join(OUT,"fig7a_fates.csv")); m=d.groupby("cisplatin_uM").mean(numeric_only=True).reset_index()
fig,ax=plt.subplots(figsize=FS)
x=np.arange(len(m)); w=0.26
for j,(sim,obs,name) in enumerate((("died","obs_died","died"),("divided","obs_divided","divided"),("survived","obs_survived","survived"))):
    pos=x+(j-1)*w
    ax.bar(pos,m[sim],width=w*0.86,color=C[j],label=name,zorder=2)
    ax.scatter(pos,m[obs],marker="_",s=900,linewidths=4.0,color=INK,zorder=4,label="observed" if j==0 else None)
ax.set_xticks(x); ax.set_xticklabels([f"{v:g} µM"+("\n(held out)" if v==10 else "") for v in m.cisplatin_uM])
ax.set_ylabel("fraction of cells"); ax.set_ylim(0,0.80)
ax.legend(ncol=4,loc="upper center",bbox_to_anchor=(0.5,-0.30),columnspacing=1.4,handlelength=1.4)
save(fig,"chart_fates.png")

d=pd.read_csv(os.path.join(OUT,"fig8b_melanoma_trajectories.csv")); d=d[d.mechanism=="partial protection"]
fig,ax=plt.subplots(figsize=FS)
for i,(s,g) in enumerate(d.groupby("schedule",sort=False)):
    ax.plot(g.t_weeks,g.N_over_N0,lw=4.0,color=C[i],label=s.replace(" (S1320)","").replace(" (50 %)",""))
ax.axhline(1.2,color=INK2,ls=":",lw=2.2); ax.axvline(8,color=INK2,ls="--",lw=1.8)
ax.set_yscale("log"); ax.set_xlabel("weeks"); ax.set_ylabel("tumour / initial")
ax.legend(ncol=3,loc="upper center",bbox_to_anchor=(0.5,-0.26),columnspacing=1.6,handlelength=1.7)
save(fig,"chart_schedules.png")

d=pd.read_csv(os.path.join(OUT,"fig7h_memory_slice.csv"))
piv=d.pivot(index="p_on",columns="memory_generations",values="rmse")
fig,ax=plt.subplots(figsize=(8.0,4.0))
im=ax.imshow(piv.values,origin="lower",aspect="auto",cmap="viridis_r",
             extent=[-0.5,piv.shape[1]-0.5,-0.5,piv.shape[0]-0.5])
best=d.rmse.min()
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        if piv.values[i,j]<=best+0.028:
            ax.add_patch(plt.Rectangle((j-0.5,i-0.5),1,1,fill=False,ec="#e34948",lw=3.0))
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels([f"{c:g}" for c in piv.columns],fontsize=21)
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels([f"{int(r*100)}%" for r in piv.index],fontsize=21)
ax.set_xlabel("memory (generations)"); ax.set_ylabel("cells resistant"); ax.grid(False)
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.set_label("error",fontsize=22); cb.ax.tick_params(labelsize=19)
save(fig,"chart_ident.png")
