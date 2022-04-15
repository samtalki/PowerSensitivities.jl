import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matrix_sensitivities import calc_rel_err
sns.set_context("paper")
sns.set_theme(style="darkgrid")
plt.style.use("bmh")
matplotlib.rc('text', usetex=True)
matplotlib.rc('text.latex', preamble=r'\usepackage{amsmath,amssymb}')
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['legend.fontsize'] = 7
matplotlib.rcParams['lines.linewidth'] = 2




def single_performance_plot(O,S_0,pct_obs,S_tilde,S_hat):
    """Plot the sensitivity matrix recovery performance"""
    #Subplots
    fig,axes = plt.subplots(nrows=3,ncols=1, 
                            figsize=(3.5,3.5/1.61828),
                            #constrained_layout=True,
                            sharex=True, sharey=True)

    #Colorbar axis and cmap
    cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)
    cbar_ax = fig.add_axes([.85, .2, .02, .5])
    cbar_ax.tick_params(labelsize=6,pad=0.5)
    cbar_ax.get_yaxis().labelpad = -1
    
    #cbar_ax.set_label(r"$\partial V_i \big/ \partial X_j$",pad=-1)
    cbar_kws = {
        'label':r"$\partial V_i \big/ \partial X_j$",
        'fraction':0.05,
        'shrink':15,
        'pad':-10}


    rel_err = calc_rel_err(M_hat=S_hat,M=S_tilde)*100
    
    #Plot observed matrix
    #Include colorbar on first heatmap
    axes[0] = sns.heatmap(np.asarray(S_0),ax=axes[0],
                            cmap=cmap,cbar=True,vmin=-0.1,vmax=0.2,
                            cbar_ax=cbar_ax,cbar_kws=cbar_kws,
                            square=True,xticklabels=False,yticklabels=False,
                            mask=(np.asarray(S_0)==0),
                           )

    axes[0].set_title(r"{pct_obs:.0f}\%".format(pct_obs=pct_obs*100) + r" Observable $\tilde{\mathbf{S}}_0$",fontsize=7,pad=3)
    
    #Plot recovered matrix
    axes[1] = sns.heatmap(S_hat,ax=axes[1],
                            cmap=cmap,cbar=False,vmin=-0.1,vmax=0.2,
                            square=True,xticklabels=False,yticklabels=False
                         )

    axes[1].set_title(r"Est. $\tilde{\mathbf{S}}^{\#}$,"+" Rel. Err={rel_err:.1f}\%".format(rel_err=rel_err),fontsize=7,pad=3)
    
    #Plot true matrix
    axes[2] = sns.heatmap(S_tilde,ax=axes[2],
                          cmap=cmap,cbar=False,vmin=-0.1,vmax=0.2,
                          square=True,xticklabels=False,yticklabels=False
                         )
    
    axes[2].set_title(r"True $\tilde{\mathbf{S}}^{\#}$",fontsize=7,pad=3)
    axes[2].set_xlabel(r"\ \ \ Active Inj. \quad Reactive Inj.",fontsize=7)

    #Figure-level formatting
    plt.suptitle(r"IEEE 123-Bus $\tilde{\mathbf{S}}$ Matrix Recovery, $\alpha=0.9$",fontsize=9)
    return fig,axes
    #plt.savefig("/home/sam/github/PowerSensitivities.jl/figures/spring_22/IEE123_recovery_multi_lamb.125_delta.006.png",dpi=400)