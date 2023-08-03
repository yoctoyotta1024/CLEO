### Author: Clara Bayley
### File: exmpl_anisrc.py
### classes and functions called by exmpl_plot.py
### to plot animations of some of the output
### data from the example run of CLEO

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation
from matplotlib.cm import ScalarMappable


def animate_me(fig, update_frame, frames, plot_init, saveani=False,
              savedir=None, savename=None, fargs=(), fps=5):
  
  print("making animation")
  ani = FuncAnimation(fig, update_frame, frames=range(0, frames),
                      init_func=plot_init, fargs=fargs)

  if saveani:
    print("saving animation as "+savedir+savename+".gif")
    ani.save(savedir+savename+".gif",
            writer=animation.PillowWriter(fps=fps, bitrate=5000, codec='h264'),
            dpi=100, savefig_kwargs={'transparent': True})
    
def animate1dprofile(gbxs, mom, time, nframes, 
                     xlabel=None, xlims=[None, None], color="black",
                     saveani=False, savedir=None, savename=None, fps=5):
  
  fig, ax, plots, txt, zkm = prepare_1dprofile(gbxs, mom, time,
                                               xlabel, color=color) 
  
  def init_1dprofile():
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.set_ylim([0, 1.5])
    yticks = [0, 0.75, 1.5]
    ax.set_yticks(yticks, yticks, fontsize=16)

    ax.set_xlim([xlims[0], xlims[1]*1.1])
    xticks = np.linspace(xlims[0], xlims[1], 3)
    xticklabels = xticklabels = ["{:.1f}".format(x) for x in xticks]
    ax.set_xticks(xticks, xticklabels, fontsize=16)

    ax.tick_params(length=10, width=1)

    fig.tight_layout()
    
    return plots, txt,

  animate_me(fig, update_1dprofileframe, nframes, init_1dprofile,
           saveani=saveani, savedir=savedir, savename=savename,
           fargs=(plots, txt, zkm, time, mom), fps=fps)
  
def prepare_1dprofile(gbxs, massmom, time, xlabel, color):
  
  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 8))

  alpha = 0.3
  zkm = gbxs.zfull / 1000 # convert m to km

  f = 0
  plots = ax.plot(massmom[f], zkm, color=color)[0]

  timetext = "t = {:.0f}min".format(time[f])
  txt = fig.text(0.7, 0.875, timetext, fontsize=16)

  ax.set_ylabel("z /km", fontsize=16)
  ax.set_xlabel(xlabel, fontsize=16)

  fig.tight_layout()

  return fig, ax, plots, txt, zkm

def update_1dprofileframe(f, plots, txt, zkm, time, massmom):

  plots.set_data(massmom[f], zkm)

  timetext = "t = {:.0f}min".format(time[f])
  txt.set_text(timetext)
  
  return plots, txt,

def animate2dcmap(gbxs, mom2ani, time, nframes, 
                  cbarlabel=None, cmapnorm=None, cmap="viridis",
                  saveani=False, savedir=None, savename=None, fps=5):
  
  fig, ax, cbar, plot, txt = prepare_2dplot(gbxs, mom2ani, time,
                                            cmap, cmapnorm)
  cbar.set_label(cbarlabel, fontsize=16)

  def init_2dcmap():
    
    ax.set_aspect("equal")
    ax.set_xlim([0, 1.5])
    ax.set_ylim([0, 1.5])
    
    ax.set_xlabel("x / km", fontsize=16) 
    ax.set_ylabel("z / km", fontsize=16) 

    ticks = [0, 0.75, 1.5]
    ax.set_xticks(ticks, ticks, fontsize=16)
    ax.set_yticks(ticks, ticks, fontsize=16)
    
    ax.tick_params(length=10, width=1)

    fig.tight_layout()
    return plot
   
  animate_me(fig, update_2dcmapframe, nframes, init_2dcmap,
            saveani, savedir, savename,
            fargs=(plot, txt, time, mom2ani), fps=fps)

def prepare_2dplot(gbxs, massmom, time, cmap, cmapnorm):
  
  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 9))
  
  f=0
  data2d = massmom[f,:,:]
  plot = ax.pcolormesh(gbxs.xxh/1000, gbxs.zzh/1000, data2d,
                      cmap=cmap, norm=cmapnorm)
  plot.cmap.set_under("w")
  timetext = "t = {:.0f}min".format(time[f])
  txt = fig.text(0.765, 0.925, timetext, fontsize=16, ha="right")
  
  cbar = fig.colorbar(ScalarMappable(norm=cmapnorm, cmap=cmap),
                      ax=ax, extend="both")
  cbar.ax.tick_params(labelsize=16, which="both")
  cbar.ax.tick_params(length=10, width=1, which="major")
  cbar.ax.tick_params(length=10/3, width=1, which="minor")

  fig.tight_layout()
  
  return fig, ax, cbar, plot, txt

def update_2dcmapframe(f, plot, txt, time, massmom):
    
  plot.set_array(np.array(massmom[f,:,:]).ravel())
  
  timetext = "t = {:.0f}min".format(time[f])
  txt.set_text(timetext)
  
  return plot, txt
