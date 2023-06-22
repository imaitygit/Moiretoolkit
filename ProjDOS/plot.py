import numpy as np
import matplotlib.pyplot as plt


## Orbitals
#E = np.load("E_9.npy")
#DOS = np.load("DOS_9.npy")
#DOS_UD = DOS[:,0] + DOS[:,1]
#plt.plot(E, DOS_UD/27142, color="r",\
#         linewidth=2, ls = "--", label=r"DOS")

plt.figure(figsize=(4.07, 4.5))
sh = 4.02+0.134
#sh = 0.0
E_BMM = np.load("E_BMM_9.npy")
DOS_BMM = np.load("DOS_BMM_9.npy")
DOS_UD_BMM = DOS_BMM[:,0] + DOS_BMM[:,1] #+ DOS_BMM[:,2] + DOS_BMM[:,3]
#plt.semilogy(E_BMM+sh, DOS_UD_BMM/2334, color="g",\
#         linewidth=2.5, label=r"$\mathrm{B^{W/W}}$")
#
E_BXX = np.load("E_BXX_9.npy")
DOS_BXX = np.load("DOS_BXX_9.npy")
DOS_UD_BXX = DOS_BXX[:,0] + DOS_BXX[:,1] #+ DOS_BXX[:,2] + DOS_BXX[:,3]
plt.semilogy(E_BXX+sh, DOS_UD_BXX/2422, color="r",\
         linewidth=2, label=r"$\mathrm{B^{S/S}}$")

E_2H = np.load("E_2H_9.npy")
DOS_2H = np.load("DOS_2H_9.npy")
DOS_UD_2H = DOS_2H[:,0] + DOS_2H[:,1] #+DOS_2H[:,2]+DOS_2H[:,3]
plt.semilogy(E_2H+sh, DOS_UD_2H/2378, color="b",\
         linewidth=2.5, label=r"$\mathrm{AA^\prime}$")

#avg_dat = (DOS_UD_2H/2378) + (DOS_UD_BXX/2422) + (DOS_UD_BMM/2334)
#plt.semilogy(E_2H+sh, avg_dat, color="c",\
#         linewidth=2.5, label=r"AVG")

plt.legend(frameon=False, loc=3, fontsize=14)
plt.axvline(0.0, linestyle="--", lw=2.5, alpha=0.6, color="k")

plt.ylim(0.0003, 0.2)
plt.xlim(-1.3, 1.8)
plt.yticks([])
plt.xticks(fontsize=14)
plt.ylabel(r"PDOS (a.u.)", fontsize=14)
plt.xlabel(r"E (eV)", fontsize=14)
plt.savefig("pdos_midpoint.png", dpi=300, bbox_inches="tight",\
            pad_inches=0.1)
plt.show() 
