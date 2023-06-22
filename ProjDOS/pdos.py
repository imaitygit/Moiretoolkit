#--------------------------------|
# Author: Indrajit Maity         |
# email: i.maity@imperial.ac.uk  |
#--------------------------------|

import numpy as np
import matplotlib.pyplot as plt
import sys


class selected_PDOS(object):
  """
  Projected DOS calculator for selected area
  aound a point "origin" 
  """
  
  def __init__(self, alat, pdosfile, origin, radius, out1, out2):
    """
    Initial attributes
    @input
      pdosfile: Projected density of states files
      origin: The center around which the circular area is
              to be constructed
      radius: Radius in Ang. around origin.
      out1: Energy output file
      out2: Extracted PDOS for the local patch      
    """
    self.alat = alat
    self.pdosfile = pdosfile
    self.origin = origin
    self.radius = radius
    self.out1 = out1
    self.out2 = out2 
    self.A = self.set_A()


  def set_A(self):
    """
    Lattice constants manually entered
    as we are not reading any other files
    """
    A = np.array([[1.0, 0.0, 0.0],\
                    [0.5, 0.8660254, 0.0],\
                    [0.0, 0.0, 0.43211452819238777]])\
                    * self.alat*1.88972687
    return A


  def get_energy(self):
    """
    Get the energies at which DOS are computed
    To DO: Add the ability to extract norbitals and energy
           points
    @output
      E_d: Eeturns the energies for the dos
    """
    f = open(self.pdosfile, "r")
    lines = f.readlines()
    f.close()
    
    E_d = np.zeros((400), dtype=float)
    for i in range(len(lines)):
      if "energy_values" in lines[i] and\
         "units" in lines[i]:
        for j in range(i+1, i+1+E_d.shape[0]):
          E_d[j-i-1] = eval(lines[j].split()[0])
        break
    return E_d 


  def get_pdos(self):
    """
    Read and extracts pdos informations
    @output
      norb: Number of orbitals
      p: Positions that include atom index, type(string),
         atomic positions in angstroms 
    """
    f = open(self.pdosfile, "r")
    lines = f.readlines()
    f.close()

    # Get indices of the lines to read
    # Based on atom indices
    info_idx = []
    for i in range(len(lines)):
      if "atom_index" in lines[i]:
        info_idx.append(i)

    diff = info_idx[1]-info_idx[0]
    atom_id = []; atom_type = []
    p_sel = []
    counter = 0
    ne = 400
    sel_dos = np.zeros((ne, 4))
    # Run everything on the fly
    for i in range(len(info_idx)):
      # For each atom index:
      # -> Extract atomic id, atomic type, positions
      # -> Find if the atom is within the specified 
      #    region using the periodic boundary cond. 
      # -> If inside the boundary then include in the pdos;
      for j in range(info_idx[i], info_idx[i]+diff-1):
        if "atom_index" in lines[j]:
          atom_id.append(int(lines[j].split('"')[1]))
        elif "species" in lines[j]:
          atom_type.append(str(lines[j].split('"')[1]))
        elif "position" in lines[j]:
          content = lines[j].split('"')[1].split()
          pos = [eval(content[0]), eval(content[1]),\
                 eval(content[2])]
          # Compute the position in crystal coordinates
          orig_c = self.get_poscrys(origin)
          pos_c_un = self.get_poscrys(pos)
          #pos_c = self.get_wrap(pos_c_un)
          # Bohr-to-Angstrom
          pos_c = pos_c_un
          # Find out the distance between two-points
          #if self.dist_c_basic(orig_c, pos_c) < self.radius:
          if self.dist_c(orig_c, pos_c) < self.radius:
            counter = counter + 1
            p_sel.append(pos)
            for k in range(info_idx[i], info_idx[i]+diff-1):
              if "<data>" in lines[k]:
                for l in range(k+1, k+1+ne):
                  tmp_dos = np.array([eval(lines[l].split()[0]),\
                                      eval(lines[l].split()[1]),\
                                      eval(lines[l].split()[2]),\
                                      eval(lines[l].split()[3])])
                  sel_dos[l-k-1] = sel_dos[l-k-1] + tmp_dos
    E_d = self.get_energy()
    print("Saving file")
    print()
    np.save(self.out1, E_d)
    np.save(self.out2, sel_dos)
    print("Coordinate:", self.origin)
    print("orbitals within the area: ", counter)
    print("-------------------")
#    plt.plot(E_d, sel_dos[:,0])
#    plt.show()
#    #print("There are %d points"%(counter)) 
#    plt.scatter(np.array(p_sel)[:,0], np.array(p_sel)[:,1])
#    zero = np.zeros((2))
#    lw = 2
#    plt.plot([zero[0],self.A[0][0]+zero[0]],\
#          [zero[1], self.A[0][1]+zero[1]],\
#          color = "k", linewidth=lw, ls='--')
#    plt.plot([zero[0],self.A[1][0]+zero[0]],\
#         [zero[1], self.A[1][1]+zero[1]],\
#         color = "k",linewidth=lw, ls ='--')
#    Sum = self.A[0] + self.A[1]
#    plt.plot([self.A[1][0]+zero[0], Sum[0]+zero[0]],\
#         [self.A[1][1]+zero[1], Sum[1]+zero[1]],\
#         color = "k",linewidth=lw, ls='--')
#    plt.plot([self.A[0][0]+zero[0], Sum[0]+zero[0]],\
#         [self.A[0][1]+zero[1], Sum[1]+zero[1]],\
#         color = "k",linewidth=lw, ls='--')    
#    plt.show() 
  
          

  def get_wrap(self, pos_c_un):
    """
    Wrap the crystal coordinates by subtracting 1 
    or adding 1
    """
    pos_c = np.copy(pos_c_un)
    for i in range(pos_c_un.shape[0]):
      if pos_c[i] >= 1.0:
        print("Why crystal coordinates have larger value")
        print("Exiting..")
        sys.exit()
        pos_c[i] = pos_c_un[i] - 1.0
      elif pos_c[i] < 0.0:
        pos_c[i] = pos_c_un[i] + 1.0
        print("Why crystal coordinates have lesser value")
        print("Exiting..")
        sys.exit()
    return pos_c


  def get_poscrys(self, pos):
    """
    Compute the position in crystal coordinates
    for all the atoms;
    @input
      A: lattice vectors
      pos: (x, y, z) for individual atom
    @output
      pos_c: Position in crystal coordinates
  """
    pos_c = np.dot(pos, np.linalg.inv(self.A))
    return pos_c


  def dist_c_basic(self, p1,p2):
    """
    Computes and returns the distance between two points
    without the periodic boundary conditions in 2D.
    @input
      p1,p2: Points in crystal coordinates
    """
    d = (p2-p1)[:2]
    dist = np.dot(d, self.A[:2,:2])
    return np.linalg.norm(dist)


  def dist_c(self, p1,p2):
    """
    Computes and returns the distance between two points
    including the periodic boundary conditions in 2D.
    @input
      p1,p2: Points in crystal coordinates
    """
    d = p2-p1
    if np.abs(d[0]) >= 0.5:
      if p2[0] >= p1[0]:
        x2_un = p2[0] - 1.0
      else :
        x2_un = p2[0] + 1.0
    else:
      x2_un = p2[0]
  
    if np.abs(d[1]) >= 0.5:
      if p2[1] >= p1[1]:
        y2_un = p2[1] - 1.0
      else :
        y2_un = p2[1] + 1.0
    else:
      y2_un = p2[1]
    d = np.array([x2_un, y2_un])-p1[:2]
    dist = np.dot(d, self.A[:2,:2])
    return np.linalg.norm(dist)


  def read_lammpspos(self):
    """
    Reads a lammps pos file and returns the necessary
    informations
    """
    f = open(self.posfile)
    lines = f.readlines()
    f.close()
    
    for i in range(len(lines)):
      if "atoms" in lines[i]:
        natom = int(lines[i].split()[0])
      elif "atom type" in lines[i]:
        at_type = int(lines[i].split()[0])
      elif "xlo xhi" in lines[i]:
        xlo = eval(lines[i].split()[0])
        xhi = eval(lines[i].split()[1]) 
      elif "ylo yhi" in lines[i]:
        ylo = eval(lines[i].split()[0])
        yhi = eval(lines[i].split()[1]) 
      elif "zlo zhi" in lines[i]:
        zlo = eval(lines[i].split()[0])
        zhi = eval(lines[i].split()[1]) 
      elif "xy xz yz" in lines[i]:
        xy = eval(lines[i].split()[0])
        xz = eval(lines[i].split()[1]) 
        yz = eval(lines[i].split()[2]) 
      # information on atoms
      elif "Atoms # atomic" in lines[i] or\
           "Atoms" in lines[i]:
        p = np.empty((natom, 5), dtype=object)
        for k in range(i+2, i+2+natom):
          for l in range(2):
            print(lines[k])
            p[k-i-2][l] = int(lines[k].split()[l])
          for l in range(2, 5):
            p[k-i-2][l] = eval(lines[k].split()[l])

    A = np.array([[(xhi-xlo), 0.0, 0.0],\
                   [xy, yhi-ylo, 0.0],\
                   [xz, yz, zhi-zlo]])
    return natom, at_type, A, p


        



  def myregion(self):
    """
    Returns all the atoms within a given radius 
    from a user-defined origin; All of these are 
    planar calculations.
    @input
      origin: Position of the origin
      radius: Find all the atoms within this circular
              area with the specified radius. Periodic
              Boundary Condition is implemented.
      posfile: Position of all the atoms with atom ids, labels;
    @output
      region: Atoms that fall within the specified area
    """
    natom, at_type, A, p = self.read_pos()
    p_c = self.get_poscrys(p, A)
    return p_c


#--------------------- RUN stuffs ----------------
alat = 57.8550323327 # in Angstrom.
## BXX
print("BXX stacking")
origin = np.array([28.9, 16.7, 0.0])*1.88972687
# Angstrom to Bohr conversion
radius = 9.0*1.88972687
mypdos = selected_PDOS(alat, "WS2.PDOS",origin, radius,\
                       "E_BXX_9.npy", "DOS_BXX_9.npy") 
mypdos.get_pdos()

# --- 2H/AA' stacking ---
origin = np.array([58.7, 35, 0.0])*1.88972687
print("2H stacking")
mypdos = selected_PDOS(alat, "WS2.PDOS",origin, radius,\
                       "E_2H_9.npy", "DOS_2H_9.npy") 
mypdos.get_pdos()

# --- BMM stacking ---- 
print("BWW")
origin = np.array([0.0, 0.0, 0.0])*1.88972687
mypdos = selected_PDOS(alat, "WS2.PDOS",origin, radius,\
                       "E_BMM_9.npy", "DOS_BMM_9.npy") 
mypdos.get_pdos()
