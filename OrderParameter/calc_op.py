#------------------------------------|
# Written by Indrajit Maity          |
# e-mail: indrajit.maity02@gmail.com |
#------------------------------------|


import sys
import numpy as np
import operator
import matplotlib.pyplot as plt
#from load_data import *


class OrderParameter(object):
  """
  Class to compute Order Parameter
  """
  def __init__(self, datafile, outfile, type1=1,\
               type2=4, r0=3.0):
    """
    Initial attributes
    @input
      Order parameters computed using two default atom types
      type1: Metal atom of first layer;
      type2: Metal atom of second layer;
    """
    self.datafile = datafile
    self.outfile = outfile
    self.type1 = type1
    self.type2 = type2
    self.r0 = r0
    # Get the positions for specific types
    self.pos_1 = self.get_specific_data(self.type1)
    self.pos_2 = self.get_specific_data(self.type2)
    self.A = self.get_lattice()


  def get_specific_data(self, atomtype):
    """
    Returns positions for specific data type
    """
    a, data = self.get_data()
    return data[np.where(data[:,0]==atomtype),1:].reshape(-1,3)


  def get_lattice(self):
    """
    Gets the lattice constants.
    """
    # Read contents of the data
    with open(self.datafile, 'r') as f:
      contents = f.readlines()
    for j in range(len(contents)):
      if "xlo xhi" in contents[j]:
        xlo = eval(contents[j].split()[0])
        xhi = eval(contents[j].split()[1])
      elif "ylo yhi" in contents[j]:
        ylo = eval(contents[j].split()[0])
        yhi = eval(contents[j].split()[1])
      elif "zlo zhi" in contents[j]:
        zlo = eval(contents[j].split()[0])
        zhi = eval(contents[j].split()[1])
      elif "xy xz yz" in contents[j]:
        xy = eval(contents[j].split()[0])
        xz = eval(contents[j].split()[1])
        yz = eval(contents[j].split()[2])
    return np.array([[xhi-xlo, 0.0, 0.0],\
                   [xy, yhi-ylo, 0.0],\
                   [xz, yz, zhi-zlo]])      


  def get_data(self):
    """
    Load the data for post-processing.
    """
    # Read contents of the data
    with open(self.datafile, 'r') as f:
      contents = f.readlines()
    for j in range(len(contents)):
      # Extract the number of atoms
      if "atoms" in contents[j]:
        natom = eval(contents[j].split()[0])
      elif "atom types" in contents[j]:
        ntype = eval(contents[j].split()[0])
      # lattice constant
      elif "xy xz yz" in contents[j]:
        a = 2*(eval(contents[j].split()[0]))
      # skip lines
      elif "Atoms # atomic" in contents[j] or\
           "Atoms" in contents[j]:
        skip = j+2
    print()
    print(f"Number of atoms in the data file: {natom}")
    print(f"Lattice constant of moiré superlattice: {a} Ang.")
    print()

    # Based on the data segregate atom types to different types
    # Quite easy to access then 
    # Type, x, y, z
    data = np.zeros((natom, 4), dtype=float) 
    for d in range(skip, natom+skip, 1):
      data[d-skip] =\
              np.array([eval(contents[d].split()[1]),\
                        eval(contents[d].split()[2]),\
                        eval(contents[d].split()[3]),\
                        eval(contents[d].split()[4])]) 
    return a, data



  def nn_tb(self):
    """
    nearest neighbor of bottom layer for the 
    two specific atom types
    """
    print(f"r0 to search for unit-cell lattice: {self.r0} Ang.")
    print("Modify r0 according to your needs")
    f = open(self.outfile, "w")
    f.write("%s %s %s %s %s %s %s\n"%("type1_x", "type1_y", "type1_z",\
                                   "type2_x", "type2_y", "type2_z",\
                                   "distance"))
    # No periodic boundary conditions implemented
    for i in range(self.pos_1.shape[0]):
      # 0:x; 1:y; 2:z; 3:distance 
      tmp = []
      for j in range(self.pos_2.shape[0]):
        rvec = self.get_distance(self.pos_1[i,:2], self.pos_2[j,:2])
        if np.linalg.norm(rvec) <= self.r0:
          tmp.append(rvec.tolist())
      if not tmp:
        print("Distance list is empty!")
        print("Increase r0. Exiting...")
        sys.exit()
      else:
        tmp = np.array(tmp)
      tmp_mod = np.sqrt(tmp[:,0]**2 + tmp[:,1]**2)
      indx, val = min(enumerate(tmp_mod[:,]), key=operator.itemgetter(1))
      f.write("%.6f %.6f %.6f %.6f %.6f %.6f\n"%(self.pos_1[i,0],\
               self.pos_1[i,1], self.pos_1[i,2], tmp[indx,0],\
               tmp[indx,1], tmp_mod[indx]))
    f.close()
    print("Written data to file") 


  def get_distance(self, r1, r2):
    """
    In-plane distance vector with PBC
    """
    Ainv = np.linalg.inv(self.A[:2,:2])
    p2 = np.dot(r1, Ainv[:2,:2])
    p1 = np.dot(r2, Ainv[:2,:2])
    d = p2-p1
    p2_un = np.zeros(p2.shape[0], dtype=float)
    for i in range(p2.shape[0]):
      if np.abs(d[i]) >= 0.5:
        if p2[i] >= p1[i]:
          p2_un[i] = p2[i] - 1.0
        else :
          p2_un[i] = p2[i] + 1.0
      else:
        p2_un[i] = p2[i]

    r1 = np.dot(p1, self.A[:2,:2])  
    r2 = np.dot(p2_un, self.A[:2,:2])  
    return r2-r1


  def plt_op(self):
    alpha = 0.25
    plt.figure(figsize=(9.6, 4))
    x0, y0, z0, rx1, ry1, d = \
                  np.loadtxt(self.outfile, skiprows=1, unpack=True)
    Q = plt.quiver(x0, y0, rx1, ry1, d,units='xy', cmap='plasma_r')
    plt.xlabel(r"x ($\AA$)", fontsize=16)
    plt.ylabel(r"y ($\AA$)", fontsize=16)
    plt.colorbar()
    plt.savefig("OP.png", dpi=300, bbox_inches="tight", pad_inches=0.1)
    plt.show()


# RUN the code 
OP = OrderParameter("lammps.dat", "S_inter_nn_new")
#OP.nn_tb()
OP.plt_op()
