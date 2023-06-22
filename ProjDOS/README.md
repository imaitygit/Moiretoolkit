## Compute Order Parameter from LAMMPS file 

* Compute the Projected Density of States (PDOS) for a circular patch by running
python `pdos.py`

* Plot the data
python `plot.py`


## Installation
This code uses python3 (tested on python3.9),
+ matplotlib,
+ numpy.
You also need the `WS2.PDOS` file for a fine grid from SIESTA. 
Due to large size of the file, we haven't included the file.


## Support
Please email me [@Imperial](mailto:i.maity@imperial.ac.uk) or
[@gmail](mailto:indrajit.maity02@gmail.com) if you have any questions or have
found any bugs.


## Authors and acknowledgment
The code was developed for the paper [arXiv](https://arxiv.org/abs/2302.11497):
```
@article{molino2023influence,
  title={Influence of atomic relaxations on the moiré flat band wavefunctions in antiparallel twisted bilayer WS2},
  author={Molino, Laurent and Aggarwal, Leena and Maity, Indrajit and Plumadore, Ryan and Lischner, Johannes and Luican-Mayer, Adina},
  journal={arXiv preprint arXiv:2302.11497},
  year={2023}
}
```
