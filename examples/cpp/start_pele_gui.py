import numpy as np

from pele.systems import LJCluster
from pele.gui import run_gui

natoms = 38
system = LJCluster(natoms)
pot = system.get_potential()
db = system.create_database()


x1 = np.genfromtxt("lj{}_m1".format(natoms)).ravel()
#x2 = np.genfromtxt("lj{}_m2".format(natoms)).ravel()
x2 = np.genfromtxt("lj{}_m3".format(natoms)).ravel()

db.addMinimum(pot.getEnergy(x1), x1)
db.addMinimum(pot.getEnergy(x2), x2)

run_gui(system, db)
