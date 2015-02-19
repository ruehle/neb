import numpy as np
from pele.systems import LJCluster
from pele.utils.xyz import write_xyz

natoms=38
system = LJCluster(natoms)
db = system.create_database()
bh = system.get_basinhopping(db)
bh.run(1000)

m1, m2 = db.minima()[:2]

if False:
    mindist = system.get_mindist()
    d, x1, x2 = mindist(m1.coords, m2.coords)
else:
    x1, x2 = m1.coords, m2.coords


#write_xyz(open("lj6_m1.xyz", "w"), x1)
#write_xyz(open("lj6_m2.xyz", "w"), x2)
np.savetxt("lj{}_m1".format(natoms), x1.reshape(-1,3))
np.savetxt("lj{}_m2".format(natoms), x2.reshape(-1,3))
