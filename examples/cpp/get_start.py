import numpy as np
from pele.systems import LJCluster
from pele.utils.xyz import write_xyz

system = LJCluster(6)
db = system.create_database()
bh = system.get_basinhopping(db)
bh.run(100)

m1, m2 = db.minima()[:2]

mindist = system.get_mindist()
d, x1, x2 = mindist(m1.coords, m2.coords)


write_xyz(open("lj6_m1.xyz", "w"), x1)
write_xyz(open("lj6_m2.xyz", "w"), x2)

np.savetxt("lj6_m1", x1)
np.savetxt("lj6_m2", x2)
