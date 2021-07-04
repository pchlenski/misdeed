from src.OmicsGenerator import OmicsGenerator
from src.visualization import plot_timecourse
import numpy as np

gen = OmicsGenerator(100, ['a', 'b'], [100, 25], init_full=True)
gen.add_intervention(
    'i1',
    'b',
    vector=0.5-np.random.rand(25),
    affects_abundance=True,
    start=50,
    end=100
)
# xc, yc, zc, xC, yC, zC = gen.case_control(100, 0.5, 'a', .1, dt=1e-1)
# print(xc)

x,y,z = gen.generate()
print(x)

# trajectory 1

n_mgx = 15
n_mbx = 15

gen = OmicsGenerator(
    250,
    ['mgx', 'mbx'],
    [n_mgx, n_mbx],
    init_full=True,
    # silent=True
)

gen.add_intervention(
    'diet1',
    'mbx',
    vector=0.5-np.random.rand(n_mbx),
    affects_abundance=True,
    start=50,
    end=100
)

gen.add_intervention(
    'diet2',
    'mbx',
    vector=np.random.rand(n_mbx),
    affects_abundance=True,
    start=150,
    end=200
)

x, y, z = gen.generate()
plot_timecourse(y)