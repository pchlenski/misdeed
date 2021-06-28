from src.OmicsGenerator import OmicsGenerator

gen = OmicsGenerator(100, ['a', 'b'], [100, 25], init_full=True)
xc, yc, zc, xC, yC, zC = gen.case_control(100, 0.5, 'a', .1, dt=1e-1)
print(xc)
