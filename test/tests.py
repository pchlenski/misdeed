from src.OmicsGenerator import OmicsGenerator
import numpy as np

def test_1():
    # UNIT TEST 1 : interaction deletion
    gen1 = OmicsGenerator(
        time_points=200, 
        nodes=["mgx", "mbx"],
        node_sizes=[30, 30]
    )

    gen1.add_interaction(
        'F', 
        "mbx", 
        "mgx", 
        0.5-np.random.rand(30, 30), 
        lag=1
    )
    f = gen1.get('F')
    print(f)
    print('inbound node', f.inbound_node)
    print('outbound node', f.outbound_node)
    print('inbound<-', f.inbound_node.outbound.keys())
    print('outbound->', f.outbound_node.inbound.keys())

    gen1.remove("F")
    print(gen1)
    for n in gen1.nodes:
        print(n)
        print("IN:  ", n.inbound)
        print("OUT: ", n.outbound)
        print()

def test_2():
    # UNIT TEST 2: intervention removal

    n_mgx = 30
    n_mbx = 30
    gen2 = OmicsGenerator(
        time_points=200, 
        nodes=["mgx", "mbx"], 
        node_sizes=[30, 30]
    )

    gen2.add_intervention(
        'B', 
        "mgx", 
        1*(0.5-np.random.rand(30)), 
        affects_abundance=False, 
        start=40, 
        end=80
    )
    print(gen2)
    print([x.name for x in gen2.get("mgx").interventions])

    gen2.remove("B")
    print(gen2)
    print([x.name for x in gen2.get("mgx").interventions])

def test_3():
    # UNIT TEST 3: node removal

    n_mgx = 30
    n_mbx = 30
    gen3 = OmicsGenerator(
        time_points=200,
        nodes=["mgx", "mbx"], 
        node_sizes=[n_mgx, n_mbx]
    )

    mgx = gen.get("mgx")
    mgx.log_noise = False

    gen3.add_interaction(
        'A', 
        "mgx", 
        "mgx", 
        0.5-np.random.rand(n_mgx, n_mgx)
    )
    gen3.add_intervention(
        'B', 
        "mgx", 
        1*(0.5-np.random.rand(n_mgx)), 
        affects_abundance=False, 
        start=40, 
        end=80
    )
    gen3.add_intervention(
        'C', 
        "mgx", 
        1*(0.5-np.random.rand(n_mgx)), 
        affects_abundance=True,  
        start=120, 
        end=160
    )

    gen3.add_interaction(
        'D', 
        "mbx", 
        "mbx",  0.5-np.random.rand(n_mbx, n_mbx)
    )
    gen3.add_interaction(
        'E', 
        "mgx", 
        "mbx",  0.5-np.random.rand(n_mgx, n_mbx)
    )
    gen3.add_interaction(
        'F', 
        "mbx", 
        "mgx",  0.5-np.random.rand(n_mbx, n_mgx), 
        lag=1
    )
    f = gen3.get('F')
    print(f)
    print('inbound node', f.inbound_node)
    print('outbound node', f.outbound_node)
    print('inbound<-', f.inbound_node.outbound.keys())
    print('outbound->', f.outbound_node.inbound.keys())

    # gen.remove("mgx")
    gen3.remove('mgx')

    print(gen3)

def test_4(n_nodes : int) -> None:
    # UNIT TEST 4: test set_interactions for n nodes

    nodes = []
    for i in range(n_nodes):
        nodes.append(f"n{i}")

    test_gen = OmicsGenerator(
        time_points=100, 
        nodes=nodes, 
        node_sizes=np.random.randint(0, 999, size=n_nodes), 
        silent=True
    )
    print(test_gen)
    print('***************************')
    test_gen.set_interactions()
    print(test_gen)

def test_5(
    t : int = 1000, 
    n : int = 10, 
    a_size : int = 10, 
    b_size : int = 20, 
    nv : float = 0, 
    dt : float = 1e-2) -> "pca plot":
    # UNIT TEST 5: test plot_pca

    gen1 = OmicsGenerator(
        time_points=t, 
        nodes=['a', 'b'], 
        node_sizes=[a_size, b_size], 
        init_full=True, 
        silent=True
    )
    x1, y1, z1 = gen1.generate_multiple(n, noise_var=nv, dt=dt)
    
    gen1.add_intervention('c', 'a', np.random.rand(a_size), start=0, end=t)
    x2, y2, z2 = gen1.generate_multiple(n, noise_var=nv, dt=dt)
    
    gen1.add_intervention('d', 'b', np.random.rand(b_size), start=0, end=t)
    x3, y3, z3 = gen1.generate_multiple(n, noise_var=nv, dt=dt)

    plot_pca([x1, x2, x3], 'a', alpha=0.4)

def test_6():
    # UNIT TEST 6: EXPORT CSV, SINGLE TIMECOURSE
    gen = OmicsGenerator(
        time_points=20, 
        nodes=['a', 'b'], 
        node_sizes=[10, 10], 
        init_full=True, 
        silent=True
    )
    x,y,z = gen.generate()
    gen.save(x)

def test_7():
    # UNIT TEST 7: EXPORT CSV, MULTIPLE TIMECOURSES
    gen = OmicsGenerator(
        time_points=20, 
        nodes=['a', 'b'], 
        node_sizes=[10, 10], 
        init_full=True, 
        silent=True
    )
    x,y,z = gen.generate_multiple(10)
    gen.save(x)

def test_8():
    # UNIT TEST 8: EXPORT CSV, MULTIPLE TIMECOURSES, CUSTOM PATH
    gen = OmicsGenerator(
        time_points=20, 
        nodes=['a', 'b'], 
        node_sizes=[10, 10], 
        init_full=True, 
        silent=True
    )
    x,y,z = gen.generate_multiple(10)
    gen.save(x, "/tmp/abc", prefix="abc_")

def test_9():
    # Check that means balance out
    gen = OmicsGenerator(
        time_points=200, 
        nodes=['a', 'b'], 
        node_sizes=[10, 10], 
        init_full=True, 
        silent=True
    )
    x,y,z = gen.generate(downsample=5)
    print("X:", x['a'].sum(axis=1))
    print("Y:", y['a'].sum(axis=1))
    print("Z:", z['a'].sum(axis=1))

