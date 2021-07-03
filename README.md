# ubiome-sim
Microbiome data simulator for power analysis

<!-- generator schematic goes here -->
<!-- example synthetic data goes here -->

The OmicsGenerator class has the following methods:
* `add_node`: Adds nodes to generator object.
* `add_interaction`: Adds interactions to generator object.
* `add_intervention`:  Adds an intervention to generator.
* `set_initial_value`: Sets a node value or growth rate.
* `get`: Gets a (node/interaction/intervention) by name.
* `remove`: Removes a node, intervention, or interaction from the generator by name.
* `generate`: Generates a single timecourse of synthetic data.
* `generate_multiple`: Generates several timecourses of synthetic data.
* `case_control`: Generates synthetic case and control timecourses.
* `copy`: Makes a deep copy of generator.

## Dependencies
Dependencies are listed in `requirements.txt` and can be installed using
```bash
pip install -r requirements.txt
```

## Examples
The corresponding jupyter notebook can be found at `notebooks/examples.ipynb`.

### Initialize generator
```python
# initialize generator:
gen = OmicsGenerator(
    100,                   # 100 time points
    ['mgx', 'mbx'],        # 2 nodes named 'mgx' and 'mbx'
    [15, 15],              # each node has 15 dimensions
    init_full=True         # set interaction matrices and growth rates randomly
)

# add intervention:
gen.add_intervention(
    'intervention1',       # intervention name
    'mgx',                 # apply to 'mgx' node
    10*np.random.rand(15), # set intervention response vector randomly
    start=50,              # start at t=50
    end=100                # go to end
)

```

### Single timecourse
```python
# run generator and plot:
x1, y1, z1 = gen.generate(dt=1e-2)
plot_timecourse(y1['mgx'])
plt.vlines(50, 0, 1)
```
![Single timecourse](./img/ex1.png)

### Multiple timecourses
```python
# run multi-generator and plot:
x2, y2, z2 = gen.generate_multiple(20)
plot_pca([y2], 'mgx')
```
![Multiple timecourses](./img/ex2.png)

### Case-control
```python
# run case-control and plot:
x3_control, y3_control, z3_control, x3_case, y3_case, z3_case = gen.case_control(100, .75, 'mgx', 1)
plot_pca([y3_control, y3_case], 'mgx', colors=['red', 'blue'], plot_trajectories=False)
```
![Case-control](./img/ex3.png)

### Using learned interaction matrices
TODO

## Citation
TODO