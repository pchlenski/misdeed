# Plotting helper functions
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from sklearn.decomposition import PCA

def plot_timecourse(
    df, 
    cols=15) -> None:

    """
    Makes a bar chart of longitudinal data from a generator

    Args:
        df:
            A pandas dataframe (or dtype that can be cast as a dataframe) with relative abundance timecourse data.
            Assumes rows = time points, columns = dimensions (e.g. taxa)
        cols:
            Int. Number of columns in legend.
        
    Returns:
        A stacked bar chart showing a colored relative abundance timecourse.
    
    Raises:
        TODO
    """

    # make dataframe automatically
    if type(df) is not pd.DataFrame:
      df = pd.DataFrame(df)

    # get colors
    n_colors = df.shape[1]
    viridis = cm.get_cmap('magma', n_colors)
    colors  = viridis(range(n_colors))
    plots = df.plot.bar(stacked=True, figsize=(40,10), width=1.0, color=colors)
    plt.legend(loc='lower left', bbox_to_anchor=(0,1,1,0), ncol=cols, mode='expand')

    # tweak params
    plots.patch.set_visible(False)               # make background transparent
    plots.spines['top'].set_visible(False)       # hide top line of plot
    plots.spines['right'].set_visible(False)     # hide right line of plot
    plt.xlim(-0.5, df.shape[0] -0.5)             # tight x axis
    plt.ylim(0, 1)                               # tight y axis
    plt.xlabel("Sample ID")
    plt.ylabel("Relative abundance")

    return plots

def plot_pca(
    trajectories, 
    node_name, 
    colors=False, 
    cmap='magma', 
    plot_endpoints=True, 
    plot_trajectories=True,
    **kwargs) -> None:

    """
    Makes PCA plots from (a list of) generator outputs.

    Args:
        trajectories:
            A list of multi-individual samples (output from generate_multiple()). Each timecourse is expected to be a
            list of dicts keyed by node names.
        node_name:
            String. The name of the node to extract from the inputs.
        colors:
            A list of colors for each input. Only used when #{inputs} > 1.
        cmap:
            A cmap to be used for plotting individual trajectories. Only used when #{inputs} = 1.
        plot_endpoints:
            Boolean. If True, marks endpoints of trajectories with an X.
        plot_trajectories:
            Boolean. If True, plots entire timecourse as a PCA trajectory.
        **kwargs:
            Passed to plt.plot() and plt.scatter().
    
    Returns:
        A PCA plot of trajectories for each individual in the input set. 
        
        If multiple inputs are presented, then individuals will be colored according to the 'colors' parameter. If only
        one input is presented, then individuals will be colored individually according to the 'cmap' parameter.
    
    Raises:
        TODO
    """

    # Flatten list of lists, get node values
    node_inputs = [sample[node_name] for samples in trajectories for sample in samples]

    # PCA on list of lists
    pca = PCA(n_components=2)
    node_size = len(node_inputs[0][0]) # dimension of first timepoint of first individual
    node_inputs = np.array(node_inputs).reshape(-1, node_size)
    node_inputs = np.nan_to_num(node_inputs)
    pca = pca.fit(node_inputs)

    # If plotting non-clustered data, generate cmap
    if len(trajectories) == 1:
        n_samples = len(trajectories[0])
        temp_cmap = cm.get_cmap(cmap, n_samples)
        colormap = temp_cmap(range(n_samples + 1)) # +1 to avoid really light colors

    # If plotting clustered data without specified clusters, generate cmap
    elif colors == False:
        n_clusters = len(trajectories)
        temp_cmap = cm.get_cmap(cmap, n_clusters)
        colormap = temp_cmap(range(n_samples + 1))
        colors = colormap.colors

    # Start plotting inputs
    input_counter = 0

    for input in trajectories:
        node_input = [sample[node_name] for sample in input]
        individual_counter = 0

        for individual in node_input:
            individual = np.array(individual).reshape(-1, node_size)
            individual = np.nan_to_num(individual)
            individual_pca = pca.transform(individual)

            # color case 1: all inputs belong to the same cluster
            if len(trajectories) == 1:
                c = colormap[individual_counter]
            
            # color case 2: different clusters
            else:
                c = colors[input_counter]
            
            # plot PCA trajectories we need
            if plot_trajectories == True:
                plt.plot(individual_pca[:,0], individual_pca[:,1], color=c, **kwargs)

            # draw Xes
            if plot_endpoints == True:
                plt.scatter(individual_pca[-1,0], individual_pca[-1,1], color=c, marker='x', **kwargs)
            
            individual_counter += 1

        input_counter += 1
