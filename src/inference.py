# Inferring gLV parameters
import numpy as np

def infer_glv_params(
    abundances : np.ndarray,
    interventions : np.ndarray = None,
    interaction_reg : float = 0,
    growth_reg : float = 0,
    intervention_reg : float = 0,
    dt : float = 1e-3,
    pseudocount : float = 1e-5) -> (np.ndarray, np.ndarray, np.ndarray):
    """
    Infers interaction matrix, growth rates, and intervention responses from absolute abundance data.

    Args:
    -----
    abundances:
    interventions:
    interactions_reg:
        Float. L1 penalty for interaction matrix coefficients.
    growth_reg:
        Float. L1 penalty for growth rate coefficients.
    intervention_reg:
        Float. L1 penalty for intervention response coefficients.
    dt:
        If float: size of (uniform) timestep between observations. Equivalent to "dt" in OmicsGenerator methods.
        If ndarray: vector of time-steps for samples.
    pseudocount:
        Float. Size of pseudocount to assign to all observations to ensure log-transform works.

    Returns:
    --------
    inferred_M:
        Array. The inferred interaction matrix.
    inferred_u
        Array. The inferred growth rates.
    inferred_E
        Array. The inferred response(s) to intervention(s), if intervention vectors are provided.
    
    Raises:
    -------
    TODO

    Based on the paper "Ecological Modeling from Time-Series Inference: Insight into Dynamics and Stability of 
    Intestinal Microbiota" by Stein, Bucci, et al (2013).
    """

    # TODO: add dimension checks, make intervention size handling more robust
    # TODO: write unit tests

    # Set inferred variables
    no_interventions = False
    n_times, n_clades = abundances.shape

    # Need dummy interventions
    if interventions is None:
        interventions = np.zeros(n_times)
        no_interventions = True

    # Reshape vector-valued interventions
    elif interventions.ndim == 1:
        interventions = interventions.reshape(1, -1)

    # Build up Y matrix
    Y1 = abundances.T # read in ABSOLUTE abundances (Z matrix)
    Y2 = np.ones((1, n_times)) # 1s for each time point
    Y3 = interventions
    Y = np.concatenate((Y1, Y2, Y3), axis=0)
    Y = Y[:,:-1] # drop last time point

    # Build up F matrix
    if type(dt) is not np.ndarray:
        dt = dt * np.ones((n_times - 1, 1)) # change in times
    logs = np.log(abundances + pseudocount)
    dx = np.diff(logs, axis=0) # changes in log-abundances
    F = dx / dt
    F = F.T # need transpose here

    # TODO: test that this works without interventions
    # TODO: test that this works for custom time-steps

    # Build up lambda matrix
    lambdas = np.zeros(Y.shape[0]) # set lambda values
    lambdas[:n_clades] = interaction_reg
    lambdas[n_clades] = growth_reg
    lambdas[n_clades+1:] = intervention_reg
    lambdas = np.diag(lambdas)

    # Infer M, u, E
    MuE = F @ Y.T @ np.linalg.inv(Y @ Y.T + lambdas)
    inferred_M = MuE[:,:n_clades]
    inferred_u = MuE[:,n_clades]
    if no_interventions == True:
        inferred_E = np.zeros(n_clades)
    else:
        inferred_E = MuE[:,n_clades+1:]

    return inferred_M, inferred_u, inferred_E

