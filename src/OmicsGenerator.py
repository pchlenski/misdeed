import numpy as np
from scipy.integrate import solve_ivp
from copy import deepcopy

class OmicsGenerator:

    """
    Handles all omics generation.

    This class is used to specify omics generation parameters and generate synthetic data. Typical workfolow is:
        Initialize generator -> set interactions -> set interventions -> generate synthetic data

    Attributes:
    -----------
    nodes:
        List of nodes.
    """

    def __init__(
        self, 
        time_points : int, 
        nodes : list = [], 
        node_sizes : list = [],
        discard_first : int = 0,
        init_full : bool = False,
        silent : bool = False
        **kwargs) -> None:

        """
        Initializes generator.

        Args:
        -----
        time_points:
            Integer. How many total time points to generate. Not to be confused with downsampling coefficient (applied 
            later).
        nodes:
            List of strings. (Unique) node names for each node.
        node_sizes:
            List of ints. Node sizes for each node.
        discard_first:
            Integer. How many initial time points to discard. Setting higher discard_first values generally ensures 
            samples closer to equilibrium.
        init_full:
            Boolean. If True, initializes all interactions, growth rates,and initial abundances at random.
        silent:
            Boolean. If True, suppresses all print statements.
        **kwargs:
            C, d, sigma, rho for AT-Normal matrix

        Returns:
        --------
        OmicsGenerator object.
        
        Raises:
        -------
        TODO
        """

        if len(nodes) != len(node_sizes):
            raise Exception(f"Node lengths and node sizes do not match: {len(nodes)} != {len(node_sizes)}")

        self._interactions = []
        self._interventions = []
        self._time_points = time_points + discard_first
        self._T = np.array(range(self._time_points))
        self._namespace = set()
        self._discard_first = discard_first
        self._silent = silent

        # Process nodes
        self.nodes = []
        for node_name, size in zip(nodes, node_sizes):
            self.add_node(node_name, size)
        
        if init_full:
            self._init_full(**kwargs)

        if self._silent == False:
            print("Initialized")

        return None

    class _OmicsNode:
        """
        PRIVATE METHOD. Call with self.add_node() instead.

        A class for omics nodes. Contains pointers to interactions, interventions.

        Attributes:
        -----------
        inbound:
            A dict of (node name, matrix) tuples representing matrix interactions of the type Ax --> y, where y is 
            another node. Maintained by self.add_interaction().
        outbound:
            A dict of (node name, matrix) tuples representing matrix interactions of the type Ay --> x, where y is 
            another node. Maintained by self.add_interaction().
        interventions: 
            A list of interventions which affect this node. Maintained by self.add_intervention().
        """

        def __init__(
            self, 
            name : str, 
            size : int, 
            initial_value : np.ndarray,
            growth_rates : np.ndarray, 
            log_noise : bool, 
            verbose : bool = True) -> None:

            """
            Initializes node.

            Args:
            -----
            name:
                String. The node name. Must be unique.
            size:
                Integer: How many elements does this node have?
            initial_value:
                A vector of initial abundances for node elements. Length must be equal to size. Generally not called 
                on initialization - use self.add_initial_value() instead.
            growth_rates:
                Intrinsic growth/death rates for node elements. Length must be equal to size. Generally not called on 
                initialization - use self.add_initial_value() with 'growth_rate = True' instead.
            log_noise:
                Boolean. If True, noise will be added to log-relative abundances. True by default.
            verbose:
                Boolean. If False, suppresses print statements.

            Returns:
            --------
            _OmicsNode object.

            Raises:
            -------
            None (fails silently, use add_node() instead.)
            """
            
            self.name = name
            self.size = size
            self.initial_value = initial_value
            self.growth_rates = growth_rates
            self.log_noise = log_noise
            self.outbound = {}
            self.inbound = {}
            self.interventions = []

            if verbose == True:
                print(f"Node '{name}' initialized")

            return None
        
        def __str__(self):
            return f"{self.name}\t{self.size}"

    class _OmicsInteraction:

        """
        PRIVATE METHOD. Call with self.add_interaction() instead.

        A class for omics interactions. This has the general form of an m x n matrix representing interactions between 
        one set (e.g. taxa) and another set (e.g. other taxa, metabolites, whatever)

        Attributes:
        -----------
        nrows:
            Number of rows (e.g. taxa) in matrix.
        ncols:
            Number of columns (e.g. metabolites) in matrix.
        """

        def __init__(
            self, 
            name : str, 
            outbound_node : None, 
            inbound_node : None, 
            matrix : np.ndarray, 
            lag : int, 
            verbose : bool = True) -> None:

            """
            Initializes interaction.

            Args:
            -----
            name:
                String. A name for this interaction. Must be unique.
            outbound_node:
                Node from which the edge originates
            inbound_node:
                Node at which the edge terminates
            matrix:
                A matrix-type object with interactions
            lag:
                Integer. How much delay to put into dependencies. For instance, a lag of 1 on an interaction means we 
                compute Ax_t = y_(t+1)
            verbose:
                Boolean. If False, suppresses print statements.
            
            Returns:
            --------
            _OmicsInteraction object.

            Raises:
            -------
            None (fails silently, use add_interaction() instead).
            """

            self.name = name
            self.outbound_node = outbound_node
            self.inbound_node = inbound_node
            self.matrix = np.array(matrix)
            self.lag = lag
            self.nrows = matrix.shape[0] # e.g. number of taxa
            self.ncols = matrix.shape[1] # e.g. number of metabolites

            if verbose == True:
                print(f"Interaction '{name}' added")

            return None
        
        def __str__(self):
            return f"{self.name}:\t({self.outbound_node.name})-->({self.inbound_node.name})\tLag: {self.lag}"

    class _OmicsIntervention:

        """
        PRIVATE METHOD. Call with self.add_intervention() instead.

        A class for omics interventions. This has the general form of an n-length matrix which describes the reactions 
        of some set (e.g. taxa) to this particular intervention.
        """

        def __init__(
            self, 
            name : str, 
            vector : np.ndarray, 
            node_name : str, 
            U : np.ndarray, 
            affects_abundance : bool, 
            verbose : bool = True) -> None:

            """
            Initializes an intervention.

            Args:
            -----
            name:
                String. A name for our intervention. Only used for printing and other bookkeeping.
            vector:
                A vector-type object with reactions to the intervention.
            node_name:
                String. Name of node affected by this intervention/matrix.
            U:
                An indicator vector which is 1 for time points when the intervention is active, 0 otherwise.
            affects_abundance:  
                Boolean. If True, intervention vector will be applied directly to the abundance vector rather than 
                growth rates.
            verbose:
                Boolean. If False, suppresses print statements.
            
            Returns:
            --------
            _OmicsIntevention object.
            
            Raises:
            -------
            None (fails silently, use add_intervention() instead).
            """

            self.name = name
            self.vector = vector
            self.node_name = node_name
            self.U = np.array(U)
            self.affects_abundance = affects_abundance

            if verbose == True:
                print(f"Intervention '{name}' added")

            return 

        def __str__(self):
            end = ""
            if self.affects_abundance == True:
                end = "\taffects abundance"
            return f"{self.name}\t{self.node_name}{end}"

    def add_node(
        self, 
        name : str, 
        size : int, 
        initial_value : np.ndarray = None, 
        growth_rates : np.ndarray = None, 
        log_noise : bool = True, 
        verbose : bool = True) -> None:

        """
        Adds nodes to generator object.

        Args:
        -----
        name:           
            String. Used to identify node. Must be unique.
        size:           
            Length of vector associated with a time point of this node. For instance, for a metagenomics node, this 
            would correspond to the number of taxa.
        initial_value:  
            Value of this node at t = 0. Must be same length as node size.
        growth_rates:   
            Element-wise growth/death rates for this node. Must be same length as node size.
        log_noise:      
            Boolean. If True, noise will be added to log-relative abundance.If False, noise will be added to relative 
            abundances.
        verbose:
            Boolean. If False, suppresses print statements.
    
        Returns:
        --------
        None (modifies generator in place).

        Raises:
        -------
        TODO
        """

        # Check sizes
        if initial_value != None:
            if len(initial_value) != size:
                raise Exception(f"initial_value is wrong size: {len(initial_value)} != {size}")
        if growth_rates != None:
            if len(growth_rates) != size:
                raise Exception(f"growth_rates is wrong size: {len(growth_rates)} != {size}")
        
        # Check namespace
        if name in self._namespace:
            raise Exception(f"Name {name} already in use. Please use a unique name")

        # Check verbosity
        if self._silent == True:
            verbose = False

        # Generate node and append to object
        node = self._OmicsNode(
            name,
            size,
            initial_value, 
            growth_rates, 
            log_noise, 
            verbose
        )
        self.nodes.append(node)
        self._namespace.add(name)

        return None
    
    def add_interaction(
        self, 
        name : str, 
        outbound_node_name : str, 
        inbound_node_name : str, 
        matrix : np.ndarray, 
        lag : int = 0, 
        verbose : bool = True) -> None:

        """
        Adds interactions to generator object.

        Edges look like this:
            Graphical:          (OUTBOUND NODE)--->(INBOUND NODE)
            Linear algebra:     [inbound] = [matrix] @ [outbound] + [...]

        Args:
        -----
        name:               
            String. A name for this interaction.
        outbound_node_name: 
            String. Name of node from which the edge originates
        inbound_node_name:  
            String. Name of node at which the edge terminates
        matrix:             
            A matrix-type object with interactions
        lag:                
            Integer. How much delay to put into dependencies. For instance, a lag of 1 on an interaction means we 
            compute Ax_t = y_(t+1)
        verbose:
            Boolean. If False, suppresses print statements.

        Returns:
        --------
        None (modifies generator in place).

        Raises:
        -------
        TODO
        """

        # Check namespace
        if name in self._namespace:
            raise Exception(f"Name {name} already in use. Please use a unique name")

        # Check verbosity
        if self._silent == True:
            verbose = False

        # Get nodes
        outbound_node = self.get(outbound_node_name, "node")
        if outbound_node == None:
            raise Exception("Outbound node is invalid")

        inbound_node = self.get(inbound_node_name, "node")
        if inbound_node == None:
            raise Exception("Inbound node is invalid")

        # Check that matrix dimensions match
        if matrix.shape[1] != inbound_node.size:
            raise Exception(f"Matrix dimension mismatch: matrix shape[1] = {matrix.shape[1]} != {inbound_node.size} (size of inbound node '{inbound_node.name}')")
        if matrix.shape[0] != outbound_node.size:
            raise Exception(f"Matrix dimension mismatch: matrix shape[0] = {matrix.shape[0]} != {outbound_node.size} (size of outbound node '{outbound_node.name}')")

        interaction = self._OmicsInteraction(
            name, 
            outbound_node, 
            inbound_node, 
            matrix, 
            lag, 
            verbose
        )
        self._interactions.append(interaction)

        # Append to nodes
        outbound_node.inbound[inbound_node_name]  = interaction
        inbound_node.outbound[outbound_node_name] = interaction

        self._namespace.add(name)

        return None

    def add_intervention(
        self, 
        name : str, 
        node_name : str, 
        vector : np.ndarray, 
        affects_abundance : bool = False, 
        U : np.ndarray = None, 
        start : int = None, 
        end : int = None, 
        verbose : bool = True) -> None:

        """
        Adds an intervention to generator.

        Must have either U or (start, end) set to specify timeframe.

        Args:
        -----
        name:               
            String. A name for our intervention. Only used for printing and other bookkeeping.
        node_name:          
            String. Name of node affected by this intervention/matrix.
        vector:             
            A vector-type object detailing, elementwise, the reactions of each node coordinate to an intervention.
        affects_abundance:  
            Boolean. If True, intervention vector will be applied directly to the abundance vector rather than to growth 
            rates.
        U:                  
            An indicator vector which is 1 for time pointswhen the intervention is active, 0 otherwise.
        start:              
            First time point when interaction begins. Use only for interactions of the form 0*1+0*. Otherwise, use U 
            variable instead.
        end:                
            Last node when interaction is active. Use only for interactions of the form 0*1+0*. Otherwise, use U 
            variable instaed.
        verbose:
            Boolean. If False, suppresses print statements.
            
        Returns:
        --------
        None (modifies generator in place).

        Raises:
        -------
        TODO
        """

        # Check namespace
        if name in self._namespace:
            raise Exception(f"Name {name} already in use. Please use a unique name")

        # Check U vector is correct length
        if U != None:
            if len(U) != self._time_points:
                raise Exception(f"U vector is different size from number of time points: {len(U)} != {self._time_points}")

        # Check verbosity
        if self._silent == True:
            verbose = False

        # Process node
        node = self.get(node_name, "node")
        if node == None:
            raise Exception("Invalid node! Please try again")

        # A bunch of control flow to make a boolean vector called U
        if U == None and (start == None or end == None):
            raise Exception("Need to supply a (start,end) pair or a U vector") 
        elif start != None and end != None:
            U = np.array([0] * self._time_points)
            U[start:end] = 1
        else:
            raise Exception("Other exception in OmicsIntervention")
        
        # Make the intervention and add it to self
        intervention = self._OmicsIntervention(
            name, 
            vector, 
            node_name, 
            U, 
            affects_abundance, 
            verbose=False
        )
        if len(intervention.U) == self._time_points:
            self._interventions.append(intervention)
        else:
            raise Exception("Intervention vector is not the same length at time vector")

        # Modify node accordingly
        node.interventions.append(intervention)

        self._namespace.add(name)

        return None

    def set_initial_value(
        self, 
        node_name : str, 
        values : np.ndarray, 
        growth_rate : bool = False, 
        verbose : bool = True) -> None:

        """
        Sets a node value or growth rate.

        Args:
        -----
        node_name:      
            Name of node being altered
        values:         
            Vector. Initial values for node. Must be same length as node size.
        growth_rate:    
            Boolean. If True, affects the growth_rate parameter of the node. Otherwise, affects initial values of node.
        verbose:
            Boolean. If False, suppresses print statements.

        Returns:
        --------
        None (modifies generator in place).

        Raises:
        -------
        TODO
        """

        node = self.get(node_name, "node")

        # Check node exists
        if node == None:
            raise Exception(f"Invalid node name: {node_name} does not exist")

        # Check dimensions match
        if len(values) != node.size:
            raise Exception(f"Size mismatch with node size: {len(values)} != {node.size}")
        
        # Set values
        if growth_rate == False:
            node.initial_value = values
        elif growth_rate == True:
            node.growth_rates = values
        
        # Print output
        if verbose == True and self._silent == False:
            print(f"Added growth rates to node {node_name}")
        
        return None
        
    def get(
        self, 
        name : str, 
        type : str in ["node", "interaction", "intervention"] = None) -> "generator element":

        """
        Gets a (node/interaction/intervention) by name.

        Args:
        -----
        name:
            String. Name of node/interaction/intervention.
        type:
            String. One of ["node", "interaction", "intervention"]. Specifies the type of generator element to look for.

        Returns:
        --------
        _OmicsNode, _OmicsInteraction, _OmicsIntervention, or None.

        Raises:
        -------
        None
        """

        if type == None or type == "node":
            for node in self.nodes:
                if node.name == name:
                    return node

        if type == None or type == "interaction":
            for interaction in self._interactions:
                if interaction.name == name:
                    return interaction

        if type == None or type == "intervetion":
            for intervention in self._interventions:
                if intervention.name == name:
                    return intervention

        return None

    def remove(
        self, 
        name : str, 
        verbose : bool = True) -> None:

        """
        Removes a node, intervention, or interaction from the generator by name.

        Args:
        -----
        name:   
            A string specifying the (unique) name of the element to be removed.
        
        Returns:
        --------
        None (modifies generator in place).
        
        Raises:
        -------
        TODO
        """

        obj = self.get(name)

        if obj is None:
            raise Exception(f"Cannot find object named {name} to remove")

        if type(obj) is self._OmicsNode:
            for interaction in reversed(self._interactions): # reversed so we can remove interactions as we go
                if interaction.inbound_node == obj or interaction.outbound_node == obj:
                    self._interactions.remove(interaction)
            for intervention in reversed(self._interventions): 
                if intervention.node_name == name:
                    self._interventions.remove(intervention)
            for node in self.nodes:
                node.inbound.pop(name, None)
                node.outbound.pop(name, None)
            self.nodes.remove(obj)
            if verbose == True:
                print(f"Removed node '{name}'")

        elif type(obj) is self._OmicsInteraction:
            # Remove interaction from inbound node
            obj.inbound_node.outbound.pop(obj.outbound_node.name, None)

            # Remove interaction from outbound node
            obj.outbound_node.inbound.pop(obj.inbound_node.name, None)  

            # Remove interaction from list
            self._interactions.remove(obj)

            if verbose == True:
                print(f"Removed interaction '{name}'")

        elif type(obj) is self._OmicsIntervention:
            node = self.get(obj.node_name)
            node.interventions.remove(obj)
            self._interventions.remove(obj)
            if verbose == True:
                print(f"Removed intervention '{name}'")

        else:
            raise Exception(f"Cannot remove '{name}': unknown type. Is the name correct?")

        self._namespace.remove(name)

        return None

    def generate(
        self, 
        noise_var : float = 1e-2, 
        n_reads : int = 1e5, 
        dt : float = 1e-2, 
        downsample : int = 1) -> (dict, dict, dict):

        """
        Generates a single timecourse of synthetic data. 
        
        Args:
        -----
        noise_var:   
            Float. variance parameter for gaussian noise term.
        n_reads:
            Integer. Number of reads to draw from the unsampled distribution.
        dt:         
            Float. time step size which gets passed to IVP solver
        downsample: 
            Integer. fraction of outputs to keep (1/n). By default, keeps all samples. downsample=4 means every 4th 
            sample is kept, etc. Downsample is deprecated. Simply modify "dt" instead.

        Returns:
        --------
        The following three dicts (in order):

            //======================================================\\
            ||Name:   Sampling:   Normalization:  Number of samples:||
            ||======================================================||
            ||X       unsampled   normalized      downsampled       ||
            ||Y       sampled     normalized      downsampled       ||
            ||Z       unsampled   unnormalized    full              ||
            \\======================================================//
        
        Each X/Y/Z dict contains (node, timecourse) pairs. The timecourse is a numpy array with shape (number of time 
        points, node size).

        Raises:
        -------
        TODO
        """

        def _grad_fn(
            node : None, 
            X : list, 
            growth_rates : np.ndarray, 
            t : int) -> None:

            """
            This gets passed to the solver. It's just the vector f used in GLV calculations.
            """

            # Interactions:
            interaction_coef = np.zeros(node.size)
            for node_name in node.outbound:
                interaction = node.outbound[node_name]

                # Adjust for lag
                idx = -1 - interaction.lag
                try:
                    # Get interaction matrix
                    M = interaction.matrix

                    # Get last value (modulo lag term) of node abundance
                    y = X[node_name][idx]
                    
                     # f += yM (GLV equation)
                    interaction_coef += y @ M

                except IndexError:
                    # Happens when lag is larger than number of values already generated
                    pass
            
            # Interventions:
            intervention_coef = np.zeros(node.size)
            for intervention in node.interventions:
                if intervention.affects_abundance == False:
                    intervention_coef += intervention.vector.dot(intervention.U[t])
            
            # Self
            xt = X[node.name][-1]

            # The function itself:
            def fn(t, x):
                return xt * (growth_rates + interaction_coef + intervention_coef)
            return fn

        # Initialization steps
        X = {} # Abundances
        Y = {} # Final outputs
        Z = {} # Latent absolute abundances

        for node in self.nodes:
            X[node.name] = [node.initial_value]
 
        # Generalized Lotka-Volterra steps, plus bells and whistles
        for t in range(self._time_points - 1):
            X_temp = {} # Use this so that all values are updated at once

            for node in self.nodes:
                # Get values from dicts
                x = X[node.name]
                g = node.growth_rates

                # Initialize values
                Xprev = np.copy(x[-1])  # last time point, X_(t-1)

                # Pass to solver
                # TODO: possible to do this all in one shot rather than looping?
                grad = _grad_fn(node, X, g, t)
                ivp = solve_ivp(grad, (0,dt), Xprev, method="RK45")
                Xt = ivp.y[:,-1]

                # Tweak abundances on a per-node basis
                # TODO: Maybe this would be better if it were size-adjusted?
                for intervention in node.interventions:
                    if intervention.affects_abundance == True:
                        Xt += intervention.vector * intervention.U[t]

                # Add biological noise:
                noise = np.random.normal(scale=noise_var, size=node.size) 

                # No noise for missing taxa
                noise = noise * (Xt > 0)

                # Equivalent to log->add noise->exp
                if node.log_noise == True:
                    Xt *= np.exp(noise)                                   
                else:
                    Xt += noise

                # Push to results
                Xt = np.clip( Xt, 0, None )
                X_temp[node.name] = Xt
            
            # Push all values for this time point to X at once
            for key in X_temp:
                X[key] += [X_temp[key]]

        # Simulate sampling noise
        for node in self.nodes:
            x = np.array(X[node.name])

            # Save latent state
            x_old = x.copy()

            # Discard first couple elements (ensure values are near attractor)
            x = x[self._discard_first:]

            # Take every nth element
            # Negative coefficient ensures we sample from the end
            x = x[::-downsample]

            # Need to un-reverse the data now
            x = x[::-1]

            # Relative abundances         
            x = np.apply_along_axis(lambda a: a/sum(a), 1, x)

            # Draw samples
            y = []
            for idx in range(x.shape[0]):
                try:
                    Yt = np.random.multinomial(n_reads, x[idx]) / n_reads
                    y += [Yt]
                except ValueError:
                    # TODO: circle back and figure out what was breaking this
                    # print("ERROR: check self._weird for more info")
                    # self._weird = X[node.name][idx] # debugging variable
                    y += [np.zeros(node.size)]

            # Push to output
            X[node.name] = x
            Y[node.name] = np.array(y)
            Z[node.name] = x_old
        
        return X, Y, Z

    def generate_multiple(
        self, 
        n : int, 
        extinct_fraction : float = 0, 
        **generate_args) -> (list, list, list):
        
        """
        Generates several timecourses of synthetic data. 

        This is essentially a wrapper around a loop of generate() calls, with the added element of reinitializing 
        individuals. The extinct_fraction parameter gives some degree of control over re-initialization.
        
        Args:
        -----
        n:
            Integer. Number of individuals for whom to generate synthetic data timecourses.
        extinct_fraction:
            Float in [0, 1) range. Fraction of abundances that should be extinct for each individual.

        Additional args (same as generate()):
        -------------------------------------
        noise_var:   
            Float. variance parameter for gaussian noise term.
        n_reads:
            Integer. Number of reads to draw from the unsampled distribution.
        dt:         
            Float. time step size which gets passed to IVP solver
        downsample: 
            Integer. fraction of outputs to keep (1/n). By default, keeps all samples. downsample=4 means every 4th 
            sample is kept, etc. Downsample is deprecated. Simply modify "dt" instead.

        Returns:
        --------
        The following three arrays (in order):
        
            //======================================================\\
            ||Name:   Sampling:   Normalization:  Number of samples:||
            ||======================================================||
            ||X       unsampled   normalized      downsampled       ||
            ||Y       sampled     normalized      downsampled       ||
            ||Z       unsampled   unnormalized    full              ||
            \\======================================================//
        
        Each X/Y/Z array contains n dicts, each of which contains (node, timecourse) pairs. The timecourse is a numpy 
        array with shape (number of time points, node size).

        Raises:
        -------
        TODO
        """

        # Initialize:
        old_nodes = self.nodes # store old initial values
        out_X = []
        out_Y = []
        out_Z = []

        # Generation loop
        for i in range(n):
            # Set new initial values for each node
            for node in self.nodes:
                # TODO: allow passing of any function to generate this
                abundances = np.random.exponential(size=node.size) * np.random.binomial(1, 1-extinct_fraction, size=node.size)
                self.set_initial_value(node.name, abundances, verbose=False)
            
            X,Y,Z = self.generate(**generate_args)
            out_X.append(X)
            out_Y.append(Y)
            out_Z.append(Z)
        
        # return nodes to old values
        self.nodes = old_nodes

        return out_X, out_Y, out_Z

    def _allesina_tang_normal_matrix(
        self, 
        n : int, 
        C : float, 
        d : float, 
        sigma : float, 
        rho : float) -> np.ndarray:

        """
        Generates an Allesina-Tang normal matrix.

        Inspired by https://stefanoallesina.github.io/Sao_Paulo_School/intro.html#multi-species-dynamics.

        How this works:
        ---------------
        1. Creates covariance matrix has the following form:
            1   rho rho  ...
            rho 1   rho  ...
            rho rho 1    ...
            ... (you get the idea)
        2. Draws multivariate normal pairs from this covariance matrix
        3. Populates non-diagonal entries of matrix with drawn pairs
        4. Symmetrically sparsifies matrix, keeping only ~C% of entries
        5. Sets diagonals of matrix to -d

        Args:
        -----
        n:
            Integer. Number of rows/columns in square matrix.
        C:
            Float in (0,1]: Sparsity parameter. Higher C = less sparse.
        d:
            Float. Negative self-interaction size.
        sigma:
            Float. Variance used to generate multivariate normal covariance matrix.
        rho: 
            Float in [-1, 1]. Correlation term of covariance matrix. Higher rho = positive connectance = mutualism = 
            harder to stabilize. Lower rho = predator-prey--type relationships = easier to stabilize.
        
        Returns:
        --------
        A matrix M that can be used as an interaction matrix.
        
        Raises:
        -------
        None (fails silently).
        """

        # sample coefficients
        mu = np.zeros(2)
        cov = sigma ** 2 * np.array([[1, rho], [rho, 1]])
        n_samples = int(n * (n-1) / 2)
        pairs = np.random.multivariate_normal(mu, cov, n_samples)

        # completely filled matrix
        M = np.ndarray((n, n))
        M[np.triu_indices(n, 1)] = pairs[:,0]
        M = M.transpose()
        M[np.triu_indices(n, 1)] = pairs[:,1]

        # winnow down
        connections = np.random.rand(n, n) <= C
        connections = connections * 1 # binarize
        connections[np.tril_indices(n,1)] = 0
        connections += connections.transpose() # symmetric
        M *= connections

        # set negative self-interactions
        M[np.diag_indices(n)] = -d

        return M
    
    def _set_interactions(
        self,
        C : float = 0.5, 
        d : float = 5, 
        sigma : float = 1, 
        rho : float = -0.4) -> None:

        """
        Sets all interaction matrices from one big AT-normal matrix

        Args:
        -----
        C:
            Float in (0,1]: Sparsity parameter. Higher C = less sparse.
        d:
            Float. Negative self-interaction size.
        sigma:
            Float. Variance used to generate multivariate normal covariance matrix.
        rho: 
            Float in [-1, 1]. Correlation term of covariance matrix. Higher rho = positive connectance = mutualism = 
            harder to stabilize. Lower rho = predator-prey--type relationships = easier to stabilize.
        
        Returns:
        --------
        None (modifies generator in place).
        
        Raises:
        -------
        None (fails silently).
        """

        # Generate master matrix
        sizes = [node.size for node in self.nodes]
        n = np.sum(sizes)
        n_nodes = len(sizes)
        m0 = self._allesina_tang_normal_matrix(n, C, d, sigma, rho)

        # Carve up master matrix
        i = 0 # row
        for node1 in self.nodes:
            j = 0 # col
            for node2 in self.nodes:
                m_ij = m0[i:i + node1.size, j:j + node2.size]
                self.add_interaction(
                    f"{node1.name}->{node2.name}", 
                    node1.name, 
                    node2.name, 
                    m_ij
                )

                if self._silent == False:
                    print(f"set m:({node1.name})->({node2.name}):   {i}:{i + node1.size}    {j}:{j + node2.size}")

                j += node2.size
            i += node1.size

        return None
    
    def _init_full(
        dist : None = np.random.exponential,
        **kwargs) -> None:

        """
        A fully random initialization of all generator parameters.

        Args:
        -----
        dist:
            A function to draw initial distributions (e.g. np.random.exponential, np.random.rand, etc)

        Returns:
        --------
        None (modifies generator in place)

        Raises:
        -------
        None
        """

        self._set_interactions(**kwargs)
        for node in self.nodes:
            self.set_initial_value(
                node.name, 
                np.random.exponential(size=node.size)
            )
            self.set_initial_value(
                node.name, 
                2 * (0.5 - np.random.rand(node.size)),
                growth_rate=True
            )
        
        return None

    def case_control(
        self,
        participants : int, 
        case_frac : float,
        node_name: str,
        effect_size : float, 
        **generate_args) -> (list, list, list, list, list, list):
        """
        Generates synthetic case and control timecourses.

        Args:
        -----
        participants:
            Integer. The total number of participants in the study.
        case_frac:
            Float in [0,1]. Fraction of total participants belonging to the case group.
        node_name:
            String. Name of node to which the intervention is applied.
        effect_size:
            Float. Magnitude of intervention.
        **kwargs:
            Arguments that get passed to generate_multiple().
        
        Returns:
        --------
        X_control:
            X-list like generate_multiple() for control group.
        Y_control:
            Y-list like generate_multiple() for control group.
        Z_control:
            Z-list like generate_multiple() for control group.
        X_case:
            X-list like generate_multiple() for case group.
        Y_case:
            Y-list like generate_multiple() for case group.
        Z_case:
            Z-list like generate_multiple() for case group.
        
        Raises:
        -------
        TODO
        """

        # inferred settings
        n_cases = int(participants * case_frac)
        n_controls = int(participants * (1-case_frac))

        # get control values
        x_control, y_control, z_control = self.generate_multiple(n_controls, **generate_args)

        # get case values
        case_gen = self.copy()
        node_size = self.get(node_name).size
        case_gen.add_intervention(
            'CASE', 
            node_name, 
            effect_size * np.random.rand(node_size), 
            start=0, 
            end=self._time_points
        )
        x_case, y_case, z_case = case_gen.generate_multiple(n_cases, **generate_args)

        return x_control, y_control, z_control, x_case, y_case, z_case

    def copy(self) -> None:

        """
        Makes a deep copy of generator.

        Args:
        -----
        None

        Returns:
        --------
        OmicsGenerator copy

        Raises:
        -------
        None
        """

        return deepcopy(self)

    def __str__(self):
        # TODO: Rewrite this more cleanly with f-strings
        out = "\n=========================GENERATOR=========================\n\nTime_points:\t"
        out += str(self._time_points)
        out += "\n\nNodes:\n\t"
        out += "\n\t".join([ str(x) for x in self.nodes ] )
        out += "\n\nInteractions:\n\t"
        out += "\n\t".join([ str(x) for x in self._interactions ] )
        out += "\n\nInterventions:\n\t"
        out += "\n\t".join([ str(x) for x in self._interventions ] )
        return out