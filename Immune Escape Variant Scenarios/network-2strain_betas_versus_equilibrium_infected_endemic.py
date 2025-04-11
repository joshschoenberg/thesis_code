#For multiple beta values:

import numpy as np
import random
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from numba import njit
from copy import deepcopy
import networkx as nx
import scipy.stats as stats
from multiprocessing import Pool
import seaborn as sns
from itertools import chain


def create_graph(connection_bias, num_nodes=600, hesitant_percentage=0.3):
    """Create a graph with specified connection bias and assortativity."""
    hesitant_count = int(num_nodes * hesitant_percentage)
    non_hesitant_count = num_nodes - hesitant_count
    types = np.array(["hesitant"] * hesitant_count + ["non_hesitant"] * non_hesitant_count)
    np.random.shuffle(types)

    p_diff = (10 / connection_bias) / num_nodes
    p_hh = (10 - non_hesitant_count * p_diff) / hesitant_count
    p_nn = (10 - hesitant_count * p_diff) / non_hesitant_count

    # Assuming `types` is a NumPy array for efficient comparison
    types = np.array(types)  # Ensure types is a NumPy array

    # Create a random matrix of probabilities
    rand_probs = np.random.rand(num_nodes, num_nodes)

    # Create a boolean mask for pairs (i, j) where i < j (upper triangular part)
    upper_tri_mask = np.triu(np.ones((num_nodes, num_nodes), dtype=bool), k=1)

    # Pairwise type comparison
    same_type = types[:, None] == types

    # Assign probabilities based on type conditions
    prob_matrix = np.full((num_nodes, num_nodes), p_diff)  # Default to `p_diff`
    prob_matrix[(types[:, None] == "hesitant") & same_type] = p_hh
    prob_matrix[(types[:, None] == "non_hesitant") & same_type] = p_nn

    # Generate adjacency matrix
    adj_matrix = (rand_probs < prob_matrix) & upper_tri_mask  # Apply probability threshold

    # Make the matrix symmetric
    adj_matrix = adj_matrix | adj_matrix.T

    # Convert to a NetworkX graph
    G = nx.from_numpy_array(adj_matrix.astype(int))
    for node in G.nodes():
        G.nodes[node]['type'] = 1 if types[node] == "hesitant" else 0

    assortativity = nx.attribute_assortativity_coefficient(G, 'type')
    return G, assortativity

def simulate_one_run(params):
    num_nodes = params['N']
    connection_bias = params['connection_bias']
    initial_infected_number = params['initial_infected_number']
    dt = params['dt']
    time = params['time']
    time_steps = time/dt
    gamma = params['gamma']
    avg_connectivity = params['avg_connectivity']
    beta = params['beta']
    vaccine_effectiveness = params['sigma']
    frac_vaccpos = params['frac_vaccpos']
    hesitant_percentage = 1 - frac_vaccpos
    p_vax = params['p_vax']
    p_mutation = params['p_mutation']
    contacts_per_infected = params['contacts_per_infected']
    var2_outbreak_threshold = params['var2_outbreak_threshold']
    var2_outbreak_time_threshold = params['var2_outbreak_time_threshold']
    var2_outbreaks = 0

    S_t = []
    I1_t = []
    I2_t = []
    R1_t = []
    R2_t = []
    times = []

    G, assortativity = create_graph(connection_bias, num_nodes, hesitant_percentage=hesitant_percentage)

    S = np.ones(num_nodes)
    I1 = np.zeros(num_nodes)  # First strain
    I2 = np.zeros(num_nodes)  # Mutated strain
    R1 = np.zeros(num_nodes)
    R2 = np.zeros(num_nodes)
    V = np.zeros(num_nodes)
  
    # initial_infected_nodes = np.random.choice(np.arange(num_nodes), initial_infected_number, replace=False)
    # I1[initial_infected_nodes] = 1
    # S[initial_infected_nodes] = 0

    ## Vaccinate people before simulation starts
    susceptible_nodes = np.where(S == 1)[0]
    if len(susceptible_nodes) > 0:
        for i in susceptible_nodes:
            if G.nodes[i]['type'] == 0 and np.random.rand() < p_vax:
                S[i] = 0
                V[i] = 1
    
    # Seed exactly one initially infected person from the hesitant group
    hesitant_indices = np.where((V == 0) & (S == 1))[0]
    if len(hesitant_indices) > 0:
        initial_infected = np.random.choice(hesitant_indices)
        I1[initial_infected] = 1
        S[initial_infected] = 0

    # -------------------------
    # Simulation Loop
    # -------------------------
    t = 0
    while t < time_steps and np.sum(I1 + I2) > 0:
        infected_nodes_var_1 = np.where(I1 == 1)[0]
        infected_nodes_var_2 = np.where(I2 == 1)[0]
        all_infected_nodes = np.concatenate((infected_nodes_var_1, infected_nodes_var_2))
        infected_nodes_random_order = np.random.permutation(all_infected_nodes)
        for i in infected_nodes_random_order:
            if I1[i] == 1:
                ## Choose infected person's contacts
                neighbors = list(G.neighbors(i))  # Get the neighbors of i
                if len(neighbors) > contacts_per_infected:
                    contacts = np.random.choice(neighbors, contacts_per_infected, replace=False)  # Choose random contacts
                else:
                    contacts = neighbors  # If not enough, take all
                ## Infect the contacts
                for contact in contacts:
                    if S[contact] == 1:
                        if random.random() < beta * dt: 
                            # Risk of mutation to strain 2
                            if np.random.rand() < p_mutation:
                                I2[contact] = 1
                            else:
                                I1[contact] = 1
                            S[contact] = 0
                    elif (R1[contact] == 1 or V[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                        effective_beta = (1.0 - vaccine_effectiveness) * beta
                        if random.random() < effective_beta * dt:
                            R1[contact] = 0
                            if random.random() < p_mutation:
                                I2[contact] = 1
                            else:
                                I1[contact] = 1
                          
            elif I2[i] == 1:
                neighbors = list(G.neighbors(i))  # Get the neighbors of i
                if len(neighbors) > contacts_per_infected:
                    contacts = np.random.choice(neighbors, contacts_per_infected, replace=False)  # Choose random contacts
                else:
                    contacts = neighbors  # If not enough, take all
                ## Infect the contacts
                for contact in contacts:
                    # ## 1. Vaccine Universal, recovered different for the two strains
                    if (S[contact] == 1):
                        if random.random() < beta * dt:
                            I2[contact] = 1
                            S[contact] = 0
                    elif (R1[contact] == 1 and V[contact] == 0) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                        # Recovered/protected; if just recovered from strain 1, transmission probability is NOT reduced.
                        effective_beta = beta 
                        if random.random() < effective_beta * dt:
                            R1[contact] = 0
                            I2[contact] = 1
                    elif (V[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                        # If vaccinated, transmission probability IS reduced.
                        effective_beta = beta * (1.0 - vaccine_effectiveness)
                        if random.random() < effective_beta * dt:
                            R1[contact] = 0
                            I2[contact] = 1

                    ## 2. Vaccine not universal, recovered different for the two strains
                    # if (S[contact] == 1):
                    #     if random.random() < beta * dt:
                    #         I2[contact] = 1
                    #         S[contact] = 0
                    # elif (R1[contact] == 1 or V[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # Recovered/protected; if just recovered from strain 1, transmission probability is NOT reduced.
                    #     effective_beta = beta
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1

                    # ## 3. Vaccine universal, recovered same for the two strains
                    # if (S[contact] == 1):
                    #     if random.random() < beta * dt:
                    #         I2[contact] = 1
                    #         S[contact] = 0

                    # elif (V[contact] == 1 or R1[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # If vaccinated, transmission probability IS reduced.
                    #     effective_beta = beta * (1.0 - vaccine_effectiveness)
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1

                    ## 4. Vaccine not universal, recovered same for the two strains
                    # if (S[contact] == 1):
                    #     if random.random() < beta * dt:
                    #         I2[contact] = 1
                    #         S[contact] = 0
                    # elif (V[contact] == 1 and R1[contact] == 0) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # Vaccinated; if just vaccinated , transmission probability is NOT reduced.
                    #     effective_beta = beta 
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1
                    # elif (R1[contact] == 1) and (R2[contact] + I2[contact] + I1[contact]) == 0:
                    #     # If recpvered, transmission probability IS reduced.
                    #     effective_beta = beta * (1.0 - vaccine_effectiveness)
                    #     if random.random() < effective_beta * dt:
                    #         R1[contact] = 0
                    #         I2[contact] = 1



        ## Recoveries
        for i in infected_nodes_var_1:
            if random.random() < gamma * dt:
                I1[i] = 0
                R1[i] = 1
        for i in infected_nodes_var_2:
            if random.random() < gamma * dt:
                I2[i] = 0
                R2[i] = 1
        
        if np.sum(I2) > var2_outbreak_threshold:
            if times[-1] >= var2_outbreak_time_threshold: # Only if transient has passed!
                var2_outbreaks += 1
            infected_indices_2 = np.where(I2 == 1)[0]
            for i in infected_indices_2:
                # Move them to recovered from strain 1 so as to
                # disturb the dynamics as little as possible ...
                S[i] = 0
                I1[i] = 0
                R1[i] = 1
                I2[i] = 0
                R2[i] = 0 
            recovered_indices_2 = np.where(R2 == 1)[0]
            for i in recovered_indices_2:
                # Move them to recovered from strain 1 so as to
                # disturb the dynamics as little as possible ...
                S[i] = 0
                I1[i] = 0
                R1[i] = 1
                I2[i] = 0
                R2[i] = 0   
        I1_t.append(np.sum(I1))
        I2_t.append(np.sum(I2))
        R1_t.append(np.sum(R1))
        R2_t.append(np.sum(R2))
        S_t.append(np.sum(S))
        times.append(t*dt)     
        t += 1

    total_infected = np.sum(R1) + np.sum(R2) + np.sum(I1) + np.sum(I2)
    total_infected_var_1 = np.sum(R1) + np.sum(I1)
    total_infected_var_2 = np.sum(R2) + np.sum(I2)
    # time_of_peak_infection = np.argmax(i_over_time)
    percent_considered_outbreak = 0.01
    # outbreak_occurred_var_1 = total_infected_var_1 > (percent_considered_outbreak * num_nodes)
    # outbreak_occurred_var_2 = total_infected_var_2 > (percent_considered_outbreak * num_nodes)
    # outbreak_occurred = total_infected > (percent_considered_outbreak * num_nodes)

    times = np.array(times)
    peak_infected = max(I1_t + I2_t)
    peak_infected_var_1 = max(I1_t)  # Peak infected for strain 1 only
    outbreak = peak_infected >= num_nodes * percent_considered_outbreak
    # print(var2_outbreaks)
    I_eq = I1_t[-1]
    return times, I1_t, I2_t, peak_infected, outbreak, var2_outbreaks, I_eq, peak_infected_var_1


def run_simulations_for_beta(params, num_simulations=200, max_workers=20):
    """
    Runs multiple simulation realizations in parallel for a given beta.
    Returns:
        results (list of tuples): Each tuple is (peak_infected, outbreak).
    """
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(simulate_one_run, params) for _ in range(num_simulations)]
        for future in as_completed(futures):
            results.append(future.result())
    return results

def main():

    # Define a structured array type for parameters
    param_dtype = np.dtype([
        ('beta', np.float64),
        ('connection_bias', np.float64),
        ('sigma', np.float64),
        ('gamma', np.float64),
        ('N', np.int64),
        ('p_mutation', np.float64),
        ('initial_infected_number', np.int64),
        ('dt', np.float64),
        ('time', np.float64),
        ('avg_connectivity', np.float64),
        ('frac_vaccpos', np.float64),
        ('p_vax', np.float64),
        ('contacts_per_infected', np.int64),
        ('var2_outbreak_threshold', np.float64),
        ('var2_outbreak_time_threshold', np.float64),
    ])

    params = np.zeros(1, dtype=param_dtype)[0]

    params["beta"] = 0.8
    params["gamma"] = 0.2
    params["N"] = 3000
    params["p_mutation"] = 0
    # params["p_mutation"] = 0.05

    params["initial_infected_number"] = 1
    params["dt"] = 0.5
    params["time"] = 50
    params["avg_connectivity"] = 10
    params["p_vax"] = 1
    params["frac_vaccpos"] = 0.7
    params["contacts_per_infected"] = 1
    params["var2_outbreak_threshold"] = 0.01 * params["N"]
    params["var2_outbreak_time_threshold"] = params["time"] * params["dt"]
    params["sigma"] = 0.5

    n_sims = 300

    params["connection_bias"] = 1
    a_1 = 0.00

    outbreak_probabilities_1 = []
    equilibrium_infected_1 = []
    peak_infected_1 = []
    betas = np.linspace(0.1, 2.5, 15)
    # Loop over each sigma value.
    for beta in betas:
        params_loc = deepcopy(params)
        params_loc['beta'] = beta
        results = run_simulations_for_beta(params_loc, num_simulations=n_sims, max_workers=8)
        #outbreak_flags = [res[1] for res in results]
        #fraction_outbreak = sum(outbreak_flags) / len(outbreak_flags)
        #outbreak_probabilities_1.append(fraction_outbreak)
        n_var2_outbreaks = [res[5] for res in results]
        # Calculate number of infected at equilibrium 
        eq_infected = [res[6] for res in results]
        peak_infected_var_1 = [res[7] for res in results]
        peak_infected_1.append(np.mean(peak_infected_var_1))  # Store the average peak infected for strain 1
        equilibrium_infected_1.append(np.mean(eq_infected))
        outbreak_probabilities_1.append(np.mean(n_var2_outbreaks))

        print(f"beta = {beta:.2f} -> Outbreak rate: {np.mean(n_var2_outbreaks):.3f}, equilibrium infected: {np.mean(eq_infected):.3f}, peak infected: {np.mean(peak_infected_var_1):.3f}")

        # -------------------------
    # Plotting the Results 1 # recall that the sigma parameter controls the strength of immunity (vaccinated/recovered)
    # -------------------------
    params["connection_bias"] = 90
    a_2 = 0.99
    outbreak_probabilities_2 = []
    equilibrium_infected_2 = []
    peak_infected_2 = []
    
    # Loop over each sigma value.
    for beta in betas:
        params_loc = deepcopy(params)
        params_loc['beta'] = beta
        results = run_simulations_for_beta(params_loc, num_simulations=n_sims, max_workers=25)
        #outbreak_flags = [res[1] for res in results]
        #fraction_outbreak = sum(outbreak_flags) / len(outbreak_flags)
        #outbreak_probabilities_1.append(fraction_outbreak)
        # n_var2_outbreaks = [res[5] for res in results]
        eq_infected = [res[6] for res in results]
        peak_infected_var_1 = [res[7] for res in results]
        peak_infected_2.append(np.mean(peak_infected_var_1))  # Store the average peak infected for strain 1
        equilibrium_infected_2.append(np.mean(eq_infected))
        outbreak_probabilities_2.append(np.mean(n_var2_outbreaks))
        print(f"beta = {beta:.2f} -> Outbreak rate: {np.mean(n_var2_outbreaks):.3f}, equilibrium infected: {np.mean(eq_infected):.3f}, peak infected: {np.mean(peak_infected_var_1):.3f}")
    

    plt.figure(figsize=(8, 5))
    plt.plot(np.array(betas), equilibrium_infected_1, marker='o', linestyle='-', linewidth=4, label=f"Assortativity low ({a_1})")
    plt.plot(np.array(betas), equilibrium_infected_2, marker='o', linestyle='-', linewidth=4,label=f"Assortativity high ({a_2})")
    #plt.plot(np.array(betas)/params["gamma"], outbreak_probabilities_2, marker='o', linestyle='-', label=f"Assortativity low")
    
    plt.xlabel(r"$\beta$ (Infection Rate)", fontsize=20)
    plt.ylabel("Strain 1 Prevalence at Equilibrium (People Infected)", fontsize=20)
    plt.grid(True)
    # plt.ylim([-0.01, 1.02])
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title("Viral Prevalence vs. Beta for Different Assortativities", fontsize=20)
    #plt.ylim([-0.01, 1.02])
    plt.legend()
    plt.show()

    ## Plot of peak infected versus beta for both assortativity levels
    plt.figure(figsize=(8, 5))
    plt.plot(np.array(betas), peak_infected_1, marker='o', linestyle='-', linewidth=4, label=f"Assortativity low ({a_1})")
    plt.plot(np.array(betas), peak_infected_2, marker='o', linestyle='-', linewidth=4, label=f"Assortativity high ({a_2})")
    plt.xlabel(r"$\beta$ (Infection Rate)", fontsize=20)
    plt.ylabel("Peak Infected (Strain 1)", fontsize=20)
    plt.grid(True)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title("Peak Infected vs. Beta for Different Assortativities", fontsize=20)
    plt.legend()
    plt.show()


    
if __name__ == "__main__":
    main()

