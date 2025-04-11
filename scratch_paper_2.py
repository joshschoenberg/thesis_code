### EXPERIMENTING WITH FASTER RUNTIME
### Violin Plot for Size of Outbreak v.s. Assortativity Level (2000 nodes, gamma=0.2, R0=3.5, Assortativity = 0 and Assortativity = 1)

import numpy as np
import networkx as nx
import scipy.stats as stats
import matplotlib.pyplot as plt
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

def run_simulation(sim_index, params):
    num_nodes = params.get('num_nodes', 600)
    connection_bias = params.get('connection_bias', 1)
    initial_infected_number = params.get('initial_infected_number', 1)
    time_steps = params.get('time_steps', 100)
    gamma = params.get('gamma', 0.2)
    avg_connectivity = params.get('avg_connectivity', 10)
    beta = params.get('beta', gamma*2.5/avg_connectivity)
    p_vax = params.get('p_vax', 1)
    p_mutation = params.get('p_mutation', 0.01)  # 1% chance of mutation per transmission
    p_breakthrough = params.get('p_breakthrough', 0.5)  # 50% chance of breakthrough infection

    G, assortativity = create_graph(connection_bias, num_nodes)

    S = np.ones(num_nodes)
    I1 = np.zeros(num_nodes)  # First strain
    I2 = np.zeros(num_nodes)  # Mutated strain
    R1 = np.zeros(num_nodes)
    R2 = np.zeros(num_nodes)
    V = np.zeros(num_nodes)

    initial_infected_nodes = np.random.choice(range(num_nodes), initial_infected_number, replace=False)
    
    I1[initial_infected_nodes] = 1
    S[initial_infected_nodes] = 0
    for t in range(time_steps):
        # Convert 'type' attributes of nodes to a NumPy array for efficient indexing
        types = np.array([G.nodes[i]['type'] for i in range(num_nodes)])

        # Find indices where conditions are met
        eligible_nodes = (S == 1) & (types == 0) & (np.random.rand(num_nodes) < p_vax)

        # Apply updates in one step
        S[eligible_nodes] = 0
        V[eligible_nodes] = 1
        
        # Process Recovery
        new_recovered_1 = np.where((I1 == 1) & (np.random.rand(num_nodes) < gamma))[0]
        new_recovered_2 = np.where((I2 == 1) & (np.random.rand(num_nodes) < gamma))[0]
    
        # Identify currently infected individuals
        infected_1_nodes = np.where(I1 == 1)[0]
        infected_2_nodes = np.where(I2 == 1)[0]

        # # Find susceptible neighbors of infected individuals
        # all_neighbors = {i: list(G.neighbors(i)) for i in np.concatenate((infected_1_nodes, infected_2_nodes))}
        
        # neighbor_nodes = np.array([n for neighbors in all_neighbors.values() for n in neighbors], dtype=int)
        # # Find susceptible neighbors
        # susceptible_neighbors = neighbor_nodes[S[neighbor_nodes] == 1]
        # Collect all susceptible neighbors of infected individuals
        susceptible_neighbors = set()
        vaccinated_neighbors = set()
        variant_1_recovered_neighbors = set()

        for i in np.concatenate((infected_1_nodes, infected_2_nodes)):
            for neighbor in G.neighbors(i):
                if S[neighbor] == 1:  # Check if the neighbor is susceptible
                    susceptible_neighbors.add(neighbor)
                if V[neighbor] == 1:
                    vaccinated_neighbors.add(neighbor)  
                if R1[neighbor] == 1:
                    variant_1_recovered_neighbors.add(neighbor)

        susceptible_neighbors = np.array(list(susceptible_neighbors))  # Convert to numpy array
        vaccinated_neighbors = np.array(list(vaccinated_neighbors))
        variant_1_recovered_neighbors = np.array(list(variant_1_recovered_neighbors))
        vax_recovered_neighbors = np.concatenate((vaccinated_neighbors, variant_1_recovered_neighbors))

        # Determine exposure
        # Find neighbors exposed to each infection type
        exposed_to_I1 = np.array(
    [any(neigh in infected_1_nodes for neigh in G.neighbors(node)) for node in susceptible_neighbors], dtype=bool
)
        exposed_to_I2 = np.array(
    [any(neigh in infected_2_nodes for neigh in G.neighbors(node)) for node in susceptible_neighbors], dtype=bool
)
        exposed_vax_recovered_to_I2 = np.array(
    [any(neigh in infected_2_nodes for neigh in G.neighbors(node)) for node in vax_recovered_neighbors], dtype=bool
)
      # Infection probability check
        infection_probs = np.random.rand(len(susceptible_neighbors)) < beta
    
        # Mutation probability check (Virus 1 mutating into Virus 2)
        mutations = np.random.rand(len(susceptible_neighbors)) < p_mutation
    
        # Handle dual exposure
        dual_exposure = exposed_to_I1 & exposed_to_I2
    
        # Assign new infections
        new_infected_I1 = susceptible_neighbors[infection_probs & exposed_to_I1 & ~mutations & ~dual_exposure]
        new_infected_I2 = susceptible_neighbors[infection_probs & exposed_to_I2]
        new_infected_mutated = susceptible_neighbors[infection_probs & exposed_to_I1 & mutations]
    
        # Randomly assign dual exposures
        dual_infected_I1 = susceptible_neighbors[dual_exposure & infection_probs & (np.random.rand(len(susceptible_neighbors)) < 0.5)]
        dual_infected_I2 = susceptible_neighbors[dual_exposure & infection_probs & (np.random.rand(len(susceptible_neighbors)) >= 0.5)]
    
        # Infect vaccinated/recovered individuals 
        # Breakthrough infections for vaccinated and recovered individuals
        breakthrough_probs = np.random.rand(len(vax_recovered_neighbors)) < (beta * p_breakthrough)
        new_breakthrough_I2 = np.where(breakthrough_probs & exposed_vax_recovered_to_I2)[0]

        # Update infection states
    
        I1[list(new_infected_I1) + list(dual_infected_I1)] = 1  # Standard I1 infections
        I2[list(new_infected_I2) + list(new_infected_mutated) + list(dual_infected_I2) + list(new_breakthrough_I2)] = 1  # Include breakthrough infections
        S[list(new_infected_I1) + list(new_infected_I2) + list(new_infected_mutated) + list(dual_infected_I1) + list(dual_infected_I2)] = 0

        # Remove vaccinated or recovered status for breakthrough cases
        V[new_breakthrough_I2] = 0
        R1[new_breakthrough_I2] = 0
        R2[new_breakthrough_I2] = 0
                # Move recovered individuals to R1 or R2
        R1[new_recovered_1] = 1
        I1[new_recovered_1] = 0
        R2[new_recovered_2] = 1
        I2[new_recovered_2] = 0

    total_infected = np.sum(R1) + np.sum(R2) + np.sum(I1) + np.sum(I2)
    # time_of_peak_infection = np.argmax(i_over_time)
    total_infected_var_1 = np.sum(R1) + np.sum(I1)
    total_infected_var_2 = np.sum(R2) + np.sum(I2)
    percent_considered_outbreak = 0.01
    outbreak_occurred = total_infected > (percent_considered_outbreak * num_nodes)
    outbreak_occurred_var_1 = total_infected_var_1 > (percent_considered_outbreak * num_nodes)
    outbreak_occurred_var_2 = total_infected_var_2 > (percent_considered_outbreak * num_nodes)
    return total_infected, total_infected_var_1, total_infected_var_2, outbreak_occurred, outbreak_occurred_var_1, outbreak_occurred_var_2, assortativity

def run_parallel_simulations(param_set):
    """Run multiple simulations for given parameter set in parallel."""
    ## Retrieve parameters, or use default values
    connection_bias = param_set.get('connection_bias', 1)
    avg_connectivity = param_set.get('avg_connectivity', 10)
    gamma = param_set.get('gamma', 0.2)
    beta = param_set.get('beta', gamma * 2.5 / avg_connectivity)
    p_vax = param_set.get('p_vax', 1)
    num_nodes = param_set.get('num_nodes', 600)
    initial_infected_number = param_set.get('initial_infected_number', 1)
    time_steps = param_set.get('time_steps', 100)
    num_simulations = param_set.get('num_simulations', 100)
    params = {'connection_bias': connection_bias, 'avg_connectivity': avg_connectivity,
              'gamma': gamma, 'beta': beta, 'p_vax': p_vax, 'num_nodes':num_nodes,
              'initial_infected_number': initial_infected_number,
              'time_steps':time_steps}
    ## Run simulations in parallel
    with Pool() as pool:
        # results = pool.starmap(run_simulation, [(i, connection_bias, beta, gamma, p_vax, num_nodes, initial_infected_number, time_steps) for i in range(num_simulations)])
        results = pool.starmap(run_simulation, [(i, params) for i in range(num_simulations)])
    return results

# Violin plot for total infections comparison
def plot_violin(total_infected_1, total_infected_2, labels, title):
    data = [total_infected_1, total_infected_2]
    flattened_data = [val for sublist in data for val in sublist]
    group_labels = [labels[i] for i, sublist in enumerate(data) for _ in sublist]

    # Custom colors for the violin plot
    custom_palette = {
        'Assortativity 0': '#00778B',  # Deep Teal for assortativity 0
        'Assortativity 1': '#E66101'   # Burnt Orange for assortativity 1
    }

    plt.figure(figsize=(10, 6))
    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sns.violinplot(x=group_labels, y=flattened_data, density_norm="area", cut=0, palette=custom_palette)
    # Overlay strip plot to show individual points
    sns.stripplot(x=group_labels, y=flattened_data, jitter=True, dodge=True, size=5, color="black", alpha=0.6)

    # plt.xlabel('Assortativity Level', fontsize=20)
    plt.ylabel('Size of Outbreak', fontsize=20)
    plt.title(title, fontsize=20, pad=15)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=18)
    plt.grid(axis='y', alpha=0.3)
    plt.show()


# Example: Run simulations for 15 beta values
def main():
    num_simulations = 600
    total_infected_assort_0 = []
    var_1_infected_assort_0 = []
    var_2_infected_assort_0 = []
    total_infected_assort_1 = []
    var_1_infected_assort_1 = []
    var_2_infected_assort_1 = []
    peak_infection_times_0 = []
    peak_infection_times_1 = []
    outbreak_likelihoods_0 = []
    outbreak_likelihoods_1 = []
    outbreak_likelihoods_var_1_0 = []
    outbreak_likelihoods_var_1_1 = []
    outbreak_likelihoods_var_2_0 = []
    outbreak_likelihoods_var_2_1 = []
    outbreak_var_2_likelihood_given_var_1_outbreak_0 = []
    outbreak_var_2_likelihood_given_var_1_outbreak_1 = []

    avg_connectivity = 10
    num_nodes = 2000
    gamma_value = 0.2
    r_values = np.linspace(0.5, 8, 20)
    beta_values = r_values * gamma_value / avg_connectivity

    # Run simulations for Assortativity = 0
    for beta in beta_values:
        param_set = {'beta': beta, 'connection_bias': 1, 'num_simulations': num_simulations, 'num_nodes': num_nodes}
        results = run_parallel_simulations(param_set)

        total_infected = [result[0] for result in results]
        var_1_infected = [result[1] for result in results]
        var_2_infected = [result[2] for result in results]
        outbreak_occurred = [result[3] for result in results]
        outbreak_var_1_occurred = [result[4] for result in results]
        outbreak_var_2_occurred = [result[5] for result in results]
        # Create a new list that only includes values of strain2_outbreaks where strain1_outbreaks is 1
        outbreak_var_2_occurred_given_outbreak_var_1 = [s2 for s1, s2 in zip(outbreak_var_1_occurred, outbreak_var_2_occurred) if s1 == 1]
        # peak_times = [result[6] for result in results]

        total_infected_assort_0.append(np.mean(total_infected))
        var_1_infected_assort_0.append(np.mean(var_1_infected))
        var_2_infected_assort_0.append(np.mean(var_2_infected))
        # peak_infection_times_0.append(np.mean(peak_times))
        outbreak_likelihoods_0.append(np.mean(outbreak_occurred))
        outbreak_likelihoods_var_1_0.append(np.mean(outbreak_var_1_occurred))
        outbreak_likelihoods_var_2_0.append(np.mean(outbreak_var_2_occurred))
        outbreak_var_2_likelihood_given_var_1_outbreak_0.append(np.mean(outbreak_var_2_occurred_given_outbreak_var_1))
        print(f"finished {beta} value")

    # Run simulations for Assortativity = 1
    for beta in beta_values:
        param_set = {'beta': beta, 'connection_bias': 90, 'num_simulations': num_simulations,
                     'num_nodes': num_nodes}
        results = run_parallel_simulations(param_set)

        total_infected = [result[0] for result in results]
        var_1_infected = [result[1] for result in results]
        var_2_infected = [result[2] for result in results]
        outbreak_occurred = [result[3] for result in results]
        outbreak_var_1_occurred = [result[4] for result in results]
        outbreak_var_2_occurred = [result[5] for result in results]
        outbreak_var_2_occurred_given_outbreak_var_1 = [s2 for s1, s2 in zip(outbreak_var_1_occurred, outbreak_var_2_occurred) if s1 == 1]

        # peak_times = [result[6] for result in results]

        total_infected_assort_1.append(np.mean(total_infected))
        var_1_infected_assort_1.append(np.mean(var_1_infected))
        var_2_infected_assort_1.append(np.mean(var_2_infected))
        # peak_infection_times_1.append(np.mean(peak_times))
        outbreak_likelihoods_1.append(np.mean(outbreak_occurred))
        outbreak_likelihoods_var_1_1.append(np.mean(outbreak_var_1_occurred))
        outbreak_likelihoods_var_2_1.append(np.mean(outbreak_var_2_occurred))
        outbreak_var_2_likelihood_given_var_1_outbreak_1.append(np.mean(outbreak_var_2_occurred_given_outbreak_var_1))

        print(f"finished {beta} value")

    colors = {
        'assort_0': '#E66101',
        'assort_1': '#1B3A6D'
    }
    # **Plot 1a: Beta vs. Total Infections**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, total_infected_assort_0, label='Assortativity 0', color=colors['assort_0'], marker='o', linewidth=4)
    plt.plot(beta_values, total_infected_assort_1, label='Assortativity 1',color=colors['assort_1'],  marker='o', linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Average Total Infected', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Total Infections vs Beta (All Infections)', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    ## Plot 1b: Beta v.s. Strain 1 Infections
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, var_1_infected_assort_0, label='Assortativity 0', color=colors['assort_0'], marker='o', linewidth=4)
    plt.plot(beta_values, var_1_infected_assort_1, label='Assortativity 1',color=colors['assort_1'],  marker='o', linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Average Total Infected', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Total Infections vs Beta (Strain 1)', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    ## Plot 1c: Beta v.s. Strain 2 Infections
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, var_2_infected_assort_0, label='Assortativity 0', color=colors['assort_0'], marker='o', linewidth=4)
    plt.plot(beta_values, var_2_infected_assort_1, label='Assortativity 1',color=colors['assort_1'],  marker='o', linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Average Total Infected', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Total Infections vs Beta (Strain 2)', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()


    # # **Plot 2: Beta vs. Time of Peak Infection**
    # plt.figure(figsize=(10, 6))
    # plt.plot(beta_values, peak_infection_times_0, label='Assortativity 0', color=colors['assort_0'], marker='s', linewidth=4)
    # plt.plot(beta_values, peak_infection_times_1, label='Assortativity 1', color=colors['assort_1'], marker='s', linewidth=4)

    # plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    # plt.ylabel('Time of Peak Infection', fontsize=20)
    # plt.xticks(fontsize=18)
    # plt.yticks(fontsize=18)
    # plt.legend(fontsize=18)
    # plt.title('Beta vs. Time of Peak Infection for Different Assortativities', fontsize=20)

    # plt.grid(False)
    # ax = plt.gca()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # plt.show()

    # **Plot 3a: Beta vs. Likelihood of Outbreak**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, outbreak_likelihoods_0, label='Assortativity 0', color=colors['assort_0'], marker='D', linewidth=4)
    plt.plot(beta_values, outbreak_likelihoods_1, label='Assortativity 1', marker='D', color=colors['assort_1'], linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Likelihood of Outbreak', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Beta vs. Likelihood of Outbreak', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    # **Plot 3b: Beta vs. Likelihood of Outbreak (Strain 1)**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, outbreak_likelihoods_var_1_0, label='Assortativity 0', color=colors['assort_0'], marker='D', linewidth=4)
    plt.plot(beta_values, outbreak_likelihoods_var_1_1, label='Assortativity 1', marker='D', color=colors['assort_1'], linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Likelihood of Outbreak', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Beta vs. Likelihood of Outbreak (of Strain 1)', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    # **Plot 3c: Beta vs. Likelihood of Outbreak (Strain 2)**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, outbreak_likelihoods_var_2_0, label='Assortativity 0', color=colors['assort_0'], marker='D', linewidth=4)
    plt.plot(beta_values, outbreak_likelihoods_var_2_1, label='Assortativity 1', marker='D', color=colors['assort_1'], linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Likelihood of Outbreak', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Beta vs. Likelihood of Outbreak (of Strain 2)', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    # **Plot 3d: Beta vs. Likelihood of Outbreak Strain 2, Given Outbreak of Strain 1**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, outbreak_var_2_likelihood_given_var_1_outbreak_0, label='Assortativity 0', color=colors['assort_0'], marker='D', linewidth=4)
    plt.plot(beta_values, outbreak_var_2_likelihood_given_var_1_outbreak_1, label='Assortativity 1', marker='D', color=colors['assort_1'], linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Likelihood of Outbreak', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Beta vs. Likelihood of Outbreak of Strain 2, Given Outbreak of Strain 1', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()


if __name__ == "__main__":
    main()