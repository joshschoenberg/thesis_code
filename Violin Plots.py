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

import numpy as np

def run_simulation(sim_index, params):
    num_nodes = params.get('num_nodes', 600)
    connection_bias = params.get('connection_bias', 1)
    initial_infected_number = params.get('initial_infected_number', 1)
    time_steps = params.get('time_steps', 100)
    gamma = params.get('gamma', 0.2)
    avg_connectivity = params.get('avg_connectivity', 10)
    beta = params.get('beta', gamma*2.5/avg_connectivity)
    vaccine_effectiveness = params.get('vaccine_effectiveness', 1)
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
    ## Vaccinate people before simulation starts
    susceptible_nodes = np.where(S == 1)[0]
    if len(susceptible_nodes) > 0:
        for i in susceptible_nodes:
            if G.nodes[i]['type'] == 0 and np.random.rand() < p_vax:
                S[i] = 0
                V[i] = 1

    for t in range(time_steps):
        new_infected_1 = set()
        new_infected_2 = set()
        new_recovered_1 = set()
        new_recovered_2 = set()
       
        infected_nodes_var_1 = np.where(I1 == 1)[0]
        infected_nodes_var_2 = np.where(I2 == 1)[0]
        all_infected_nodes = np.concatenate((infected_nodes_var_1, infected_nodes_var_2))
        infected_nodes_random_order = np.random.permutation(all_infected_nodes)
        for i in infected_nodes_random_order:
            if I1[i] == 1:
                if np.random.rand() < gamma:
                    new_recovered_1.add(i)
                for neighbor in G.neighbors(i):
                    if S[neighbor] == 1 and np.random.rand() < beta:
                        if np.random.rand() < p_mutation:
                            if (neighbor not in new_infected_1) and (neighbor not in new_recovered_1) and (neighbor not in new_recovered_2):
                              new_infected_2.add(neighbor)
                        else:
                            if neighbor not in new_infected_2 and (neighbor not in new_recovered_1) and (neighbor not in new_recovered_2):
                              new_infected_1.add(neighbor)
            elif I2[i] == 1:
                if np.random.rand() < gamma:
                    new_recovered_2.add(i)
                for neighbor in G.neighbors(i):
                    if S[neighbor] == 1 and np.random.rand() < beta:
                        if (neighbor not in new_infected_1) and (neighbor not in new_recovered_1) and (neighbor not in new_recovered_2):
                          new_infected_2.add(neighbor)

                    elif (V[neighbor] == 1 or R1[neighbor] == 1) and (np.random.rand() < beta * p_breakthrough):
                      if (neighbor not in new_infected_1) and (neighbor not in new_recovered_1) and (neighbor not in new_recovered_2):
                        new_infected_2.add(neighbor)      

        for node in new_infected_1:
            I1[node] = 1
            V[node] = 0
            S[node] = 0
        for node in new_infected_2:
            I2[node] = 1
            V[node] = 0
            S[node] = 0
        for node in new_recovered_1:
            R1[node] = 1
            I1[node] = 0
        for node in new_recovered_2:
            R2[node] = 1
            I2[node] = 0

    total_infected = np.sum(R1) + np.sum(R2) + np.sum(I1) + np.sum(I2)
    total_infected_var_1 = np.sum(R1) + np.sum(I1)
    total_infected_var_2 = np.sum(R2) + np.sum(I2)
    # time_of_peak_infection = np.argmax(i_over_time)
    percent_considered_outbreak = 0.1
    outbreak_occurred_var_1 = total_infected_var_1 > (percent_considered_outbreak * num_nodes)
    outbreak_occurred_var_2 = total_infected_var_2 > (percent_considered_outbreak * num_nodes)
    outbreak_occurred = total_infected > (percent_considered_outbreak * num_nodes)
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
    p_mutation = param_set.get('p_mutation', 0.01)
    vaccine_effectiveness = param_set.get('vaccine_effectiveness', 1)

    params = {'connection_bias': connection_bias, 'avg_connectivity': avg_connectivity,
              'gamma': gamma, 'beta': beta, 'p_vax': p_vax, 'num_nodes':num_nodes,
              'initial_infected_number': initial_infected_number,
              'time_steps':time_steps, 'p_mutation': p_mutation,
              'vaccine_effectiveness':vaccine_effectiveness}
    ## Run simulations in parallel
    with Pool() as pool:
        # results = pool.starmap(run_simulation, [(i, connection_bias, beta, gamma, p_vax, num_nodes, initial_infected_number, time_steps) for i in range(num_simulations)])
        results = pool.starmap(run_simulation, [(i, params) for i in range(num_simulations)])
    return results

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Violin plot for total infections comparison
def plot_violin(total_infected, assortativity_values, title):
    flattened_data = [val for sublist in total_infected for val in sublist]
    group_labels = [f'Assortativity {assortativity_values[i]:.2f}' for i, sublist in enumerate(total_infected) for _ in sublist]

    # Define custom colors that blend well with Deep Teal and Burnt Orange
    custom_palette = {
        f'Assortativity {assortativity_values[0]:.2f}': '#00778B',  # Deep Teal
        f'Assortativity {assortativity_values[1]:.2f}': '#009E73',  # Teal-Green
        f'Assortativity {assortativity_values[2]:.2f}': '#F0E442',  # Yellow
        f'Assortativity {assortativity_values[3]:.2f}': '#D55E00',  # Dark Orange
        # f'Assortativity {assortativity_values[4]:.2f}': '#E66101'   # Burnt Orange
    }
    
    plt.figure(figsize=(10, 6))
    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sns.violinplot(x=group_labels, y=flattened_data, density_norm="area", cut=0, palette=custom_palette)
    # Overlay strip plot to show individual points
    sns.stripplot(x=group_labels, y=flattened_data, jitter=True, dodge=True, size=5, color="black", alpha=0.6)

    plt.ylabel('Size of Outbreak', fontsize=20)
    plt.title(title, fontsize=20, pad=15)
    plt.xticks(fontsize=20, rotation=45)
    plt.yticks(fontsize=18)
    plt.grid(axis='y', alpha=0.3)
    plt.show()

# Example: Run simulations for different assortativity values
def main():
    num_simulations = 600
    num_nodes = 2000
    gamma = 0.2
    avg_connections = 10
    beta = gamma * 3.5 / avg_connections
    assortativity_levels = [1, 1.42, 3.5, 60]  # Different connection bias values
    
    total_infected = []
    var_1_infected = []
    var_2_infected = []
    avg_assortativity_values = []

    for connection_bias in assortativity_levels:
        param_set = {'connection_bias': connection_bias, 'num_nodes': num_nodes, 'beta': beta, 'num_simulations': num_simulations}
        results = run_parallel_simulations(param_set)

        total_infected.append([result[0] for result in results])
        var_1_infected.append([result[1] for result in results])
        var_2_infected.append([result[2] for result in results])
        avg_assortativity_values.append(np.abs(np.mean([result[6] for result in results])))
    
    print(f"Avg. assortativity values: {avg_assortativity_values}")
    
    plot_violin(total_infected, avg_assortativity_values, 'Size of Outbreak by Assortativity Level (Both Variants)')
    plot_violin(var_1_infected, avg_assortativity_values, 'Size of Outbreak by Assortativity Level (Variant 1)')
    plot_violin(var_2_infected, avg_assortativity_values, 'Size of Outbreak by Assortativity Level (Variant 2)')

if __name__ == "__main__":
    main()
