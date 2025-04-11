### Analyzing Beta Value v.s. Number of Total Infections (2000 nodes, gamma=0.2, R0=3.5, Assortativity = 0 and Assortativity = 1)

import numpy as np
import networkx as nx
import scipy.stats as stats
import matplotlib.pyplot as plt
from multiprocessing import Pool

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
    vaccine_efficacy = params.get('vaccine_efficacy', 1)

    """Run a single SIR simulation."""
    G, assortativity = create_graph(connection_bias, num_nodes)

    S = np.ones(num_nodes)
    I = np.zeros(num_nodes)
    R = np.zeros(num_nodes)
    V = np.zeros(num_nodes)

    initial_infected_nodes = np.random.choice(range(num_nodes), initial_infected_number, replace=False)
    I[initial_infected_nodes] = 1
    S[initial_infected_nodes] = 0

    s_over_time = [np.sum(S)]
    i_over_time = [np.sum(I)]
    r_over_time = [np.sum(R)]
    v_over_time = [np.sum(V)]
    for i in range(num_nodes):
            if S[i] == 1 and G.nodes[i]['type'] == 0 and np.random.rand() < p_vax:
                S[i] = 0
                V[i] = 1

    pos = None
    for t in range(time_steps):
        new_infected = []
        new_recovered = []
        infected_nodes = np.where(I == 1)[0]
        for i in infected_nodes:
            if np.random.rand() < gamma:
                new_recovered.append(i)
            for neighbor in G.neighbors(i):
                if S[neighbor] == 1 and np.random.rand() < beta:
                    new_infected.append(neighbor)
                # if V[neighbor] == 1 and np.random.rand() < beta and np.random.rand() < vaccine_efficacy:
                #     new_infected.append(neighbor)
        for node in new_infected:
            I[node] = 1
            S[node] = 0
            V[node] = 0
        for node in new_recovered:
            R[node] = 1
            I[node] = 0
            V[node] = 0
        s_over_time.append(np.sum(S))
        i_over_time.append(np.sum(I))
        r_over_time.append(np.sum(R))
        v_over_time.append(np.sum(V))

    total_infected = np.sum(R) + np.sum(I)
    time_of_peak_infection = np.argmax(i_over_time)
    percent_considered_outbreak = 0.01
    outbreak_occurred = total_infected > (percent_considered_outbreak * num_nodes)
    return total_infected, assortativity, time_of_peak_infection, outbreak_occurred

def run_parallel_simulations(param_set):
    """Run multiple simulations for given parameter set in parallel."""
    ## Retrieve parameters, or use default values
    connection_bias = param_set.get('connection_bias', 1)
    avg_connectivity = param_set.get('avg_connectivity', 10)
    gamma = param_set.get('gamma', 0.2)
    beta = param_set.get('beta', gamma * 2.5 / avg_connectivity)
    p_vax = param_set.get('p_vax', 1)
    num_nodes = param_set.get('num_nodes', 600)
    initial_infected_number = param_set.get('inital_infected_number', 1)
    time_steps = param_set.get('time_steps', 100)
    num_simulations = param_set.get('num_simulations', 100)
    vaccine_efficacy = param_set.get('vaccine_efficacy', 1)
    params = {'connection_bias': connection_bias, 'avg_connectivity': avg_connectivity,
              'gamma': gamma, 'beta': beta, 'p_vax': p_vax, 'num_nodes':num_nodes,
              'initial_infected_number': initial_infected_number,
              'time_steps':time_steps, 'vaccine_efficacy':vaccine_efficacy}
    ## Run simulations in parallel
    with Pool() as pool:
        # results = pool.starmap(run_simulation, [(i, connection_bias, beta, gamma, p_vax, num_nodes, initial_infected_number, time_steps) for i in range(num_simulations)])
        results = pool.starmap(run_simulation, [(i, params) for i in range(num_simulations)])
    return results

# Example: Run simulations for 15 beta values
def main():
    num_simulations = 800
    num_nodes = 2000
    total_infected_assort_0 = []
    total_infected_assort_1 = []
    total_infected_assort_05 = []
    total_infected_assort_03 = []
    assortativities_0 = []
    assortativities_1 = []
    assortativities_05 = []
    assortativities_03 = []
    peak_infection_times_0 = []
    peak_infection_times_1 = []
    peak_infection_times_05 = []
    peak_infection_times_03 = []
    outbreak_likelihoods_0 = []
    outbreak_likelihoods_1 = []
    outbreak_likelihoods_05 = []
    outbreak_likelihoods_03 = []
    avg_connectivity = 10
    gamma_value = 0.2
    r_values = np.linspace(0.5, 8, 20)
    beta_values = r_values * gamma_value / avg_connectivity
    vaccine_efficacy = 1
    # Run simulations for Assortativity = 0
    for beta in beta_values:
        param_set = {'beta': beta, 'connection_bias': 1, 'num_simulations': num_simulations, 
                     'num_nodes':num_nodes, 'vaccine_efficacy':vaccine_efficacy}
        results = run_parallel_simulations(param_set)

        total_infected = [result[0] for result in results]
        assortativity = [result[1] for result in results]
        peak_times = [result[2] for result in results]
        outbreak_occurred = [result[3] for result in results]

        assortativities_0.append(np.mean(assortativity))
        total_infected_assort_0.append(np.mean(total_infected))
        peak_infection_times_0.append(np.mean(peak_times))
        outbreak_likelihoods_0.append(np.mean(outbreak_occurred))
        print(f"finished beta {beta} for assortativity 0")

    # Run simulations for Assortativity = 1
    for beta in beta_values:
        param_set = {'beta': beta, 'connection_bias': 90, 'num_simulations': num_simulations, 
                     'num_nodes':num_nodes, 'vaccine_efficacy:':vaccine_efficacy}
        results = run_parallel_simulations(param_set)

        total_infected = [result[0] for result in results]
        assortativity = [result[1] for result in results]
        peak_times = [result[2] for result in results]
        outbreak_occurred = [result[3] for result in results]

        total_infected_assort_1.append(np.mean(total_infected))
        peak_infection_times_1.append(np.mean(peak_times))
        outbreak_likelihoods_1.append(np.mean(outbreak_occurred))
        assortativities_1.append(np.mean(assortativity))

        print(f"finished beta {beta} for assortativity 1")
    
    # Run simulations for Assortativity = 0.3
    for beta in beta_values:
        param_set = {'beta': beta, 'connection_bias': 1.5, 'num_simulations': num_simulations, 
                     'num_nodes':num_nodes, 'vaccine_efficacy':vaccine_efficacy}
        results = run_parallel_simulations(param_set)

        total_infected = [result[0] for result in results]
        assortativity = [result[1] for result in results]
        peak_times = [result[2] for result in results]
        outbreak_occurred = [result[3] for result in results]

        total_infected_assort_03.append(np.mean(total_infected))
        peak_infection_times_03.append(np.mean(peak_times))
        outbreak_likelihoods_03.append(np.mean(outbreak_occurred))
        assortativities_03.append(np.mean(assortativity))
        print(f"finished beta {beta} for assortativity 0.3")
    
    # Run simulations for Assortativity = 0.5
    for beta in beta_values:
        param_set = {'beta': beta, 'connection_bias': 2, 'num_simulations': num_simulations, 
                     'num_nodes':num_nodes, 'vaccine_efficacy':vaccine_efficacy}
        results = run_parallel_simulations(param_set)

        total_infected = [result[0] for result in results]
        assortativity = [result[1] for result in results]
        peak_times = [result[2] for result in results]
        outbreak_occurred = [result[3] for result in results]

        total_infected_assort_05.append(np.mean(total_infected))
        peak_infection_times_05.append(np.mean(peak_times))
        outbreak_likelihoods_05.append(np.mean(outbreak_occurred))
        assortativities_05.append(np.mean(assortativity))
        print(f"finished beta {beta} for assortativity 0.5")

    colors = {
        'assort_0': '#E66101',
        'assort_1': '#1B3A6D',
        'assort_05': '#5D9E3E',
        'assort_03': '#F1C232'


    }

    assort_0 = np.abs(np.mean(assortativities_0))
    assort_1 = np.mean(assortativities_1)
    assort_05 = np.mean(assortativities_05)
    assort_03 = np.mean(assortativities_03)
    # **Plot 1: Beta vs. Total Infections**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, total_infected_assort_0, label=f'Assortativity {assort_0:.2f}', color=colors['assort_0'], marker='o', linewidth=4)
    plt.plot(beta_values, total_infected_assort_03, label=f'Assortativity {assort_03:.2f}',color=colors['assort_03'],  marker='o', linewidth=4)
    plt.plot(beta_values, total_infected_assort_05, label=f'Assortativity {assort_05:.2f}', color=colors['assort_05'], marker='o', linewidth=4)
    plt.plot(beta_values, total_infected_assort_1, label=f'Assortativity {assort_1:.2f}',color=colors['assort_1'],  marker='o', linewidth=4)
    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Average Total Infected', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Total Infections vs Beta for Different Assortativities', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    # **Plot 2: Beta vs. Time of Peak Infection**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, peak_infection_times_0, label=f'Assortativity {assort_0:.2f}', color=colors['assort_0'], marker='s', linewidth=4)
    plt.plot(beta_values, peak_infection_times_03, label=f'Assortativity {assort_03:.2f}', color=colors['assort_03'], marker='s', linewidth=4)
    plt.plot(beta_values, peak_infection_times_05, label=f'Assortativity {assort_05:.2f}', color=colors['assort_05'], marker='s', linewidth=4)
    plt.plot(beta_values, peak_infection_times_1, label=f'Assortativity {assort_1:.2f}', color=colors['assort_1'], marker='s', linewidth=4)
    
    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Time of Peak Infection', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Beta vs. Time of Peak Infection for Different Assortativities', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

    # **Plot 3: Beta vs. Likelihood of Outbreak**
    plt.figure(figsize=(10, 6))
    plt.plot(beta_values, outbreak_likelihoods_0, label=f'Assortativity {assort_0:.2f}', color=colors['assort_0'], marker='D', linewidth=4)
    plt.plot(beta_values, outbreak_likelihoods_03, label=f'Assortativity {assort_03:.2f}', marker='D', color=colors['assort_03'], linewidth=4)    
    plt.plot(beta_values, outbreak_likelihoods_05, label=f'Assortativity {assort_05:.2f}', color=colors['assort_05'], marker='D', linewidth=4)
    plt.plot(beta_values, outbreak_likelihoods_1, label=f'Assortativity {assort_1:.2f}', marker='D', color=colors['assort_1'], linewidth=4)

    plt.xlabel('Beta Value (Infection Rate)', fontsize=20)
    plt.ylabel('Likelihood of Outbreak', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.title('Beta vs. Likelihood of Outbreak for Different Assortativities', fontsize=20)

    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()


if __name__ == "__main__":
    main()
