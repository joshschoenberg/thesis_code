import numpy as np
import matplotlib.pyplot as plt
# Low Assortativity (~0.0001 - 0.0017)
# beta_low = [0.10, 0.27, 0.44, 0.61, 0.79, 0.96, 1.13, 1.30, 1.47, 1.64, 1.81, 1.99, 2.16, 2.33, 2.50]
# prob_low = [0.000, 0.000, 0.000, 0.000, 0.000, 0.003, 0.023, 0.120, 0.253, 0.323, 0.450, 0.510, 0.530, 0.507, 0.503]

# # Medium-Low Assortativity (~0.33)
# beta_medlow = [0.10, 0.27, 0.44, 0.61, 0.79, 0.96, 1.13, 1.30, 1.47, 1.64, 1.81, 1.99, 2.16, 2.33, 2.50]
# prob_medlow = [0.000, 0.000, 0.000, 0.150, 0.460, 0.490, 0.680, 0.713, 0.773, 0.807, 0.820, 0.843, 0.873, 0.860, 0.837]

# # Medium-High Assortativity (~0.67)
# beta_medhigh = [0.10, 0.27, 0.44, 0.61, 0.79, 0.96, 1.13, 1.30, 1.47, 1.64, 1.81, 1.99, 2.16, 2.33, 2.50]
# prob_medhigh = [0.000, 0.000, 0.293, 0.563, 0.683, 0.797, 0.840, 0.870, 0.887, 0.923, 0.923, 0.963, 0.967, 0.963, 0.953]

# # High Assortativity (~0.99)
# beta_high = [0.10, 0.27, 0.44, 0.61, 0.79, 0.96, 1.13, 1.30, 1.47, 1.64, 1.81, 1.99, 2.16, 2.33, 2.50]
# prob_high = [0.000, 0.030, 0.533, 0.710, 0.840, 0.877, 0.907, 0.920, 0.957, 0.947, 0.983, 0.993, 1.000, 1.000, 1.000]

# assort_low = 0.000
# assort_medlow = 0.33
# assort_medhigh = 0.67
# assort_high = 0.99

# plt.figure(figsize=(10, 6))

# plt.plot(beta_low, prob_low, marker='o', linestyle='-', linewidth=4, label=f"{assort_low}")
# plt.plot(beta_medlow, prob_medlow, marker='o', linestyle='-', linewidth=4, label=f"{assort_medlow}")
# plt.plot(beta_medhigh, prob_medhigh, marker='o', linestyle='-', linewidth=4, label=f"{assort_medhigh}")
# plt.plot(beta_high, prob_high, marker='o', linestyle='-', linewidth=4, label=f"{assort_high}")

#     # -------------------------
#     # Plot Outbreak Probability vs. Beta
#     # -------------------------
# plt.xlabel(r"$\beta$ (Infection Rate)", fontsize=20)
# plt.ylabel("Probability of Outbreak", fontsize=20)
# plt.title("Outbreak Probability vs. Beta for Different Assortativities", fontsize=20)
# plt.legend(title="Assortativity", fontsize=18)
# plt.grid(True)
# # plt.ylim([-0.01, 1.02])
# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.show()

# betas = [0.10, 0.27, 0.44, 0.61, 0.79, 0.96, 1.13, 1.30, 1.47, 1.64, 1.81, 1.99, 2.16, 2.33, 2.50]

# outbreak_probabilities_1 = [0.000, 0.000, 0.327, 1.437, 2.110, 2.797, 3.063, 3.363, 3.490, 3.593, 3.857, 3.843, 3.637, 3.337, 2.960] 
# outbreak_probabilities_2 = [0.000, 0.007, 0.633, 2.320, 3.153, 3.583, 4.020, 4.267, 4.470, 4.680, 4.637, 4.790, 4.387, 3.950, 3.520]

# a_1 = 0.00
# a_2 = 0.99
# plt.figure(figsize=(8, 5))
# plt.plot(np.array(betas), outbreak_probabilities_1, marker='o', linestyle='-',  linewidth=4, label=f"Assortativity low ({a_1})")
# plt.plot(np.array(betas), outbreak_probabilities_2, marker='o', linestyle='-', linewidth=4, label=f"Assortativity high ({a_2})")
# #plt.plot(np.array(betas)/params["gamma"], outbreak_probabilities_2, marker='o', linestyle='-', label=f"Assortativity low")

# plt.xlabel(r"$\beta$ (Infection Rate)", fontsize=20)
# plt.ylabel("Adaptation rate", fontsize=20)
# plt.grid(True)
# # plt.ylim([-0.01, 1.02])
# # ax = plt.gca()
# # ax.spines['top'].set_visible(False)
# # ax.spines['right'].set_visible(False)
# plt.title("Adaptation rate vs. Beta for Different Assortativities", fontsize=20)
# #plt.ylim([-0.01, 1.02])
# plt.legend()
# plt.show()


## Plotting the adaptation rate vs. beta for different assortativities for Case 4 Endemic

# betas = [0.10, 0.27, 0.44, 0.61, 0.79, 0.96, 1.13, 1.30, 1.47, 1.64, 1.81, 1.99, 2.16, 2.33, 2.50]
# outbreak_probabilities_1 = [0.000, 0.000, 0.057, 0.150, 0.300, 0.493, 0.577, 0.747, 0.850, 0.900, 1.000, 0.950, 1.000, 1.013, 1.033]
# outbreak_probabilities_2 = [0.000, 0.000, 0.110, 0.270, 0.407, 0.520, 0.607, 0.773, 0.810, 0.917, 0.890, 0.973, 1.030, 1.057, 1.123]
# a_1 = 0.00
# a_2 = 0.99
# # -------------------------
# # Plotting the Results 2 # recall that the sigma parameter controls the strength of immunity (vaccinated/recovered)
# # -------------------------
# plt.figure(figsize=(8, 5))
# plt.plot(np.array(betas), outbreak_probabilities_1, marker='o', linestyle='-', linewidth=4, label=f"Assortativity low ({a_1})")
# plt.plot(np.array(betas), outbreak_probabilities_2, marker='o', linestyle='-', linewidth=4,label=f"Assortativity high ({a_2})")
# #plt.plot(np.array(betas)/params["gamma"], outbreak_probabilities_2, marker='o', linestyle='-', label=f"Assortativity low")

# plt.xlabel(r"$\beta$ (Infection Rate)", fontsize=20)
# plt.ylabel("Adaptation rate", fontsize=20)
# plt.grid(True)
# # plt.ylim([-0.01, 1.02])
# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.title("Adaptation rate vs. Beta for Different Assortativities", fontsize=20)
# #plt.ylim([-0.01, 1.02])
# plt.legend()
# plt.show()


## Phylodyanic Curve Case 4
# 

# outbreak_probabilities_1 = [0.745,0.650,0.720,0.605,0.570,0.535,0.550,0.460,0.365,0.335,0.240,0.185,0.140,0.085,0.235,0.285,0.140,0.075,0.005,0.000]
# outbreak_probabilities_2 = [0.780,0.750,0.695,0.710,0.705,0.655,0.515,0.480,0.470,0.345,0.295,0.305,0.275,0.835,0.635,0.045,0.000,0.000,0.000,0.000]
# sigmas = [0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56, 0.61, 0.66, 0.71, 0.76, 0.81, 0.86, 0.91, 0.96]
# a_1 = 0.00
# a_2 = 0.99

# plt.figure(figsize=(8, 5))
# plt.plot(np.array(sigmas), outbreak_probabilities_1, marker='o', linestyle='-', linewidth=4, label=f"Assortativity low ({a_1})")
# plt.plot(np.array(sigmas), outbreak_probabilities_2, marker='o', linestyle='-', linewidth=4, label=f"Assortativity high ({a_2})")
# #plt.plot(np.array(betas)/params["gamma"], outbreak_probabilities_2, marker='o', linestyle='-', label=f"Assortativity low")
# plt.xlabel(r"$\sigma$ (strength of immunity)", fontsize=20)
# plt.ylabel("Adaptation rate", fontsize=20)
# plt.grid(True)
# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# #plt.ylim([-0.01, 1.02])
# plt.legend()
# plt.show()

## Phylodyanimc Curve Case 2
outbreak_probabilities_1 = [0.840,1.170,1.570,1.725,2.090,2.320,2.590,2.845,2.930,3.345,3.780,3.695,3.530,2.670,1.590,0.580,0.230,0.080,0.015,0.000]
outbreak_probabilities_2 = [0.810,1.065,1.375,1.745,2.210,2.590,2.705,3.265,3.405,3.930,4.420,4.490,4.535,3.960,1.400,0.135,0.015,0.005,0.000,0.000]

# outbreak_probabilities_1 = [0.585,0.750,1.100,1.345,1.605,1.880,2.005,2.095,2.160,2.070,2.010,1.815,1.165,0.610,0.295,0.055,0.010,0.000,0.000,0.000]
# outbreak_probabilities_2 = [0.560,0.760,1.075,1.345,1.725,1.950,2.365,2.435,2.915,2.760,2.910,2.415,1.285,0.285,0.085,0.010,0.010,0.000,0.000,0.000]
sigmas = [0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56, 0.61, 0.66, 0.71, 0.76, 0.81, 0.86, 0.91, 0.96]
a_1 = 0.00
a_2 = 0.99

plt.figure(figsize=(8, 5))
plt.plot(np.array(sigmas), outbreak_probabilities_1, marker='o', linestyle='-', linewidth=4, label=f"Assortativity low ({a_1})")
plt.plot(np.array(sigmas), outbreak_probabilities_2, marker='o', linestyle='-', linewidth=4, label=f"Assortativity high ({a_2})")
#plt.plot(np.array(betas)/params["gamma"], outbreak_probabilities_2, marker='o', linestyle='-', label=f"Assortativity low")
plt.xlabel(r"$\sigma$ (strength of immunity)", fontsize=20)
plt.ylabel("Adaptation rate", fontsize=20)
plt.grid(True)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#plt.ylim([-0.01, 1.02])
plt.legend()
plt.show()

## Philodynaimc Curve Case 1
