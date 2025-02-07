import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Inputs
dose = 5  # mg per dose
half_life = 120  # hours (~5 days)
lambda_decay = np.log(2) / half_life  # Elimination rate constant

# Time array (simulation is for for 50 days in 1-hour steps)
time_hours = np.arange(0, 1200, 1)
concentration = dose * np.exp(-lambda_decay * time_hours)

# Simulate weekly injections
weekly_times = np.arange(0, 1200, 7 * 24)  # Every 7 days (168 hours)
cumulative_concentration = np.zeros_like(time_hours, dtype=float)

for t_dose in weekly_times:
    cumulative_concentration += dose * np.exp(-lambda_decay * (time_hours - t_dose)) * (time_hours >= t_dose)

# Plot the results
plt.figure(figsize=(10, 5))
plt.plot(time_hours / 24, concentration, label="Single Dose", linestyle="dashed")
plt.plot(time_hours / 24, cumulative_concentration, label="Weekly Injections", linewidth=2)
plt.axhline(y=max(cumulative_concentration), color="red", linestyle="dotted", label="Steady-State Max")
plt.xlabel("Time (Days)")
plt.ylabel("Drug Concentration (mg)")
plt.title("GLP-1 Drug Concentration Over Time")
plt.legend()
plt.grid(True)
plt.show()
