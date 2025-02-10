import numpy as np
import matplotlib.pyplot as plt

half_life = 120  # hours (~5 days)
lambda_decay = np.log(2) / half_life  # Elimination rate constant
total_hours = 1200  # Simulate for 50 days (~7 weeks after full titration)

# Titration schedule (dose increases every 4 weeks)
titration_schedule = {
    0: 2.5,  # Weeks 1-4
    4: 5.0,  # Weeks 5-8
    8: 7.5,  # Weeks 9-12
    12: 10.0  # Weeks 13+
}

time_hours = np.arange(0, total_hours, 1)

weekly_doses = np.zeros_like(time_hours, dtype=float)

for start_week, dose in titration_schedule.items():
    start_time = start_week * 7 * 24 
    weekly_doses[time_hours >= start_time] = dose 

weekly_times = np.arange(0, total_hours, 7 * 24) 
cumulative_concentration = np.zeros_like(time_hours, dtype=float)

# Store individual dose decays
dose_decay_curves = []

for t_dose in weekly_times:
    dose_at_t = weekly_doses[t_dose]  
    decay_curve = dose_at_t * np.exp(-lambda_decay * (time_hours - t_dose)) * (time_hours >= t_dose)
    
    dose_decay_curves.append(decay_curve)
    
    cumulative_concentration += decay_curve

# Plot the results
plt.figure(figsize=(10, 5))

# Cumulative concentration (actual titration model)
plt.plot(time_hours / 24, cumulative_concentration, label="Weekly Injections with Titration", linewidth=2)

# Individual dose decay curves
for decay_curve in dose_decay_curves:
    plt.plot(time_hours / 24, decay_curve, linestyle="dashed", alpha=0.5, color="gray")

# Steady-state maximum
plt.axhline(y=max(cumulative_concentration), color="red", linestyle="dotted", label="Steady-State Max")

plt.xlabel("Time (Days)")
plt.ylabel("Drug concentration (mg)")
plt.title("GLP-1 concentration over time (with dose titration & decay")
plt.legend()
plt.grid(True)
plt.show()
