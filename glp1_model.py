import numpy as np
import matplotlib.pyplot as plt


HALF_LIFE = 120  # hours (~5 days)
LAMBDA_DECAY = np.log(2) / HALF_LIFE  # Elimination rate constant
TOTAL_HOURS = 1176  # Simulate for 50 days (7 weeks after full titration)
OUTPUT_PATH = "./glp_sim.png"
FULL_TITRATION_TIME = 12*7*24
TOTAL_HOURS_POST_TITRATION = FULL_TITRATION_TIME + TOTAL_HOURS

# Titration schedule (dose increases every 4 weeks)
TITRATION_SCHEDULE = {
    0: 2.5,  # Weeks 1-4
    4: 5.0,  # Weeks 5-8
    8: 7.5,  # Weeks 9-12
    12: 10.0  # Weeks 13+
}

def main():
    time_hours = np.arange(0, TOTAL_HOURS_POST_TITRATION, 1)
    weekly_doses = np.zeros_like(time_hours, dtype=float)

    for start_week, dose in sorted(TITRATION_SCHEDULE.items()):
        start_time = start_week * 7 * 24 
        weekly_doses[time_hours >= start_time] = dose 

    weekly_times = np.arange(0, TOTAL_HOURS_POST_TITRATION, 7 * 24) 
    cumulative_concentration = np.zeros_like(time_hours, dtype=float)

    # Store individual dose decays
    dose_decay_curves = []
    diffs = []

    for t_dose in weekly_times:
        dose_at_t = weekly_doses[t_dose]
        decay_curve = dose_at_t * np.exp(-LAMBDA_DECAY * (time_hours - t_dose)) * (time_hours >= t_dose)
        
        diff=np.zeros_like(cumulative_concentration)
        for dose_decay in dose_decay_curves:
            diff[t_dose:] += dose_decay[t_dose:]
        diff = diff+decay_curve
        diffs.append(diff)

            
        dose_decay_curves.append(decay_curve)
        
        cumulative_concentration += decay_curve
    

    # Plot the results
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Cumulative concentration (actual titration model)
    ax.plot(time_hours / 24, cumulative_concentration, label="Weekly Injections with Titration", linewidth=2)

    # Individual dose decay curves
    for decay_curve in diffs:
        ax.plot(time_hours / 24, decay_curve, linestyle="dashed", alpha=0.5, color="gray")

    # Steady-state maximum
    ax.axhline(y=max(cumulative_concentration), color="red", linestyle="dotted", label="Steady-State Max")

    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Drug concentration (mg)")
    ax.set_title("GLP-1 concentration over time (with dose titration & decay)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{OUTPUT_PATH}", dpi=300, transparent=True, bbox_inches="tight")
    plt.show()

    return None


if __name__ == "__main__":
    main()
