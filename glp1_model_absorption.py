# I added the absorption rate for this simulation but the script us crazy ineffiecnet so need to adjust it.
# The glp1_model.py is a better version currently 

import numpy as np
import matplotlib.pyplot as plt


BIOAVAILABILITY = 0.8  # Fraction absorbed
Vd = 10  # Volume of distribution in L
HALF_LIFE_ABSORPTION = 1.5  # Absorption half-life in hours
HALF_LIFE_ELIMINATION = 120  # Elimination half-life in hours (~5 days)
ka = np.log(2) / HALF_LIFE_ABSORPTION  # Absorption rate constant
ke = np.log(2) / HALF_LIFE_ELIMINATION  # Elimination rate constant
TOTAL_HOURS = 1176  # Simulate for 50 days (7 weeks after full titration)
OUTPUT_PATH = "./glp_sim2.png"
FULL_TITRATION_TIME = 12 * 7 * 24
TOTAL_HOURS_POST_TITRATION = FULL_TITRATION_TIME + TOTAL_HOURS

DOSE_TIMES = np.arange(0, TOTAL_HOURS_POST_TITRATION, 7 * 24)

# hourly time steps 
TIME_POINTS = np.arange(0, TOTAL_HOURS_POST_TITRATION, 1) 

# Dose increases every 4 weeks (simulate titrations)
DOSE_SCHEDULE = {
    0: 2.5,  # Weeks 1-4
    4: 5.0,  # Weeks 5-8
    8: 7.5,  # Weeks 9-12
    12: 10.0  # Weeks 13+
}

def get_dose_at_time(t):
    week = int(t // (7 * 24))
    if week < 4:
        return DOSE_SCHEDULE[0]
    elif week < 8:
        return DOSE_SCHEDULE[4]
    elif week < 12:
        return DOSE_SCHEDULE[8]
    else:
        return DOSE_SCHEDULE[12]

# Updated to store individual dose decay contributions
def plasma_concentration(t, dose_times):
    """Compute total plasma concentration at time t more efficiently."""
    total_concentration = 0
    dose_contributions = []

    for t_dose in dose_times:
        if t < t_dose:  # Skip future doses
            break
        dose = get_dose_at_time(t_dose)
        conc = (dose * BIOAVAILABILITY / Vd) * (ka / (ka - ke)) * (
            np.exp(-ke * (t - t_dose)) - np.exp(-ka * (t - t_dose))
        )
        total_concentration += conc
        dose_contributions.append((t_dose, conc))

    return total_concentration, dose_contributions


def main():
    concentration = []
    dose_decay_curves = []

    for t in TIME_POINTS:
        total_conc, dose_contributions = plasma_concentration(t, DOSE_TIMES)
        concentration.append(total_conc)
        dose_decay_curves.append(dose_contributions)

    concentration = np.array(concentration)

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 5))

    # cumulative concentration
    ax.plot(TIME_POINTS / 24, concentration, label="Plasma Concentration (mg/L)", color="b")

    # Individual dose decay curves
    # Individual dose decay curves (dashed gray lines)
    for dose_contrib in dose_decay_curves:
     for t_dose, conc in dose_contrib:
        decay_curve = conc * np.exp(-ke * (TIME_POINTS - t_dose)) * (TIME_POINTS >= t_dose)
        ax.plot(TIME_POINTS[TIME_POINTS >= t_dose] / 24, decay_curve[TIME_POINTS >= t_dose], linestyle="dashed", alpha=0.5, color="gray")

  
    # steady-state max line
    steady_state_max = max(concentration)
    ax.axhline(y=steady_state_max, color="red", linestyle="dotted", label="Steady-State Max")

    # Plottin' in a good way
    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Drug concentration (mg)")
    ax.set_title("GLP-1 concentration over time (Absorption Model)")
    ax.legend()
    fig.tight_layout()

    fig.savefig(OUTPUT_PATH, dpi=300, transparent=True, bbox_inches="tight")
    plt.show()

    return None


if __name__ == "__main__":
    main()
