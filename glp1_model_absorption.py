import numpy as np
import matplotlib.pyplot as plt


BIOAVAILABILITY = 0.8  # Fraction absorbed
Vd = 10  # Volume of distribution in L
HALF_LIFE_ABSORPTION = 1.5
HALF_LIFE_ELIMINATION = 120  # hours (~5 days)
ka = np.log(2) / HALF_LIFE_ABSORPTION  # Absorption rate constant
ke = np.log(2) / HALF_LIFE_ELIMINATION  # Elimination rate constant
TOTAL_HOURS = 1176  # Simulate for 50 days (7 weeks after full titration)
OUTPUT_PATH = "./glp_sim2.png"
FULL_TITRATION_TIME = 12*7*24
TOTAL_HOURS_POST_TITRATION = FULL_TITRATION_TIME + TOTAL_HOURS


DOSE_TIMES = np.arange(0, TOTAL_HOURS_POST_TITRATION, 7*24)
NUM_POINTS = 1000
TIME_POINTS = np.linspace(0, TOTAL_HOURS_POST_TITRATION, NUM_POINTS)
CONCENTRATION = np.zeros_like(TIME_POINTS)

# Titration schedule (dose increases every 4 weeks)
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

def plasma_concentration(t, dose_times):
    total_concentration = 0
    for t_dose in dose_times:
        dose = get_dose_at_time(t_dose)
        if t >= t_dose:
            total_concentration += (dose * BIOAVAILABILITY / Vd) * (ka / (ka - ke)) * (
                np.exp(-ke * (t - t_dose)) - np.exp(-ka * (t - t_dose))
            )
    return total_concentration

def main():
    concentration = np.array([plasma_concentration(t, DOSE_TIMES) for t in TIME_POINTS])

    # Plot the results
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Cumulative concentration (actual titration model)
    ax.plot(TIME_POINTS / 24, concentration, label='Plasma Concentration (mg/mL)', color='b')

    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Drug concentration (mg / mL)")
    ax.set_title("Monjoura (Tirzepatide) Plasma Levels - Titration Regimen")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{OUTPUT_PATH}", dpi=300, transparent=True, bbox_inches="tight")
    plt.show()

    return None


if __name__ == "__main__":
    main()
