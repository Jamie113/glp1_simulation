import numpy as np
from glp1_model import lambda_decay, dose, half_life

def test_exponential_decay():
    assert np.isclose(dose * np.exp(-lambda_decay * 0), dose), "Decay incorrect at t=0"
    assert np.isclose(dose * np.exp(-lambda_decay * half_life), dose / 2), "Decay incorrect at half-life"

if __name__ == "__main__":
    test_exponential_decay()
    print("✅ All tests passed!")
