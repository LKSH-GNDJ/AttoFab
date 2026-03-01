import math

class LithoSimulator:
    def __init__(self):
        # Constants for EUV/DUV
        self.wavelengths = {"EUV": 13.5, "DUV": 193.0}

    def calculate_yield(self, technique: str, na: float, k1: float, dose: float):
        """
        Calculates theoretical Critical Dimension (CD) and Yield Probability.
        Formula: CD = k1 * (lambda / NA)
        """
        wavelength = self.wavelengths.get(technique, 193.0)
        
        # Rayleigh Criterion for Resolution
        critical_dimension = k1 * (wavelength / na)
        
        # Mock logic for Yield based on Dose sensitivity (Stochastic effects)
        # In a real fab, too low a dose leads to photon shot noise (defects)
        base_yield = 98.5
        dose_penalty = abs(100 - dose) * 0.5 # Assuming 100 is optimal dose
        
        final_yield = max(0, min(100, base_yield - dose_penalty))
        
        return {
            "critical_dimension_nm": round(critical_dimension, 2),
            "theoretical_yield": f"{round(final_yield, 2)}%",
            "status": "Success" if final_yield > 90 else "Risk Detected",
            "message": self._generate_feedback(final_yield, critical_dimension)
        }

    def _generate_feedback(self, yield_val, cd):
        if yield_val < 90:
            return f"CD of {cd}nm is achievable, but Dose fluctuations are causing stochastic defects."
        return "Process parameters are within optimal Window."