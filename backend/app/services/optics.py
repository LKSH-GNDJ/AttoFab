import math

class AttoOptics:
    def calculate_metrics(self, tech, na, k1, dose):
        wavelength = 13.5 if tech == "EUV" else 193.0
        cd = k1 * (wavelength / na)
        dof = 0.5 * (wavelength / (na ** 2))
        
        # Yield drops as dose moves away from 100mJ
        yield_val = max(0, 99.0 - (abs(100 - dose) * 0.8))
        
        return {
            "cd": round(cd, 2),
            "dof": round(dof, 2),
            "yield_score": f"{round(yield_val, 1)}%",
            "status": "PROCESS_STABLE" if yield_val > 90 else "STOCHASTIC_RISK"
        }