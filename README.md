#  AttoFab: 3D Semiconductor Fab Digital Twin


AttoFab is a high-fidelity, edge-optimized 3D Digital Twin of an end-to-end semiconductor fabrication plant. Designed to democratize semiconductor education for India's Semiconductor Mission, AttoFab visualizes the entire manufacturing lifecycle—from bare Silicon (Si) substrate to advanced Fan-Out Wafer-Level Packaging (FOWLP)—entirely on local hardware.

##  Core Features
* **End-to-End 3D Visualization:** Interactive WebGL models demonstrating sequential wafer transformations (Oxidation, Nanolithography, Etching, Metallization, and Packaging).
* **YieldGuard SPC Engine:** Real-time Statistical Process Control dashboard utilizing EWMA and CUSUM methods to monitor manufacturing yield stability across sequential nodes.
* **SECOM ML Defect Classifier:** An integrated intelligence layer utilizing gradient boosting models (XGBoost/LightGBM) on simulated sensor data to predict and classify stochastic defects before they render.
* **AMD Edge-Optimized:** Built to push local multi-threading to its limits. By running React Three Fiber rendering and Python ML inference locally, AttoFab delivers zero-latency simulation without expensive cloud GPU costs, perfectly tailored for AMD Ryzen™ AI processors.

##  Architecture Stack
* **Frontend:** React 18, Three.js (WebGL), Tailwind CSS, Recharts, Axios
* **Backend:** FastAPI, Uvicorn (Python)
* **Physics & Intelligence:** Scikit-Learn, LightGBM, custom Rayleigh Criterion lithography math engines.

---

##  Local Setup & Installation

AttoFab requires two terminal windows to run simultaneously (one for the AI/Physics engine, one for the 3D UI).

### Prerequisites
* Node.js (v16+)
* Python (3.8+)
* Git

### 1. Start the Backend (Atto-Engine)
Open your first terminal and run the following commands to initialize the Python engine on port 8080:

```bash
# Clone the repository
git clone [https://github.com/YOUR_USERNAME/AttoFab.git](https://github.com/YOUR_USERNAME/AttoFab.git)
cd AttoFab/backend

# Create and activate a virtual environment
python -m venv atto_env
# On Windows:
atto_env\Scripts\activate
# On Mac/Linux:
# source atto_env/bin/activate

# Install dependencies
pip install -r requirements.txt

# Start the FastAPI server (Bypassing Windows .exe wrappers)
python -m uvicorn app.main:app --host 127.0.0.1 --port 8080 --reload
