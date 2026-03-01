from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.services.optics import AttoOptics

app = FastAPI()

# This tells the browser to allow React to communicate with FastAPI
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], 
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/status")
async def get_status():
    return {"engine": "AttoFab v1.0", "status": "online"}

@app.post("/simulate")
async def simulate(tech: str, na: float, k1: float, dose: float):
    optics = AttoOptics()
    return optics.calculate_metrics(tech, na, k1, dose)