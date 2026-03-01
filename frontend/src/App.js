import React, { useState } from 'react';
import { Microscope } from 'lucide-react';
import axios from 'axios';
import YieldGraph from './components/YieldGraph';

export default function App() {
  const [dose, setDose] = useState(100);
  const [data, setData] = useState({ cd: 16.36, yield_score: "98.4%", status: "IDLE" });

  const runSim = async () => {
    try {
      // THIS is where that line of code lives! Notice it says 8080 now.
      const res = await axios.post(`http://127.0.0.1:8080/simulate?tech=EUV&na=0.33&k1=0.4&dose=${dose}`);
      setData(res.data);
    } catch (err) {
      console.error("Connection lost to Atto-Engine.");
      setData({ cd: "ERROR", yield_score: "0%", status: "ENGINE_OFFLINE" });
    }
  };

  return (
    <div className="min-h-screen bg-zinc-950 text-sky-400 font-mono p-8">
      <header className="border-b border-sky-900 pb-4 mb-8 flex justify-between items-center">
        <h1 className="text-2xl font-black italic flex items-center gap-2">
          <Microscope className="text-white" /> ATTOFAB
        </h1>
        <div className="text-[10px] text-zinc-500 bg-zinc-900 px-2 py-1 border border-zinc-800 rounded">
          AI ACCELERATED
        </div>
      </header>

      <div className="grid grid-cols-12 gap-8">
        <div className="col-span-4 bg-zinc-900/50 p-6 border border-zinc-800 rounded-lg">
          <label className="text-xs uppercase text-zinc-500 block mb-4">Exposure Dose: {dose} mJ</label>
          <input 
            type="range" min="70" max="130" value={dose} 
            onChange={(e)=>setDose(e.target.value)} 
            className="w-full h-1 bg-zinc-800 rounded-lg appearance-none cursor-pointer mb-8 accent-sky-500" 
          />
          <button onClick={runSim} className="w-full bg-sky-600 text-white py-3 rounded font-bold hover:bg-sky-500 transition">
            EXECUTE_RUN
          </button>
        </div>

        <div className="col-span-8 space-y-6">
          <div className="grid grid-cols-3 gap-4">
            <div className="bg-black border border-zinc-800 p-4 rounded text-center">
              <p className="text-[10px] text-zinc-500">CD (RESIST)</p>
              <p className="text-2xl font-bold text-white">{data.cd} nm</p>
            </div>
            <div className="bg-black border border-zinc-800 p-4 rounded text-center">
              <p className="text-[10px] text-zinc-500">YIELD</p>
              <p className="text-2xl font-bold text-emerald-400">{data.yield_score}</p>
            </div>
            <div className="bg-black border border-zinc-800 p-4 rounded text-center">
              <p className="text-[10px] text-zinc-500">ENGINE STATUS</p>
              <p className="text-[10px] font-bold text-sky-500 mt-2 tracking-tighter">{data.status}</p>
            </div>
          </div>
          <YieldGraph dose={dose} />
        </div>
      </div>
    </div>
  );
}