import React from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area } from 'recharts';

const YieldGraph = ({ dose }) => {
  // Generate a bell curve of data centered around Dose = 100
  const generateData = () => {
    let points = [];
    for (let i = 70; i <= 130; i += 5) {
      const yieldVal = Math.max(0, 99 - (Math.abs(100 - i) * 0.8));
      points.push({ name: i, yield: yieldVal });
    }
    return points;
  };

  const data = generateData();

  return (
    <div className="h-64 w-full bg-black/40 p-4 border border-zinc-800 rounded-lg">
      <p className="text-[10px] text-zinc-500 uppercase mb-4 tracking-widest">Process Window Analysis</p>
      <ResponsiveContainer width="100%" height="100%">
        <AreaChart data={data}>
          <defs>
            <linearGradient id="colorYield" x1="0" y1="0" x2="0" y2="1">
              <stop offset="5%" stopColor="#0ea5e9" stopOpacity={0.3}/>
              <stop offset="95%" stopColor="#0ea5e9" stopOpacity={0}/>
            </linearGradient>
          </defs>
          <CartesianGrid strokeDasharray="3 3" stroke="#1f2937" vertical={false} />
          <XAxis dataKey="name" stroke="#4b5563" fontSize={10} label={{ value: 'Dose (mJ/cm²)', position: 'insideBottom', offset: -5, fontSize: 10, fill: '#4b5563' }} />
          <YAxis stroke="#4b5563" fontSize={10} domain={[60, 100]} />
          <Tooltip 
            contentStyle={{ backgroundColor: '#000', border: '1px solid #0ea5e9', fontSize: '10px' }}
            itemStyle={{ color: '#0ea5e9' }}
          />
          <Area type="monotone" dataKey="yield" stroke="#0ea5e9" fillOpacity={1} fill="url(#colorYield)" />
          {/* Vertical line representing the current user dose */}
          <Line type="monotone" dataKey="yield" stroke="transparent" dot={false} />
        </AreaChart>
      </ResponsiveContainer>
    </div>
  );
};

export default YieldGraph;