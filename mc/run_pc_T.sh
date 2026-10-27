#!/bin/bash
# chain A: T=0.20 series ; chain B: T=0.30 series
cd "$(dirname "$0")"
if [ "$1" = A ]; then
  for P in 9.5 10.5 11.5; do python run_npt_cool.py ${START:-../data/configs/fluid_start.npz} $P cT20_P$P 150000 $((100+${P%.*})) 0.40,0.32,0.26,0.22,0.20; done
else
  for P in 10.5 11.5 12.5; do python run_npt_cool.py ${START:-../data/configs/fluid_start.npz} $P cT30_P$P 150000 $((200+${P%.*})) 0.40,0.34,0.30; done
fi
