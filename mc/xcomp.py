import sys
from collections import Counter
from analyze import tiles
from famtri import famstats
for fn in sys.argv[1:]:
    X,ID,N,polys,lam=tiles(fn); c=Counter(len(p) for p in polys); tot=sum(c.values()); big=sum(v for k,v in c.items() if k>=5)
    pt,Mt,pr,Mr=famstats(fn)
    print(f"{fn:26s} N_t/N_q={c.get(3,0)/max(c.get(4,1),1):5.2f}  defects {100*big/tot:4.1f}%  M_rh={Mr:.2f}")
