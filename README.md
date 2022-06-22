#orphan fraction from UniverseMachone

orphan are galaxies which with objid>1e15 (similarly do for satellites) 
While creating catalog for orphan galaxies we select galaxies from clusters massiver than h\*1e14.
We also remove galaxies orphans whose parent is not central galaxy.
How to make catalog:
  select clusters with M>h1e14
  search for galaxies within 'projected' 1h^(-1)Mpc  (we dont need to include h here since cordinates are itself given in hinvMpc)
  Now out of searched galaxies exclude galaxies whose dz from cluster>50h^{-1}Mpc
  write these in a file and remove the galaxies which belong to more than one cluster.
  make sure each cluster has >20 satellites and >1 orphan.
