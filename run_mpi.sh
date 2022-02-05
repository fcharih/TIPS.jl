/home/fcharih/.julia/bin/mpiexecjl -n 2 julia --project -t 6 scripts/evolve_mpi.jl \
  -p 2000 \
  -i Q9H7B4 \
  -o data/smyd3-3 \
  -f data/proteome.fasta \
  -s data/proteome_sprintjl_processed.hsp \
  -t data/interactions.txt \
  --convergence_iterations 50
