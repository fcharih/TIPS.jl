/home/fcharih/.julia/bin/mpiexecjl -n 4 julia --project -t 20 scripts/evolve_mpi.jl \
  -p 2000 \
  -i Q9H7B4 \
  -o data/smyd3-3 \
  -f ../SPRINT.jl/data/sequences/proteome.fasta \
  -s ../SPRINT.jl/data/hsps/proteome_sprintjl_processed.hsp \
  -t ../SPRINT.jl/data/interactions.txt \
  --convergence_iterations 50
