julia --machine-file mpi-hostfile -t 40 scripts/evolve.jl \
  -i Q9H7B4 \
  -o data/smyd3 \
  -f ../SPRINT.jl/data/sequences/proteome.fasta \
  -s ../SPRINT.jl/data/hsps/proteome_sprintjl_processed.hsp \
  -t ../SPRINT.jl/data/interactions.txt \
  --convergence_iterations 50
