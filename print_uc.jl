include("julia_codes/framework.jl")

structures = readdlm("structurelist.txt")
@printf("Writing unit cell replication files with cutoff radius %f A...\n", float64(ARGS[1]))

for k = 1:length(structures)
    write_uc_replication_file(structures[k], float64(ARGS[1]), "once")
    write_uc_replication_file(structures[k], float64(ARGS[1]), "twice")
end
