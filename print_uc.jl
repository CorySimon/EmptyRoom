include("julia_codes/framework.jl")

if (length(ARGS) == 1)
    """
    Get UC reps for all structures in structurelist.txt
    """
    structures = readdlm("structurelist.txt")
    @printf("Writing unit cell replication files with cutoff radius %f A...\n", float64(ARGS[1]))

    for k = 1:length(structures)
        write_uc_replication_file(structures[k], float64(ARGS[1]), "once")
        write_uc_replication_file(structures[k], float64(ARGS[1]), "twice")
    end
end

if (length(ARGS) == 2)
    """
    first arg is cutoff, second arg is structure
    """
    @printf("Writing unit cell replication files for structure %s with cutoff radius %f A...\n", ARGS[2], float64(ARGS[1]))
    write_uc_replication_file(ARGS[2], float64(ARGS[1]), "once")
    write_uc_replication_file(ARGS[2], float64(ARGS[1]), "twice")
end
