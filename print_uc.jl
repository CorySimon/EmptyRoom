include("julia_codes/framework.jl")

if (length(ARGS) == 1)
    """
    Get UC reps for all structures in structurelist.txt
    """
    structures = readdlm("structurelist.txt")
    cutoff = parse(Float64, ARGS[1])
    @printf("\tWriting unit cell replication files with cutoff radius %f A...\n", cutoff)

    for k = 1:length(structures)
        framework = Framework(structures[k])
        framework.check_for_atom_overlap(0.01)
        
        # get rep factors for write grid
        rep_factors = get_replication_factors(framework, cutoff)
        f = open("data/uc_replications/" * framework.structurename * "_once.uc", "w")
        @printf(f, "%d %d %d", rep_factors[1], rep_factors[2], rep_factors[3])
        close(f)

        # get rep factors for gcmc
        rep_factors = get_replication_factors(framework, 2.0 * cutoff)
        f = open("data/uc_replications/" * framework.structurename * "_twice.uc", "w")
        @printf(f, "%d %d %d", rep_factors[1], rep_factors[2], rep_factors[3])
        close(f)
    end
end

if (length(ARGS) == 2)
    """
    first arg is cutoff, second arg is structure
    """
    cutoff = parse(Float64, ARGS[1])
    @printf("\tWriting unit cell replication files for structure %s with cutoff radius %f A...\n", ARGS[2], cutoff)
    framework = Framework(ARGS[2])
    framework.check_for_atom_overlap(0.01)
        
    # get rep factors for write grid
    rep_factors = get_replication_factors(framework, cutoff)
    f = open("data/uc_replications/" * framework.structurename * "_once.uc", "w")
    @printf(f, "%d %d %d", rep_factors[1], rep_factors[2], rep_factors[3])
    close(f)

    # get rep factors for gcmc
    rep_factors = get_replication_factors(framework, 2.0 * cutoff)
    f = open("data/uc_replications/" * framework.structurename * "_twice.uc", "w")
    @printf(f, "%d %d %d", rep_factors[1], rep_factors[2], rep_factors[3])
    close(f)
end
