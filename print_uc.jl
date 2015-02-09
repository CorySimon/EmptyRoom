include("julia_codes/framework.jl")

structures = readdlm("structurelist.txt")
println("Writing unit cell replication files...")

for k = 1:length(structures)
    write_uc_replication_file(structures[k], 12.5, "once")
    write_uc_replication_file(structures[k], 12.5, "twice")
end
