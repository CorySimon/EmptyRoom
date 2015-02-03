###
#   Author: Cory M. Simon (CoryMSimon@gmail.com)
###
using DataFrames


type Framework
    """
    Stores info of a crystal structure
    """
    structurename::String
    # unit cell sizes (Angstrom)
    a::Float64
    b::Float64
    c::Float64

    # unit cell angles
    alpha::Float64
    beta::Float64
    gamma::Float64

    # volume of the unit cell (Angstrom^3)
    v_unitcell::Float64

    # atom count (in primitive unit cell)
    natoms::Int
    
    # constructor
    Framework() = new()

    # fractional coordinates
    xf::Array{Float64}
    yf::Array{Float64}
    zf::Array{Float64}

    # atom identites
    atoms::Array{String}

    # transformation matrix from fractional to cartesian
    f_to_cartesian_mtrx::Array{Float64}
end

function constructframework(structurename::String; check_coords=true)
    """
    Construct and fill in information on Framework type from a .cssr crystal structure file

    :param String structurename: name of structure. Will try to import file structurename.cssr
    :param Bool check_coords: check that fractional coords are in [0.,1]
    """
    # construct the framework object
    framework = Framework()
    framework.structurename = structurename

    # open crystal structure file
    if ~ isfile("data/structures/" * structurename * ".cssr")
        @printf("Crystal structure file %s not present in data/structures/", structurename * ".cssr")
    end
    f = open("data/structures/" * structurename * ".cssr")

    # get unit cell sizes
    line = split(readline(f)) # first line is (a, b, c)
    framework.a = float(line[1])
    framework.b = float(line[2])
    framework.c = float(line[3])

    # get unit cell angles. Convert to radians
    line = split(readline(f))
    framework.alpha = float(line[1]) * pi / 180.0
    framework.beta = float(line[2]) * pi / 180.0
    framework.gamma = float(line[3]) * pi / 180.0

    # write transformation matrix from fractional to cartesian coords
    v_unit = sqrt(1.0 - cos(framework.alpha)^2 - cos(framework.beta)^2 - cos(framework.gamma)^2 +
            2 * cos(framework.alpha) * cos(framework.beta) * cos(framework.gamma)) # volume of unit parallelpiped
    framework.v_unitcell = v_unit * framework.a * framework.b * framework.c  # volume of unit cell
    framework.f_to_cartesian_mtrx = Array(Float64, (3,3))


    framework.f_to_cartesian_mtrx[1,1] = framework.a
    framework.f_to_cartesian_mtrx[1,2] = framework.b * cos(framework.gamma)
    framework.f_to_cartesian_mtrx[1,3] = framework.c * cos(framework.beta)
    framework.f_to_cartesian_mtrx[2,1] = 0.0
    framework.f_to_cartesian_mtrx[2,2] = framework.b * sin(framework.gamma)
    framework.f_to_cartesian_mtrx[2,3] = framework.c * (cos(framework.alpha) - cos(framework.beta) * cos(framework.gamma)) / sin(framework.gamma)
    framework.f_to_cartesian_mtrx[3,1] = 0.0
    framework.f_to_cartesian_mtrx[3,2] = 0.0
    framework.f_to_cartesian_mtrx[3,3] = framework.c * v_unit / sin(framework.gamma)
    
    # get atom count, initialize arrays holding coords
    framework.natoms = int(split(readline(f))[1])
    framework.atoms = Array(String, framework.natoms)
    framework.xf = zeros(Float64, framework.natoms)  # fractional coordinates
    framework.yf = zeros(Float64, framework.natoms)
    framework.zf = zeros(Float64, framework.natoms)

    # read in atoms and fractional coordinates
    readline(f) # waste a line
    for a = 1:framework.natoms
        line = split(readline(f))

        framework.atoms[a] = line[2]

        framework.xf[a] = float(line[3]) % 1.0 # wrap to [0,1]
        framework.yf[a] = float(line[4]) % 1.0
        framework.zf[a] = float(line[5]) % 1.0
        
        if check_coords # assert fractional coords in [0,1]
            @assert ((framework.xf[a] >= 0.0) & (framework.xf[a] <= 1.0)) "Fraction coords not in [0,1]!\n"
            @assert ((framework.yf[a] >= 0.0) & (framework.yf[a] <= 1.0)) "Fraction coords not in [0,1]!\n"
            @assert ((framework.zf[a] >= 0.0) & (framework.zf[a] <= 1.0)) "Fraction coords not in [0,1]!\n"
        end
    end
    
    close(f) # close file

    return framework
end

function crystaldensity(framework::Framework)
    """
    Computes crystal density of a framework (kg/m3)
    """
    if ! isfile("data/atomicmasses.csv")
        print("Atomic masses file data/atomicmasses.csv not present")
    end
    df = readtable("data/atomicmasses.csv")

    mass = 0.0 # count mass of all atoms here
    for a = 1:framework.natoms
        if ~ (framework.atoms[a] in df[:atom])
            error(@sprintf("Framework atom %s not present in data/atomicmasses.cv file", framework.atoms[a]))
        end
        mass += df[df[:atom] .== framework.atoms[a], :][:mass][1]
    end
    
    return mass / framework.v_unitcell * 1660.53892  # --> kg/m3
end


function put_cssr_coords_in_0_1(frameworkname::String)
    """
    In an incorrect .cssr, reflect coords to [0,1]
    """
    framework = constructframework(frameworkname, check_coords=false)
    new_cssr_file = open("data/structures/" * framework.structurename * "_corrected.cssr", "w")
    @printf(new_cssr_file, "%f %f %f\n%f %f %f\n%d 0\ncorrected by PEviz\n", 
        framework.a, framework.b, framework.c,
        framework.alpha*180.0/pi, framework.beta*180.0/pi, framework.gamma*180.0/pi,
        framework.natoms)
    
    count = 1 
    for a = 1:framework.natoms
        @printf(new_cssr_file, "%d %s %f %f %f\n", 
            count, framework.atoms[a], 
            mod(framework.xf[a], 1.0), mod(framework.yf[a], 1.0),mod(framework.zf[a], 1.0))
        count += 1
    end
    close(new_cssr_file)
    @printf("Wrote correct file with coords in [0,1] in data/structures/%s\n", framework.structurename * "_corrected.cssr")
end

function get_replication_factors(f_to_cartesian_mtrx::Array{Float64}, cutoff::Float64, once_or_twice::String)
    """
    Get replication factors for the unit cell such that only directly adjacent ghost unit cells are required for consideration in applying periodic boundary conditions (PBCs).
    That is, a particle in the "home" unit cell will not be within the cutoff distance of any particles two ghost cells away from the home cell.
    This is done by ensuring that a sphere of radius cutoff can fit inside the unit cell
    
    :param Array{Float64} f_to_cartesian_mtrx: transformation matrix from fractional to Cartesian coordinates
    :param Float64 cutoff: Lennard-Jones cutoff distance, beyond which VdW interactions are approximated as zero.
    :param Bool twice_or_once: bigger than twice the cutoff? otherwise, once (GCMC vs grid)
    """
    if once_or_twice == "twice"
        cutoff = cutoff * 2
    end
    # unit cell vectors
    a = f_to_cartesian_mtrx[:,1]
    b = f_to_cartesian_mtrx[:,2]
    c = f_to_cartesian_mtrx[:,3]

    # normal vectors to the faces of the unit cell
    n_ab = cross(a, b)
    n_bc = cross(b, c)
    n_ac = cross(a, c)

    # vector at center of unit cell. i.e center of sphere that we are ensuring fits within the unit cell
    c0 = [a b c] * [0.5, 0.5, 0.5]
    
    rep_factors = [1, 1, 1] # replication factors stored here. Initialize as 1 x 1 x 1 supercell.
    # replicate unit cell until the distance of c0 to the 3 faces is greater than cutoff/2
    # distance of plane to point: http://mathworld.wolfram.com/Point-PlaneDistance.html
    # a replication.
    while abs(dot(n_bc, c0)) / vecnorm(n_bc) < cutoff 
        rep_factors[1] += 1
        a = 2 * a
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    # b replication
    while abs(dot(n_ac, c0)) / vecnorm(n_ac) < cutoff 
        rep_factors[2] += 1
        b = 2 * b
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    # c replication
    while abs(dot(n_ab, c0)) / vecnorm(n_ab) < cutoff 
        rep_factors[3] += 1
        c = 2 * c
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    
    return rep_factors
end

function write_uc_replication_file(structurename::String, cutoff::Float64, once_or_twice::String)
    """
    Writes unit cell file structurename.uc for Cory's GPU code

    :param String once_or_twice: bigger than once times cutoff or twice? write grid vs gcmc...
    """
    if ! ((once_or_twice == "once") | (once_or_twice == "twice"))
        error("not valid for once_or_twice")
    end
    framework = constructframework(structurename)
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff, once_or_twice)

    f = open("data/uc_replications/" * framework.structurename * ".uc", "w")
    @printf(f, "%d %d %d", rep_factors[1], rep_factors[2], rep_factors[3])
end
