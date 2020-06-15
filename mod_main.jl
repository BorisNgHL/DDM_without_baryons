########################################### Information ############################################
# Tested on Julia Version 1.4.0

# Packages reqiured:
# QuadGK
####################################################################################################

include("mod_myFunctions.jl")

######################################## Tunable parameters ########################################
tol_ellipseGuess = 0.01 / 100  # Tolerance for bisection method. In unit of shellThickness
tol_barGuess = 0.01 / 100

const Mo = 1.98847e30  # kg
const kpc = 3.08567758128e19  # m

const G = 6.67430e-11 / kpc ^ 3 * Mo  # m3 kg-1 s-2 to kpc3 Mo-1 s-2
const H = 70 * 1000 / 1000 / kpc  # km s-1 Mpc-1 to s-1
const rho_c = 3 * H ^ 2 / (8 * pi * G)  # Critical density of the universe

# # Parameters for Kim/XiaoXiong's halo
# c = 23.6  # Concentration parameter
# M_vir = 0.505e10  # Mo
# rho_avg = 103.4 * rho_c

# Parameters for halo
c = 6.24
M_vir = 7.41e10  # Mo
rho_avg = 200 * rho_c

v_k_in_kms = 56
v_k = v_k_in_kms * 1000 / kpc # Recoil velocity of daughter particles; km s-1 to kpc s-1
tau = 2.9  # Half-life of mother particles; Gyr

t_end = 14  # Age of the universe; Gyr
numOfSteps = 20  # 40+ is good enough

firstShellThickness = 1e-3  # Use 1e-n to see from 1e-(n-3) (a bit conservative)
shellThicknessFactor = 1.02  # Thickness of shell grows exponentially by this factor

####################################################################################################

######################################## Calculations ##############################################
R_vir = (3 * M_vir / (4 * pi * rho_avg)) ^ (1 / 3)
R_s = R_vir / c
rho_0 = M_vir / (4 * pi * R_s ^ 3) / (log(1 + c) - c / (1 + c))
NFW_params = [rho_0, R_s, c]

initNumOf_M_Shells = floor(Int, log(1 - NFW_params[2] * NFW_params[3] / firstShellThickness * (1 - shellThicknessFactor)) / log(shellThicknessFactor)) + 1   # Determines initial shellThickness
println("initNumOf_M_Shells: ", initNumOf_M_Shells, "\n")

dt = t_end / numOfSteps
####################################################################################################

########################################## Algorithm ###############################################
function dmOnly()
    functionStart = time_ns()
    stepStart = time_ns()

    folderName = "dmOnly"
    if !isdir(folderName)
        mkdir(folderName)
    end
    folderName = folderName * "/" * string(v_k_in_kms) * "_" * string(tau) * "_" * string(initNumOf_M_Shells) * "_" * string(numOfSteps)
    if !isdir(folderName)
        mkdir(folderName)
    end

    # Print parameters to a file
    paramsFileName = folderName * "/params.txt"
    f = open(paramsFileName, "w")
    
    println(f, "tol_ellipseGuess=", tol_ellipseGuess)
    println(f, "tol_barGuess=", tol_barGuess)
    println(f, "Mo=", Mo)
    println(f, "kpc=", kpc)
    println(f, "G=", G)
    println(f, "H=", H)
    println(f, "rho_c=", rho_c)
    println(f, "c=", c)
    println(f, "M_vir=", M_vir)
    println(f, "rho_avg=", rho_avg)    
    println(f, "v_k_in_kms=", v_k_in_kms)
    println(f, "v_k=", v_k)
    println(f, "tau=", tau)
    println(f, "t_end=", t_end)
    println(f, "numOfSteps=", numOfSteps)
    println(f, "firstShellThickness=", firstShellThickness)
    println(f, "shellThicknessFactor=", shellThicknessFactor)
    println(f, "R_vir=", R_vir)
    println(f, "R_s=", R_s)
    println(f, "rho_0=", rho_0)
    println(f, "NFW_params=", NFW_params)
    println(f, "initNumOf_M_Shells=", initNumOf_M_Shells)
    println(f, "dt=", dt)
    close(f)
    t = 0
    println("Initializing at t=$t Gyr...")

    # Initialize NFW mother shells
    Mshells_radii, Mshells_mass = NFW_shells(NFW_params, initNumOf_M_Shells, shellThicknessFactor)
    # Initialize daughter shells (empty)
    Dshells_radii, Dshells_mass = Mshells_radii, zeros(size(Mshells_radii, 1))
    # Combine mother and daughter shells to get total mass shells
    Tshells_radii, Tshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
    Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
    Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

    MfileName = folderName * "/M_t=$t.txt"
    printToFile(Mshells_radii, Mshells_mass, MfileName)
    DfileName = folderName * "/D_t=$t.txt"
    printToFile(Dshells_radii, Dshells_mass, DfileName)
    TfileName = folderName * "/T_t=$t.txt"
    printToFile(Tshells_radii, Tshells_mass, TfileName)
    GPEfileName = folderName * "/GPE_t=$t.txt"
    printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)

    stepResultsFileName = folderName * "/stepResults.txt"
    g = open(stepResultsFileName, "w")

    totalDMmass = sum(Tshells_mass)
    println("Total DM mass: ", totalDMmass, " Mo")
    timeTaken = (time_ns() - stepStart) / 1e9
    println("Time taken for this step: ", timeTaken, "s\n")
    println(g, t, "\t", timeTaken, "\t", totalDMmass)

    # Rolling starts
    for t in dt:dt:t_end
        stepStart = time_ns()

        println("Working on t=$t Gyr...")
        p_undecayed = exp(log(1 / 2) * t / tau) / exp(log(1 / 2) * (t - dt) / tau)
        
        # Calculate L and total E of mother from the total mass distribution
        Mshells_L = L(Mshells_radii, Tshells_enclosedMass, G)
        Mshells_totalE_afterDecay = totalE_afterDecay(Mshells_radii, Tshells_GPE, Mshells_L, v_k)

        # Solve equation to get ellipse
        Mshells_ellipseRadii = ellipseRadii(Mshells_L, Mshells_totalE_afterDecay, Mshells_radii, Tshells_radii, Tshells_enclosedMass, Tshells_GPE, G, tol_ellipseGuess)
        # Compute the bigger array to contain the new radii
        Dshells_decayedRadii = newShellsRadii(Dshells_radii, Mshells_ellipseRadii)
        # Decay the mothers in the shells, distribute the new daughters
        Mshells_mass, Dshells_decayedMass = updateShellsMass(Dshells_decayedRadii, Mshells_ellipseRadii, Mshells_mass, p_undecayed, Mshells_L, Mshells_totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass)
        
        # # For checking how much of the daughters escaped
        # println(size(Mshells_ellipseRadii, 1))
        # println(count(i -> (i < 0), Mshells_ellipseRadii) / 2)
        # println(size(Dshells_decayedRadii, 1))

        # Now we have: 
        # 1. Mshells (remaining mothers)
        # 2. Dshells (daughters from before)
        # 3. Dshells_decayed (new daughters just decayed) 
        # Two kinds of Dshells because only the old daughters (Dshells) should be expanded

        # Prepare total enclosed mass values for adiabatic expansion
        DandMshells_radii, DandMshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
        Tshells_radii_updated, Tshells_mass_updated = totalShells(DandMshells_radii, Dshells_decayedRadii, DandMshells_mass, Dshells_decayedMass)
        Tshells_enclosedMass_updated = enclosedMass(Tshells_radii_updated, Tshells_mass_updated)
        Tshells_GPE_updated = GPE(Tshells_radii_updated, Tshells_mass_updated, Tshells_enclosedMass_updated, G)  # Just for printing

        MfileName = folderName * "/M_beforeAdia_t=$t.txt"
        printToFile(Mshells_radii, Mshells_mass, MfileName)
        DfileName = folderName * "/D_beforeAdia_t=$t.txt"
        printToFile(Dshells_radii, Dshells_mass, DfileName)
        DdefileName = folderName * "/Dde_beforeAdia_t=$t.txt"
        printToFile(Dshells_decayedRadii, Dshells_decayedMass, DdefileName)
        TfileName = folderName * "/T_beforeAdia_t=$t.txt"
        printToFile(Tshells_radii_updated, Tshells_mass_updated, TfileName)
        GPEfileName = folderName * "/GPE_beforeAdia_t=$t.txt"
        printToFile_GPE(Tshells_radii_updated, Tshells_GPE_updated, GPEfileName)

        # Adiabatic expansions
        Mshells_radii, Mshells_mass = adiabaticExpansion(Mshells_radii, Mshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
        # Dshells_radii, Dshells_mass = adiabaticExpansion(Dshells_radii, Dshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
        
        # # Test for adiabatic convergence (bad)
        # adiaCon_numOfLoops = 10
        # if adiaCon_numOfLoops != 0
        #     Tshells_radii, Tshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
        #     Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
        #     TfileName = folderName * "/T_t=$t.adiaCon=0.txt"
        #     printToFile(Tshells_radii, Tshells_mass, TfileName)
        # end
        # Tshells_enclosedMass_new = Tshells_enclosedMass_updated
        # for i in 1:adiaCon_numOfLoops
        #     Tshells_enclosedMass_old = Tshells_enclosedMass_new
        #     Tshells_enclosedMass_new = Tshells_enclosedMass

        #     Mshells_radii, Mshells_mass = adiabaticExpansion(Mshells_radii, Mshells_mass, Tshells_enclosedMass_old, Tshells_enclosedMass_new)
        #     # Dshells_radii, Dshells_mass = adiabaticExpansion(Dshells_radii, Dshells_mass, Tshells_enclosedMass_old, Tshells_enclosedMass_new)

        #     Tshells_radii, Tshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
        #     Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)

        #     TfileName = folderName * "/T_t=$t.adiaCon=$i.txt"
        #     printToFile(Tshells_radii, Tshells_mass, TfileName)
        # end

        # Update total Dshells
        Dshells_radii, Dshells_mass = totalShells(Dshells_radii, Dshells_decayedRadii, Dshells_mass, Dshells_decayedMass)
        # Update total mass shells
        Tshells_radii, Tshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
        Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
        Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

        MfileName = folderName * "/M_t=$t.txt"
        printToFile(Mshells_radii, Mshells_mass, MfileName)
        DfileName = folderName * "/D_t=$t.txt"
        printToFile(Dshells_radii, Dshells_mass, DfileName)
        TfileName = folderName * "/T_t=$t.txt"
        printToFile(Tshells_radii, Tshells_mass, TfileName)
        GPEfileName = folderName * "/GPE_t=$t.txt"
        printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)

        totalDMmass = sum(Tshells_mass)
        println("Total DM mass: ", totalDMmass, " Mo")
        timeTaken = (time_ns() - stepStart) / 1e9
        println("Time taken for this step: ", timeTaken, "s\n")
        println(g, t, "\t", timeTaken, "\t", totalDMmass)
    end

        MfileName = folderName * "/M_result.txt"
        printToFile(Mshells_radii, Mshells_mass, MfileName)
        DfileName = folderName * "/D_result.txt"
        printToFile(Dshells_radii, Dshells_mass, DfileName)
        TfileName = folderName * "/T_result.txt"
        printToFile(Tshells_radii, Tshells_mass, TfileName)
        GPEfileName = folderName * "/GPE_result.txt"
        printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)

        timeTaken_total = (time_ns() - functionStart) / 1e9
        println(f, "timeTaken_total=", timeTaken_total)
        println("Total time taken: ", timeTaken_total, "s\n")

    return nothing
end

function verify_NFW()
    # Verify my potential function
    Mshells_radii, Mshells_mass = NFW_shells(NFW_params, initNumOf_M_Shells, shellThicknessFactor)
    Dshells_radii, Dshells_mass = Mshells_radii, zeros(size(Mshells_radii, 1))
    Tshells_radii, Tshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
    Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
    
    Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)
    NFWshells_GPE = NFW_GPE(Tshells_radii, NFW_params, G)
    fileName = "verify_NFW_GPE_" * string(initNumOf_M_Shells) * ".txt"
    printToFile_verify_NFW_GPE(fileName, Tshells_radii, Tshells_GPE, NFWshells_GPE)

    # Inspect the effective potential profile for a given L
    Mshells_L = L(Mshells_radii, Tshells_enclosedMass, G)
    Mshells_totalE_afterDecay = totalE_afterDecay(Mshells_radii, Tshells_GPE, Mshells_L, v_k)
    Mshells_radii = Mshells_radii[:, 3]  # Just the shell radii
    
    potentialProfile = zeros(size(Mshells_radii, 1))
    for i in 1:size(potentialProfile, 1)
        # Pick some L at small r: floor(Int, size(Mshells_radii, 1) * 1 / 3)
        potentialProfile[i] = energyEquation(Mshells_radii[i], Mshells_L[floor(Int, size(Mshells_radii, 1) * 1 / 3)], 0, Tshells_radii, Tshells_GPE, Tshells_enclosedMass)
    end
    fileName = "NFW_EffPotentialProfile.txt"
    printToFile_NFW_effPotentialProfile(fileName, Mshells_radii, potentialProfile)

    return nothing
end

# Uncomment to pick which to run
 dmOnly()
# verify_NFW()
# withBar(totalBarMass)