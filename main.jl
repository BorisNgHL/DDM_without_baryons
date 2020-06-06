########################################### Information ############################################
# Tested on Julia Version 1.4.0

# Packages reqiured:
# QuadGK
####################################################################################################

include("myFunctions.jl")

######################################## Tunable parameters ########################################
tol_ellipseGuess = 0.01 / 100  # Tolerance for bisection method. In unit of shellThickness
tol_barGuess = 0.01 / 100

const Mo = 1.98847e30  # kg
const kpc = 3.08567758128e19  # m

const G = 6.67430e-11 / kpc ^ 3 * Mo  # m3 kg-1 s-2 to kpc3 Mo-1 s-2
const H = 70 * 1000 / 1000 / kpc  # km s-1 Mpc-1 to s-1
const k = 1.380649e-23 / kpc ^ 2 / Mo  # kg m2 s−2 K-1 to Mo kpc2 s-2 K-1
const rho_c = 3 * H ^ 2 / (8 * pi * G)  # Critical density of the universe

m_molar = 0.75 * 1.00784 + 0.25 * 4.002602;  # Molar mass per particle (75% hydrogen atom + 25% helium atom); g mol-1
m = m_molar / 1000 / 6.02214076e23 / Mo  # Mass per particle (75% hydrogen atom + 25% helium atom); kg to Mo

aIndex = 1.66666667  # Adiabatic index: 5 / 3

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

totalBarMass = 2.10674e9 # totalBarMass to guess for; Mo
initBarRho_0 = 3.015260421975121e6  # Initial baryon core density; Mo kpc-3; typical: 1e9
initCoreT = 560346.2150147454  # Core temperature; K; Milky Way: 1e6
baryonEffectIter = 0  # Number of iterations to update DM (adiabatic expansion / contraction) from changed baryon profile. 2 is good enough

barStopRho = 200 * rho_c  # Typical: 200x (e.g. 103.4 for Xiaoxiong's halo / simulation)
####################################################################################################

######################################## Calculations ##############################################
R_vir = (3 * M_vir / (4 * pi * rho_avg)) ^ (1 / 3)
R_s = R_vir / c
rho_0 = M_vir / (4 * pi * R_s ^ 3) / (log(1 + c) - c / (1 + c))
NFW_params = [rho_0, R_s, c]

K = k * initCoreT / m / (initBarRho_0 ^ (aIndex - 1))  # Polytropic equation's proportionality constant

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
    println(f, "k=", k)
    println(f, "rho_c=", rho_c)
    println(f, "m_molar=", m_molar)
    println(f, "m=", m)
    println(f, "aIndex=", aIndex)
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
    println(f, "initBarRho_0=", initBarRho_0)
    println(f, "initCoreT=", initCoreT)
    println(f, "baryonEffectIter=", baryonEffectIter)
    println(f, "barStopRho=", barStopRho)
    println(f, "R_vir=", R_vir)
    println(f, "R_s=", R_s)
    println(f, "rho_0=", rho_0)
    println(f, "NFW_params=", NFW_params)
    println(f, "K=", K)
    println(f, "initNumOf_M_Shells=", initNumOf_M_Shells)
    println(f, "dt=", dt)

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
        Mshells_mass, Dshells_decayedMass = updateShellsMass(Dshells_decayedRadii, Mshells_ellipseRadii, Mshells_mass, p_undecayed)
        
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

function withBar(totalBarMass)
    functionStart = time_ns()
    stepStart = time_ns()

    folderName = "withBar"
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
    println(f, "k=", k)
    println(f, "rho_c=", rho_c)
    println(f, "m_molar=", m_molar)
    println(f, "m=", m)
    println(f, "aIndex=", aIndex)
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
    println(f, "initBarRho_0=", initBarRho_0)
    println(f, "initCoreT=", initCoreT)
    println(f, "baryonEffectIter=", baryonEffectIter)
    println(f, "barStopRho=", barStopRho)
    println(f, "R_vir=", R_vir)
    println(f, "R_s=", R_s)
    println(f, "rho_0=", rho_0)
    println(f, "NFW_params=", NFW_params)
    println(f, "K=", K)
    println(f, "initNumOf_M_Shells=", initNumOf_M_Shells)
    println(f, "dt=", dt)

    t = 0
    println("Initializing at t=$t Gyr...")

    # Initialize NFW mother shells
    Mshells_radii, Mshells_mass = NFW_shells(NFW_params, initNumOf_M_Shells, shellThicknessFactor)
    # Initialize daughter shells (empty)
    Dshells_radii, Dshells_mass = Mshells_radii, zeros(size(Mshells_radii, 1))
    # Combine mother and daughter shells to get total dark matter mass shells
    TDMshells_radii, TDMshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
    TDMshells_enclosedMass = enclosedMass(TDMshells_radii, TDMshells_mass)

    # Solve for the baryon mass profile
    B_BC, B_params = barConditions(initBarRho_0, K, G, aIndex)
    # Bshells_radii, Bshells_mass, Bshells_radii_hiRes, Bshells_rho_hiRes = barProfile(barStopRho, B_BC, B_params, TDMshells_radii, TDMshells_enclosedMass)
    Bshells_radii, Bshells_mass, newBarRho_0, newCoreT, Bshells_radii_hiRes, Bshells_rho_hiRes = barProfileUpdate(totalBarMass, barStopRho, B_BC, B_params, TDMshells_radii, TDMshells_enclosedMass, tol_barGuess, K, G)
    println("newBarRho_0: ", newBarRho_0)
    println("radius of galaxy: ", Bshells_radii_hiRes[end])
    println("coreT: ", newCoreT)

    ########################### No iteration ever required here: the form of NFW profile should be kept. The appearance of baryons alters angular momentum of the DM particles

    # Combine dark matter and baryons to get total mass shells
    Tshells_radii, Tshells_mass = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
    Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
    Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

    # # Removing the Boltzmann tail (bad)
    # newBarRho_0 = initBarRho_0
    # T = initCoreT
    # totalBarMass = sum(Bshells_mass)
    # for i in 1:100
    #     println("")
    #     Tshells_radii, Tshells_mass = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
    #     Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
    #     Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

    #     totalBarMass_updated = barEscape(T, Tshells_GPE, Bshells_mass, m, k)

    #     B_BC, B_params = barConditions(newBarRho_0, K, G, aIndex)
    #     Bshells_radii, Bshells_mass, newBarRho_0, T, Bshells_radii_hiRes, Bshells_rho_hiRes = barProfileUpdate(totalBarMass_updated, barStopRho, B_BC, B_params, TDMshells_radii, TDMshells_enclosedMass, tol_barGuess, K, G)
    #     totalBarMass = sum(Bshells_mass)
    #     println("newBarRho_0: ", newBarRho_0)
    #     println("radius of galaxy: ", Bshells_radii_hiRes[end])
    # end

    MfileName = folderName * "/M_t=$t.txt"
    printToFile(Mshells_radii, Mshells_mass, MfileName)
    DfileName = folderName * "/D_t=$t.txt"
    printToFile(Dshells_radii, Dshells_mass, DfileName)
    TDMfileName = folderName * "/TDM_t=$t.txt"
    printToFile(TDMshells_radii, TDMshells_mass, TDMfileName)
    BfileName = folderName * "/B_t=$t.txt"
    printToFile(Bshells_radii, Bshells_mass, BfileName)
    BhiResfileName = folderName * "/BhiRes_t=$t.txt"
    printToFile_BhiRes(Bshells_radii_hiRes, Bshells_rho_hiRes, BhiResfileName)
    TfileName = folderName * "/T_t=$t.txt"
    printToFile(Tshells_radii, Tshells_mass, TfileName)
    GPEfileName = folderName * "/GPE_t=$t.txt"
    printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)

    stepResultsFileName = folderName * "/stepResults.txt"
    g = open(stepResultsFileName, "w")

    totalDMmass = sum(TDMshells_mass)
    println("Total DM mass: ", totalDMmass, " Mo")
    totalBarMass = sum(Bshells_mass)
    println("Total baryon mass: ", totalBarMass, " Mo")
    barToDm = totalDMmass / totalBarMass
    println("Baryon to DM ratio: 1 : ", barToDm)
    timeTaken = (time_ns() - stepStart) / 1e9
    println("Time taken for this step: ", timeTaken, "s\n")
    println(g, t, "\t", timeTaken, "\t", totalDMmass, "\t", totalBarMass, "\t", initBarRho_0, "\t", initCoreT, "\t", Bshells_radii_hiRes[end])

    # newBarRho_0 = initBarRho_0
    # newCoreT = initCoreT
    for t in dt:dt:t_end
        stepStart = time_ns()

        println("Working on t=$t Gyr...")
        p_undecayed = exp(log(1 / 2) * t / tau) / exp(log(1 / 2) * (t - dt) / tau)

        # Calculate L and total E of mothers
        Mshells_L = L(Mshells_radii, Tshells_enclosedMass, G)
        Mshells_totalE_afterDecay = totalE_afterDecay(Mshells_radii, Tshells_GPE, Mshells_L, v_k)

        # Get the distribution of decayed mass
        Mshells_ellipseRadii = ellipseRadii(Mshells_L, Mshells_totalE_afterDecay, Mshells_radii, Tshells_radii, Tshells_enclosedMass, Tshells_GPE, G, tol_ellipseGuess)
        Dshells_decayedRadii = newShellsRadii(Dshells_radii, Mshells_ellipseRadii)
        Mshells_mass, Dshells_decayedMass = updateShellsMass(Dshells_decayedRadii, Mshells_ellipseRadii, Mshells_mass, p_undecayed)

        # Prepare enclosed mass for adiabatic expansion
        DandMshells_radii, DandMshells_mass = totalShells(Mshells_radii, Dshells_radii, Mshells_mass, Dshells_mass)
        TDMshells_radii, TDMshells_mass = totalShells(Dshells_decayedRadii, DandMshells_radii, Dshells_decayedMass, DandMshells_mass)
        Tshells_radii_updated, Tshells_mass_updated = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
        Tshells_enclosedMass_updated = enclosedMass(Tshells_radii_updated, Tshells_mass_updated)
        Tshells_GPE_updated = GPE(Tshells_radii_updated, Tshells_mass_updated, Tshells_enclosedMass_updated, G)

        MfileName = folderName * "/M_beforeAdia_t=$t.txt"
        printToFile(Mshells_radii, Mshells_mass, MfileName)
        DfileName = folderName * "/D_beforeAdia_t=$t.txt"
        printToFile(Dshells_radii, Dshells_mass, DfileName)
        DdefileName = folderName * "/Dde_beforeAdia_t=$t.txt"
        printToFile(Dshells_decayedRadii, Dshells_decayedMass, DdefileName)
        TDMfileName = folderName * "/TDM_beforeAdia_t=$t.txt"
        printToFile(TDMshells_radii, TDMshells_mass, TDMfileName)
        BfileName = folderName * "/B_beforeAdia_t=$t.txt"
        printToFile(Bshells_radii, Bshells_mass, BfileName)  # No change so far
        BhiResfileName = folderName * "/BhiRes_beforeAdia_t=$t.txt"
        printToFile_BhiRes(Bshells_radii_hiRes, Bshells_rho_hiRes, BhiResfileName)  # No change so far
        TfileName = folderName * "/T_beforeAdia_t=$t.txt"
        printToFile(Tshells_radii_updated, Tshells_mass_updated, TfileName)
        GPEfileName = folderName * "/GPE_beforeAdia_t=$t.txt"
        printToFile_GPE(Tshells_radii_updated, Tshells_GPE_updated, GPEfileName)

        # Adiabatic expansion of DM
        Mshells_radii, Mshells_mass = adiabaticExpansion(Mshells_radii, Mshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
        # Dshells_radii, Dshells_mass = adiabaticExpansion(Dshells_radii, Dshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
        
        # Job of differentiating between D and Dde is done. D' = D + Dde
        Dshells_radii, Dshells_mass = totalShells(Dshells_radii, Dshells_decayedRadii, Dshells_mass, Dshells_decayedMass)

        # Prepare latest enclosed DM mass for baryon update
        TDMshells_radii, TDMshells_mass = totalShells(Dshells_radii, Mshells_radii, Dshells_mass, Mshells_mass)
        TDMshells_enclosedMass = enclosedMass(TDMshells_radii, TDMshells_mass)
        
        # Prepare latest enclosed mass for potential DM adiabatic expansion from baryon effect
        Tshells_radii, Tshells_mass = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
        Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
        Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

        # Implement the escape mechanism of some baryon mass as a part of the response to the changed GPE (bad)
        # totalBarMass = barEscape(newCoreT, Tshells_GPE, Bshells_mass, m, k)
        # Update baryon given a total baryon mass
        B_BC, B_params = barConditions(newBarRho_0, K, G, aIndex)
        Bshells_radii, Bshells_mass, newBarRho_0, newCoreT, Bshells_radii_hiRes, Bshells_rho_hiRes = barProfileUpdate(totalBarMass, barStopRho, B_BC, B_params, TDMshells_radii, TDMshells_enclosedMass, tol_barGuess, K, G)
        totalBarMass = sum(Bshells_mass)
        println("newBarRho_0: ", newBarRho_0)
        println("radius of galaxy: ", Bshells_radii_hiRes[end])

        # Prepared updated enclosed mass for potential DM adiabatic expansion from baryon effect
        Tshells_radii_updated, Tshells_mass_updated = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
        Tshells_enclosedMass_updated = enclosedMass(Tshells_radii_updated, Tshells_mass_updated)
        Tshells_GPE_updated = GPE(Tshells_radii_updated, Tshells_mass_updated, Tshells_enclosedMass_updated, G)

        ########################### Iteration goes here
        for i in 1:baryonEffectIter
            # DM adiabatic
            Mshells_radii, Mshells_mass = adiabaticExpansion(Mshells_radii, Mshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
            # Dshells_radii, Dshells_mass = adiabaticExpansion(Dshells_radii, Dshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
            
            TDMshells_radii, TDMshells_mass = totalShells(Dshells_radii, Mshells_radii, Dshells_mass, Mshells_mass)
            TDMshells_enclosedMass = enclosedMass(TDMshells_radii, TDMshells_mass)
            Tshells_radii, Tshells_mass = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
            Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
            Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

            # Solve bar
            B_BC, B_params = barConditions(newBarRho_0, K, G, aIndex)
            Bshells_radii, Bshells_mass, newBarRho_0, newCoreT, Bshells_radii_hiRes, Bshells_rho_hiRes = barProfileUpdate(totalBarMass, barStopRho, B_BC, B_params, TDMshells_radii, TDMshells_enclosedMass, tol_barGuess, K, G)
            totalBarMass = sum(Bshells_mass)
            println("newBarRho_0: ", newBarRho_0)
            println("radius of galaxy: ", Bshells_radii_hiRes[end])

            Tshells_radii_updated, Tshells_mass_updated = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
            Tshells_enclosedMass_updated = enclosedMass(Tshells_radii_updated, Tshells_mass_updated)
            Tshells_GPE_updated = GPE(Tshells_radii_updated, Tshells_mass_updated, Tshells_enclosedMass_updated, G)
        end
        ########################### Iteration ends

        # No harm to make sure
        TDMshells_radii, TDMshells_mass = totalShells(Dshells_radii, Mshells_radii, Dshells_mass, Mshells_mass)
        Tshells_radii, Tshells_mass = totalShells(TDMshells_radii, Bshells_radii, TDMshells_mass, Bshells_mass)
        Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
        Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

        MfileName = folderName * "/M_t=$t.txt"
        printToFile(Mshells_radii, Mshells_mass, MfileName)
        DfileName = folderName * "/D_t=$t.txt"
        printToFile(Dshells_radii, Dshells_mass, DfileName)
        TDMfileName = folderName * "/TDM_t=$t.txt"
        printToFile(TDMshells_radii, TDMshells_mass, TDMfileName)
        BfileName = folderName * "/B_t=$t.txt"
        printToFile(Bshells_radii, Bshells_mass, BfileName)
        BhiResfileName = folderName * "/BhiRes_t=$t.txt"
        printToFile_BhiRes(Bshells_radii_hiRes, Bshells_rho_hiRes, BhiResfileName)
        TfileName = folderName * "/T_t=$t.txt"
        printToFile(Tshells_radii, Tshells_mass, TfileName)
        GPEfileName = folderName * "/GPE_t=$t.txt"
        printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)

        totalDMmass = sum(TDMshells_mass)
        println("Total DM mass: ", totalDMmass, " Mo")
        totalBarMass = sum(Bshells_mass)
        println("Total baryon mass: ", totalBarMass, " Mo")
        timeTaken = (time_ns() - stepStart) / 1e9
        println("Time taken for this step: ", timeTaken, "s\n")
        println(g, t, "\t", timeTaken, "\t", totalDMmass, "\t", totalBarMass, "\t", newBarRho_0, "\t", newCoreT, "\t", Bshells_radii_hiRes[end])
    end

        MfileName = folderName * "/M_result.txt"
        printToFile(Mshells_radii, Mshells_mass, MfileName)
        DfileName = folderName * "/D_result.txt"
        printToFile(Dshells_radii, Dshells_mass, DfileName)
        TDMfileName = folderName * "/TDM_result.txt"
        printToFile(TDMshells_radii, TDMshells_mass, TDMfileName)
        BfileName = folderName * "/B_result.txt"
        printToFile(Bshells_radii, Bshells_mass, BfileName)
        BhiResfileName = folderName * "/BhiRes_result.txt"
        printToFile_BhiRes(Bshells_radii_hiRes, Bshells_rho_hiRes, BhiResfileName)
        TfileName = folderName * "/T_result.txt"
        printToFile(Tshells_radii, Tshells_mass, TfileName)
        GPEfileName = folderName * "/GPE_result.txt"
        printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)

        timeTaken_total = (time_ns() - functionStart) / 1e9
        println(f, "timeTaken_total=", timeTaken_total)
        println("Total time taken: ", timeTaken_total, "s\n")

    return nothing
end

# Uncomment to pick which to run
# dmOnly()
# verify_NFW()
# withBar(totalBarMass)