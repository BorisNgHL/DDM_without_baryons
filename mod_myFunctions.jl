# using Printf
using QuadGK


######################################### For dmOnly() #############################################
# NFW_params = [rho_0, R_s, c] (see Wiki)
NFW_density(NFW_params, r) = NFW_params[1] / (r / NFW_params[2]) / (1 + r / NFW_params[2]) ^ 2
# NFW_enclosedMass(NFW_params, r) = 4 * pi * NFW_params[1] * NFW_params[2] ^ 3 * (log(1 + r / NFW_params[2]) - r / (NFW_params[2] + r))

# shellRange = [r_1, r_2] where r_1 < r_2
# NFW_shellMass(NFW_params, shellRange) = NFW_enclosedMass(NFW_params, shellRange[2]) - NFW_enclosedMass(NFW_params, shellRange[1])


function NFW_shellMass(NFW_params, shellRange)
    integrand(r) = 4 * pi * r ^ 2 * NFW_density(NFW_params, r)
    
    return quadgk(integrand, shellRange[1], shellRange[2])[1]
end    

 
######################################## For verify_NFW() ##########################################
function NFW_GPE(NFWshells_radii, NFW_params, G)
    NFWshells_GPE = zeros(size(NFWshells_radii, 1))
    for i in 1:size(NFWshells_GPE, 1)
        NFWshells_GPE[i] = -4 * pi * G * NFW_params[1] * NFW_params[2] ^ 3 / NFWshells_radii[i, 3] * log(1 + NFWshells_radii[i, 3] / NFW_params[2])
    end

    return NFWshells_GPE
end


function printToFile_verify_NFW_GPE(fileName, Tshells_radii, Tshells_GPE, NFWshells_GPE)
    f = open(fileName, "w")

    for i in 1:size(Tshells_radii, 1)
        println(f, Tshells_radii[i, 3], "\t", Tshells_GPE[i], "\t", NFWshells_GPE[i])
    end

    return nothing
end


function printToFile_NFW_effPotentialProfile(fileName, Mshells_radii, potentialProfile)
    f = open(fileName, "w")

    for i in 1:size(Mshells_radii, 1)
        println(f, Mshells_radii[i], "\t", potentialProfile[i])
    end

    return nothing
end


####################################################################################################
# Return mass array of NFW profile
# shells_radii = [inner radius, outer radius, shell radius] in the ith row
# shells_mass = [total shell mass]. Assume all mass in a shell concentrate at the position just inside the shell radius
function NFW_shells(NFW_params, numOfShells, shellThicknessFactor)
    NFW_R_vir = NFW_params[2] * NFW_params[3]

    # Exponentially increasing shellThickness
    firstShellThickness = NFW_R_vir * (1 - shellThicknessFactor) / (1 - shellThicknessFactor ^ numOfShells)

    shells_radii = zeros(numOfShells, 3)
    shells_mass = zeros(size(shells_radii, 1))
    for i in 1:size(shells_radii, 1)
        shells_radii[i, 1] = firstShellThickness * (1 - shellThicknessFactor ^ (i - 1)) / (1 - shellThicknessFactor)
        shells_radii[i, 2] = shells_radii[i, 1] + firstShellThickness * shellThicknessFactor ^ (i - 1)
        shells_radii[i, 3] = (shells_radii[i, 1] + shells_radii[i, 2]) / 2
        shells_mass[i] = NFW_shellMass(NFW_params, shells_radii[i, 1:2])
    end

    return shells_radii, shells_mass
end

function totalShells(Ashells_radii, Bshells_radii, Ashells_mass, Bshells_mass)
    len_A = size(Ashells_radii, 1)
    len_B = size(Bshells_radii, 1)
    
    if  len_A > len_B
        Tshells_radii = Ashells_radii
        
        Tshells_mass = zeros(len_A)
        for i in 1:len_B
            Tshells_mass[i] = Ashells_mass[i] + Bshells_mass[i]
        end
        for i in len_B + 1:len_A
            Tshells_mass[i] = Ashells_mass[i]
        end
    elseif len_B > len_A
        Tshells_radii = Bshells_radii

        Tshells_mass = zeros(len_B)
        for i in 1:len_A
            Tshells_mass[i] = Ashells_mass[i] + Bshells_mass[i]
        end
        for i in len_A + 1:len_B
            Tshells_mass[i] = Bshells_mass[i]
        end
    else
        Tshells_radii = Ashells_radii
        Tshells_mass = Ashells_mass + Bshells_mass
    end

    return Tshells_radii, Tshells_mass
end


function enclosedMass(shells_radii, shells_mass)
    shells_enclosedMass = zeros(size(shells_radii, 1))
    for i in 1:size(shells_enclosedMass, 1)
        shells_enclosedMass[i] = sum(shells_mass[1:i])
    end

    return shells_enclosedMass
end


# Return GPE (per mass) array from a mass array
function GPE(shells_radii, shells_mass, shells_enclosedMass, G)
    shells_GPE = zeros(size(shells_radii, 1))
    for i in 1:size(shells_GPE, 1)
        shells_GPE[i] = -G * shells_enclosedMass[i] / shells_radii[i, 3]
        
        if i < size(shells_GPE, 1)
            GPEbyOuterShells = 0
            for j in i + 1:size(shells_GPE, 1)
                GPEbyOuterShells += -G * shells_mass[j] / shells_radii[j, 3]
            end
            
            shells_GPE[i] += GPEbyOuterShells
        end
    end

    return shells_GPE
end


# Return angular momentum (per mass) array
function L(shells_radii, shells_enclosedMass, G)
    shells_L = zeros(size(shells_radii, 1))
    for i in 1:size(shells_L, 1)
        shells_L[i] = (G * shells_enclosedMass[i] * shells_radii[i, 3]) ^ (1 / 2)
    end

    return shells_L
end


# Return total energy (per mass) array of any just-decayed particle at different radii
function totalE_afterDecay(shells_radii, shells_GPE, shells_L, v_k)
    shells_totalE_afterDecay = zeros(size(shells_radii, 1))
    for i in 1:size(shells_totalE_afterDecay, 1)
        shells_totalE_afterDecay[i] = shells_GPE[i] + (shells_L[i] / shells_radii[i, 3]) ^ 2 / 2 + v_k ^ 2 / 2
    end

    return shells_totalE_afterDecay
end


function energyEquation(r, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass)    
    if r <= 0  # Rejected
        return zeros(NaN)  # Error
    elseif r <= Tshells_radii[1, 3]  # r small
        return Tshells_GPE[1] + (L / r) ^ 2 / 2 - totalE_afterDecay
    elseif r >= Tshells_radii[end, 3]  # r big
        return -G * Tshells_enclosedMass[end] / r + (L / r) ^ 2 / 2 - totalE_afterDecay
    else  # r in between; value by interpolation
        radiusIndex = -1  # Just for the definition
        for i in 2:size(Tshells_radii, 1)
            if r < Tshells_radii[i, 3]
                radiusIndex = i
                break
            end
        end
        intervalSlope = (Tshells_GPE[radiusIndex] - Tshells_GPE[radiusIndex - 1]) / (Tshells_radii[radiusIndex, 3] - Tshells_radii[radiusIndex - 1, 3])
        intervalIntercept = Tshells_GPE[radiusIndex] - intervalSlope * Tshells_radii[radiusIndex, 3]
        radiusGPE = intervalSlope * r + intervalIntercept
        return radiusGPE + (L / r) ^ 2 / 2 - totalE_afterDecay
    end
end


# Solve for r_min, r_max of the elliptical orbit of a decayed particle from an original r_0 (one of the shell radii) orbit
function ellipseSolver(r_0, L, totalE_afterDecay, shells_radii, Tshells_radii, Tshells_enclosedMass, Tshells_GPE, G, tol_ellipseGuess)
    # Search in [l1, l2] U [r1, r2] using the bisection method

    firstShellThickness = shells_radii[1, 2]  # To be used as a tolerance

    # Some initial checking
    if energyEquation(r_0, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass) >= 0
        # This should not happen unless GPE/totalE are not updated properly (= 0 occurs when v_k = 0)
        println("ellipseSolver: v_k probably too small; no solvable roots")
        # println(energyEquation(r_0, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass))
        # zeros(NaN)  # Halt program
        return r_0, r_0  # If this happens, radii just stay put (i.e. solution for v_k = 0)
    elseif totalE_afterDecay >= 0  # Escaped
        return -1, -1
    else  # If checking passed
        l2 = r_0
        r1 = r_0
    end
    
    # Setting l1 and r2
    l1 = firstShellThickness
    while energyEquation(l1, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass) <= 0
        l1 /= 2
    end
    r2 = shells_radii[end, 3]
    while energyEquation(r2, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass) <= 0
        r2 *= 2
    end

    # Bisection method
    lastDiff = 0
    while (l2 - l1 > firstShellThickness * tol_ellipseGuess) && (l2 - l1 != lastDiff)
        lastDiff = l2 - l1
        l3 = (l1 + l2) / 2
        energyEquation_value = energyEquation(l3, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass)
        if energyEquation_value < 0
            l2 = l3
        elseif energyEquation_value > 0
            l1 = l3
        else
            l1 = l3
            l2 = l3
        end
    end
    lastDiff = 0
    while (r2 - r1 > firstShellThickness * tol_ellipseGuess) && (r2 - r1 != lastDiff)
        lastDiff = r2 - r1
        r3 = (r2 + r1) / 2
        energyEquation_value = energyEquation(r3, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass)
        if energyEquation_value < 0
            r1 = r3
        elseif energyEquation_value > 0
            r2 = r3
        else
            r1 = r3 
            r2 = r3
        end
    end

    root1 = (l1 + l2) / 2
    root2 = (r1 + r2) / 2
    return root1, root2
end


# Return ellipse array
function ellipseRadii(shells_L, shells_totalE_afterDecay, shells_radii, Tshells_radii, Tshells_enclosedMass, Tshells_GPE, G, tol_ellipseGuess)
    shells_ellipseRadii = zeros(size(shells_radii, 1), 2)
    for i in 1:size(shells_ellipseRadii, 1)
        root1, root2 = ellipseSolver(shells_radii[i, 3], shells_L[i], shells_totalE_afterDecay[i], shells_radii, Tshells_radii, Tshells_enclosedMass, Tshells_GPE, G, tol_ellipseGuess)

        shells_ellipseRadii[i, 1] = root1
        shells_ellipseRadii[i, 2] = root2
    end

    return shells_ellipseRadii
end


function newShellsRadii(shells_radii, shells_ellipseRadii)
    firstShellThickness = shells_radii[1, 2]
    shellThicknessFactor = (shells_radii[2, 2] - shells_radii[2, 1]) / firstShellThickness
    maxEllipseRadius = findmax(shells_ellipseRadii)[1]

    totalLen = 0
    newNumOfShells = 0
    while totalLen <= maxEllipseRadius  # Why not use < instead? Can be justified by shell radii describe the interval [a, b), which is consistent with beginning from 0
        newNumOfShells += 1
        # totalLen += newNumOfShells * firstShellThickness
        totalLen += firstShellThickness * shellThicknessFactor ^ (newNumOfShells - 1)
    end

    newShells_radii = zeros(newNumOfShells, 3)
    for i in 1:size(newShells_radii, 1)
        newShells_radii[i, 1] = firstShellThickness * (1 - shellThicknessFactor ^ (i - 1)) / (1 - shellThicknessFactor)
        newShells_radii[i, 2] = newShells_radii[i, 1] + firstShellThickness * shellThicknessFactor ^ (i - 1)
        newShells_radii[i, 3] = (newShells_radii[i, 1] + newShells_radii[i, 2]) / 2
    end

    return newShells_radii
end


function weightFactorSolver(phi, a, e)
    integrand(theta) = (a * (1 - e ^ 2) / (1 + e * cos(theta))) ^ 2 * theta

    nominator = quadgk(integrand, 0, phi)
    denominator =  quadgk(integrand, 0, pi)
    
    if (nominator[2] / nominator[1] > 0.01 / 100) | (denominator[2] / denominator[1] > 0.01 / 100)  # Accuracy check
        println("weightFactorSolver: absolute error from quadgk too large")
    end

    return nominator[1] / denominator[1]
end


# Return a weightFactor array (weightFactor_r_ref(r_0)) given a r_ref
function weightFactorArray(r_ref, shells_ellipseRadii)
    weightFactor = zeros(size(shells_ellipseRadii, 1))
    for i in 1:size(weightFactor, 1)  # Looping each r_0
        r_max = shells_ellipseRadii[i, 2]
        r_min = shells_ellipseRadii[i, 1]
        
        if r_max == -1 && r_min == -1  # Escaped the whole system
            weightFactor[i] = 0
        elseif r_min > r_ref
            weightFactor[i] = 0
        elseif r_max <= r_ref
            weightFactor[i] = 1
        else
            a = (r_min + r_max) / 2
            e = (r_max / r_min - 1) / (r_max / r_min + 1)

            phi = acos((a * (1 - e ^ 2) / r_ref - 1) / e)
            
            weightFactor[i] = weightFactorSolver(phi, a ,e)
        end
    end

    return weightFactor
end


function updateShellsMass(newShells_radii, shells_ellipseRadii, Mshells_mass, p_undecayed)
    Mshells_decayedMass = Mshells_mass * (1 - p_undecayed)  # To be redistributed
    Mshells_mass *= p_undecayed  # Remaining mass
    
    Dshells_enclosedMass_decayedMass = zeros(size(newShells_radii, 1))
    for i in 1:size(Dshells_enclosedMass_decayedMass, 1)
        weightFactor = weightFactorArray(newShells_radii[i, 2], shells_ellipseRadii)
        Dshells_enclosedMass_decayedMass[i] = sum(Mshells_decayedMass .* weightFactor)
    end

    Dshells_decayedMass = zeros(size(Dshells_enclosedMass_decayedMass, 1))
    if Dshells_decayedMass != []  # If all mothers at all radius escape upon decay
        Dshells_decayedMass[1] = Dshells_enclosedMass_decayedMass[1]
        for i in 2:size(Dshells_decayedMass, 1)
            Dshells_decayedMass[i] = Dshells_enclosedMass_decayedMass[i] - Dshells_enclosedMass_decayedMass[i - 1]
        end
    end

    return Mshells_mass, Dshells_decayedMass
end


function adiabaticExpansion(shells_radii, shells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated) 
    # At this moment:
    # Mshells_radii is short as original
    # Dshells_radii is extended
    # Tshells_radii is short as original
    # Tshells_radii_updated is extended 
    
    # if size(Tshells_enclosedMass, 1) < size(shells_radii, 1)
    #     println("adiabaticExpansion: shell sizes do not match by", size(shells_radii, 1) - size(Tshells_enclosedMass, 1))
    #     for i in 1:size(shells_radii, 1) - size(Tshells_enclosedMass, 1)
    #         push!(Tshells_enclosedMass, Tshells_enclosedMass[end])
    #     end
    # end

    expansionRatios = Tshells_enclosedMass[1:size(shells_radii, 1)] ./ Tshells_enclosedMass_updated[1:size(shells_radii, 1)]

    # # Hotfix for expansion ratio very close to 1 (maybe not)
    # for i in 1:size(expansionRatios, 1)
    #     expansionRatios[i] = round(expansionRatios[i], digits=3)  # I just picked digits=3
    # end
    # To check if it actaully contracts instead of expanding. But this doesn't really matter
    contractionCount = count(i -> (i < 1), expansionRatios)
    if contractionCount > 0
        # println("adiabaticExpansion: expansion ratio smaller than 1, i.e. NOT expanding. Count: ", contractionCount, ", min ratio: ", findmin(expansionRatios)[1])
        # zeros(NaN)  # To cause error, halting the program
    end

    # shells_expandedRadii = shells_radii[:, 3] .* expansionRatios
    shells_expandedRadii = shells_radii[:, 2] .* expansionRatios  # Use 2 or 3? 2

    # To make sure expandedRadii is "monotonic" (never seen useful)
    violationCount = 0
    checkedEntry = -1
    while checkedEntry != size(shells_expandedRadii, 1) - 1
        checkedEntry = -1
        for i in 1:size(shells_expandedRadii, 1) - 1
            if shells_expandedRadii[i] > shells_expandedRadii[i + 1]
                violationCount += 1
                
                eR_1 = shells_expandedRadii[i]
                eR_2 = shells_expandedRadii[i + 1]
                shells_expandedRadii[i] = eR_2
                shells_expandedRadii[i + 1] = eR_1

                break
            else
                checkedEntry = i
            end
        end
    end
    if violationCount > 0
        println("adiabaticExpansion: violationCount = ", violationCount)
    end

    expandedShells_radii = newShellsRadii(shells_radii, shells_expandedRadii)
    expandedShells_mass = zeros(size(expandedShells_radii, 1))
    for i in 1:size(expandedShells_radii, 1)  # This interpolation thing should work if the relation is monotonic. Check total mass after expansion.
        e1 = expandedShells_radii[i, 1]  # Inner radius of expanded shells
        e2 = expandedShells_radii[i, 2]  # Outer radius of expanded shells
        
        e1_smallerThanID = -1
        for j in 1:size(shells_expandedRadii, 1)
            if e1 < shells_expandedRadii[j]
                e1_smallerThanID = j
                break
            end
        end
        e2_smallerThanID = -1
        for j in 1:size(shells_expandedRadii, 1)
            if e2 < shells_expandedRadii[j]
                e2_smallerThanID = j
                break
            end
        end

        if e1_smallerThanID == 1
            m = (shells_radii[e1_smallerThanID, 2] - 0) / (shells_expandedRadii[e1_smallerThanID] - 0)
            c = 0
            r1 = m * e1 + c
        elseif e1_smallerThanID != -1
            m = (shells_radii[e1_smallerThanID, 2] - shells_radii[e1_smallerThanID - 1, 2]) / (shells_expandedRadii[e1_smallerThanID] - shells_expandedRadii[e1_smallerThanID - 1])
            c = shells_radii[e1_smallerThanID, 2] - m * shells_expandedRadii[e1_smallerThanID]
            r1 = m * e1 + c
        else
            r1 = -1  # Should never happen
        end

        if e2_smallerThanID == 1
            m = (shells_radii[e2_smallerThanID, 2] - 0) / (shells_expandedRadii[e2_smallerThanID] - 0)
            c = 0
            r2 = m * e2 + c
        elseif e2_smallerThanID != -1
            m = (shells_radii[e2_smallerThanID, 2] - shells_radii[e2_smallerThanID - 1, 2]) / (shells_expandedRadii[e2_smallerThanID] - shells_expandedRadii[e2_smallerThanID - 1])
            c = shells_radii[e2_smallerThanID, 2] - m * shells_expandedRadii[e2_smallerThanID]
            r2 = m * e2 + c
        else
            r2 = -1  # Will happen once
            # println("adiabaticExpansion: r2 = -1")
        end

        firstShellThickness = shells_radii[1, 2]
        shellThicknessFactor = (shells_radii[2, 2] - shells_radii[2, 1]) / firstShellThickness
        if r1 != -1
            totalLen = 0
            r1_smallerThanID = 0
            while totalLen <= r1
                r1_smallerThanID += 1
                totalLen += firstShellThickness * shellThicknessFactor ^ (r1_smallerThanID - 1)
            end
            if r1_smallerThanID > size(shells_radii, 1)
                println("adiabatic Expansion error: r1 > outermost radius")
                continue  # Hotfix to weird boundary cases
            end
        else
            # println(e1)
            # println(shells_expandedRadii)
            println("adiabaticExpansion error: r1 = -1")  # Prompt error
            continue  # Hotfix to weird boundary cases
        end
        if r2 != -1
            totalLen = 0
            r2_smallerThanID = 0
            while totalLen <= r2
                r2_smallerThanID += 1
                totalLen += firstShellThickness * shellThicknessFactor ^ (r2_smallerThanID - 1)
            end
        else
            r2_smallerThanID = -1  # Special treatment
        end
        
        expandedShells_mass[i] += shells_mass[r1_smallerThanID] * (1 - (r1 ^ 3 - shells_radii[r1_smallerThanID, 1] ^ 3) / (shells_radii[r1_smallerThanID, 2] ^ 3 - shells_radii[r1_smallerThanID, 1] ^ 3))
        if r2_smallerThanID == -1
            expandedShells_mass[i] += shells_mass[end]  # This is why the density is always weird at the end (solved)
            r2_smallerThanID = size(shells_radii, 1)
        else
            expandedShells_mass[i] += shells_mass[r2_smallerThanID] * (1 - (shells_radii[r2_smallerThanID, 2] ^ 3 - r2 ^ 3) / (shells_radii[r2_smallerThanID, 2] ^ 3 - shells_radii[r2_smallerThanID, 1] ^ 3))
        end

        if r1_smallerThanID == r2_smallerThanID
            expandedShells_mass[i] -= shells_mass[r1_smallerThanID]
        elseif r2_smallerThanID - r1_smallerThanID > 1
            expandedShells_mass[i] += sum(shells_mass[r1_smallerThanID + 1:r2_smallerThanID - 1])
        end
    end

    return expandedShells_radii, expandedShells_mass
end


function printToFile(shells_radii, shells_mass, fileName)
    f = open(fileName, "w")

    shells_rho = zeros(size(shells_radii, 1))
    shells_enclosedMass = zeros(size(shells_radii, 1))
    shells_avgRho = zeros(size(shells_radii, 1))
    for i in 1:size(shells_rho, 1)
        shells_rho[i] = shells_mass[i] / (shells_radii[i, 2] ^ 3 - shells_radii[i, 1] ^ 3) / (4 / 3 * pi)
        
        shells_enclosedMass[i] = sum(shells_mass[1:i])
        shells_avgRho[i] = shells_enclosedMass[i] / shells_radii[i, 2] ^ 3 / (4 / 3 * pi)
    end

    for i in 1:size(shells_radii, 1)
        println(f, shells_radii[i, 1], "\t", shells_radii[i, 2], "\t", shells_radii[i, 3], "\t", shells_mass[i], "\t", shells_rho[i], "\t", shells_enclosedMass[i], "\t", shells_avgRho[i])
    end

    return nothing
end


function printToFile_GPE(Tshells_radii, Tshells_GPE, fileName)
    f = open(fileName, "w")

    for i in 1:size(Tshells_radii, 1)
        println(f, Tshells_radii[i, 1], "\t", Tshells_radii[i, 2], "\t", Tshells_radii[i, 3], "\t", Tshells_GPE[i])
    end

    return nothing
end

######################################### Not used yet #############################################
function shellTrimmer(shells_radii, shells_mass)
    numOfZeros = 0
    for i in 0:size(shells_radii, 1) - 1
        if shells_mass[end - i] == 0
            numOfZeros += 1
        else
            break
        end
    end
 
    return shells_radii[1:end - numOfZeros, :], shells_mass[1:end - numOfZeros]
end


# Removing the "Boltzmann tail" of baryon particles (error: T should be r dependent)
function barEscape(T, Tshells_GPE, Bshells_mass, m, k)
    totalBarMass = sum(Bshells_mass)
    
    Bshells_escapeV = zeros(size(Bshells_mass, 1))
    for i in 1:size(Bshells_escapeV, 1)
        Bshells_escapeV[i] = (-Tshells_GPE[i] * 2) ^ (1 / 2) 
    end

    integrand(v) = 4 * pi * v ^ 2 * (m / (2 * pi * k * T)) ^ (3 / 2) * exp(-m * v ^ 2 / (2 * k * T))
    for i in 1:size(Bshells_mass, 1)
        retainedFraction = quadgk(integrand, 0, Bshells_escapeV[i])[1]
        Bshells_mass[i] *= retainedFraction
    end

    totalBarMass_updated = sum(Bshells_mass)

    # println("v_rms = ", (3 * k * T / m) ^ (1 / 2), " kpc / s, max escapeV = ", findmax(Bshells_escapeV)[1], " kpc / s, min escapeV = ", findmin(Bshells_escapeV)[1], " kpc / s")
    println("% escaped: ", (1 - totalBarMass_updated / totalBarMass) * 100, "%")

    return totalBarMass_updated
end


function escapedRemoval(Tshells_enclosedMass, Tshells_GPE_updated, shells_radii, shells_mass, G)
    for i in 1:size(shells_radii, 1)
        KE = G * Tshells_enclosedMass[i] / (2 * shells_radii[i, 3])  # Assume circularly moving particles
        # println(Tshells_GPE_updated[i], "\t", KE)
        if Tshells_GPE_updated[i] + KE > 0
            shells_mass[i] = 0
        end
    end

    shells_radii, shells_mass = shellTrimmer(shells_radii, shells_mass)
    return shells_mass
end


function printToFile_orbitalV(Tshells_radii, Tshells_enclosedMass, G, fileName)
    f = open(fileName, "w")

    for i in 1:size(Tshells_radii, 1)
        println(f, Tshells_radii[i, 1], "\t", Tshells_radii[i, 2], "\t", Tshells_radii[i, 3], "\t", (G * Tshells_enclosedMass[i] / Tshells_radii[i, 2]) ^ 0.5)
    end

    return nothing
end