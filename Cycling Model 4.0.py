import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, radians
import scipy.misc
import scipy
from scipy.interpolate import interp1d
import API
from PowerNorm import normPower
from PowerNorm import averagePower
from PowerNorm import cumMovAv
from PowerNorm import exMovAv

courseLength = API.distances[-1] # course length in meters
g = 9.81 # gravity


#vC using 16.6 kmph, 250 W/s, gives a time of 142.2 and NP of 300.092
vCFastest = [(5.356, 300), (10.584, 197.0), (15.441, 149.0), (20.245, 94.43), (25.014, 72.15), (29.515, 149.22), (33.676, 273.68), (37.809, 388.39), (42.152, 425.0), (46.772, 395.42), (51.58, 328.38), (56.479, 247.15), (61.398, 166.95), (66.276, 106.11), (71.047, 79.93), (75.677, 90.35), (80.207, 117.0), (84.727, 137.13), (89.298, 139.6), (93.94, 123.68), (98.602, 119.98), (103.093, 183.1), (107.381, 278.97), (111.69, 340.94), (116.213, 330.96), (121.004, 257.58), (125.951, 162.23), (130.887, 88.86), (135.673, 65.83), (140.254, 97.48), (144.658, 168.66), (148.977, 248.02), (153.326, 308.15), (157.773, 338.97), (162.331, 340.12), (166.996, 314.34), (171.749, 267.67), (176.57, 204.98), (181.447, 128.21), (186.37, 42.78), (191.331, 0), (196.333, 0), (201.314, 0), (206.089, 0.82), (210.441, 115.0), (214.568, 365.0), (218.937, 398.47), (223.5, 399.04), (228.147, 382.11), (232.861, 345.69), (237.663, 282.32), (242.584, 183.48), (247.634, 61.28), (252.646, 0), (257.336, 28.78), (261.673, 124.95), (265.934, 209.56), (270.317, 250.91), (274.861, 251.2), (279.538, 220.78), (284.304, 175.74), (289.034, 161.77), (293.554, 219.06), (297.859, 307.3), (302.256, 332.69), (306.928, 287.12), (311.789, 210.84), (316.669, 148.64), (321.428, 129.94), (326.006, 158.9), (330.46, 205.57), (334.925, 235.26), (339.494, 230.47), (344.202, 189.04), (349.003, 136.42), (353.775, 104.65), (358.438, 105.48), (362.985, 133.75), (367.453, 178.51), (371.877, 232.09), (376.281, 286.55), (380.731, 316.91), (385.324, 298.99), (390.122, 225.61), (395.084, 121.06), (400.095, 18.24), (405.149, 0), (410.323, 0), (415.655, 0), (421.138, 0), (426.771, 0), (432.576, 0), (438.579, 0), (444.799, 0), (451.237, 0), (457.875, 0), (464.643, 0), (471.39, 0), (477.895, 0), (483.925, 0), (489.336, 0), (494.067, 45.05), (498.323, 158.81), (502.581, 230.64), (507.045, 248.09), (511.628, 254.54), (516.157, 292.67), (520.571, 350.87), (524.996, 391.28), (529.497, 414.14), (534.025, 437.18), (538.563, 450.67), (543.187, 432.55), (547.905, 401.77), (552.576, 401.51), (557.122, 429.06), (561.635, 449.56), (566.233, 437.33), (570.981, 380.57), (575.903, 281.43), (580.81, 240.68), (585.201, 490.68), (589.212, 740.68), (593.014, 750), (596.917, 750), (601.309, 728.03), (606.288, 478.03), (611.105, 482.95), (615.449, 595.55), (619.59, 719.87), (623.739, 750), (628.129, 746.71), (632.87, 721.74), (637.512, 707.29), (642.402, 584.49), (647.507, 462.31), (652.318, 477.71), (656.642, 582.13), (660.938, 640.04), (665.491, 623.98), (670.276, 557.4), (675.191, 469.26)]


def dragForce(groundSpeed, CdA, dist):
    airSpeed = findAirspeed(groundSpeed, dist)
    
    airDensity = 1.225
    drag = 0.5 * airDensity * (airSpeed**2) * CdA
    return(drag)


def findAirspeed(groundSpeed, dist):
    windBearing = 90 # 0 to 359 deg (90deg would mean a wind FROM the East)
    windSpeed = 0 # m/s
    windShelter = 1 # multiplier from 0 to 1 - 0 = totally shelterred, 1 = totally exposed
    
    for i in range(len(API.distances)): # find closest bearing to current speed
        if i > 0:
            if API.distances[i-1] <= dist < API.distances[i]:
                bearingIndex = i - 1
                bearing = API.bearings[bearingIndex]
            if dist > API.distances[-1]:
                bearing = API.bearings[-1]
    
    deltaBearing = windBearing - bearing
    effectiveWind = windSpeed * cos(radians(deltaBearing)) * windShelter
    
    airSpeed = groundSpeed + effectiveWind
    
    return(airSpeed)


def rollingResistance(dist, systemMass, Crr): # function to output rolling resistance drag force
    rr = Crr * systemMass * g * math.cos(math.atan(getGrad(dist)))
    return(rr)


def gravityForce(dist, systemMass): # function to output gravity drag force
    gravdrag = systemMass * g * math.sin(math.atan(getGrad(dist)))
    return(gravdrag)


def getGrad(dist): # function to return the current gradient
    if dist <= 0.1:
        grad = (API.segmentdatainterp(dist + 0.1) - API.segmentdatainterp(dist))/0.1
    elif dist >= courseLength - 0.1:
        grad = (API.segmentdatainterp(courseLength) - API.segmentdatainterp(courseLength - 0.1))/0.1
    else:
        grad = scipy.misc.derivative(API.segmentdatainterp, dist, 0.1)
    return(grad)


def getPower(dist, maxPower, pCBias): # function to combine the two power profiles
    
    # pCBias is bias towards pC or vC, 1 = totally pC, 0 = totally vC
    
    if dist == 0:
        vCPower = vCFastest[0][1]
        
                
    for i in range(len(vCFastest)):
        if i < len(vCFastest) - 1:
            if vCFastest[i][0] <= dist < vCFastest[i+1][0]:
                vCPower = vCFastest[i][1]
            if dist > vCFastest[-1][0]:
                vCPower = vCFastest[-1][1]
            if dist < vCFastest[0][0]:
                vCPower = vCFastest[0][1] 

    currentPower = (pCBias * 300) + ((1 - pCBias) * vCPower) +0.4
    
    
    if currentPower < 0:
        currentPower = 0
    if currentPower > maxPower:
        currentPower = maxPower
    return(round(currentPower, 2))


def runSimulation(pCBias, plot):
    distanceCurrent = 0
    courseLength = API.distances[-1]   # course length in meters
    currentTime = 0
    groundSpeed = 5 # starting speed (don't use exactly 0 as then infinite acceleration)

    riderMass = 65		# ride mass in kg
    bikeMass = 9	# bike and equipment in kg
    systemMass = bikeMass + riderMass # calculating total mass
    CdA = 0.30 	# Coefficient of drag Area (baseline figure needs refining)
    Crr = 0.00387
    maxPower = 750

    iteration = 0
    distances = []
    speeds = []
    elevations = []
    gradients = []
    airdrags = []
    rollingdrags = []
    gravitydrags = []
    powers = []
 
    while distanceCurrent < courseLength: # main program loop simulating acceleration second by second
        netForce = (getPower(distanceCurrent, maxPower, pCBias) / groundSpeed) - dragForce(groundSpeed, CdA, distanceCurrent) - rollingResistance(distanceCurrent, systemMass, Crr) - gravityForce(distanceCurrent, systemMass)
        acceleration = netForce / systemMass
        groundSpeed = groundSpeed + (acceleration * 0.1)
        currentTime = currentTime + 0.1
        distanceCurrent = distanceCurrent + (groundSpeed*0.1)
        iteration = iteration + 1
        if (iteration%10) == 0: # writing data every 1 second (10 * 0.1 second interations)
            distances.append(distanceCurrent)
            speeds.append(groundSpeed*3.6) # writing speeds as kmph
            gradients.append(getGrad(distanceCurrent)*100)
            if distanceCurrent < courseLength:
                elevations.append(API.segmentdatainterp(distanceCurrent))
            else:
                elevations.append(API.segmentdatainterp(courseLength))
            airdrags.append(dragForce(groundSpeed, CdA, distanceCurrent))
            rollingdrags.append(rollingResistance(distanceCurrent, systemMass, Crr))
            gravitydrags.append(gravityForce(distanceCurrent, systemMass))
            powers.append(getPower(distanceCurrent, maxPower, pCBias))

    aP = averagePower(cumMovAv(powers))
    nP = normPower(cumMovAv(powers))
    variabilityIndex = aP/nP

    powerstenth = []
    for i in powers:
        powerstenth.append(i/10)

    if plot == 1:
        plt.plot(distances, elevations, label="Elevation (m)", color="C0")
        plt.plot(distances, speeds, label="Speed (km/h)", color="C1")
        plt.plot(distances, gradients, label="Gradient (%)", color="C2")
        plt.plot(distances, powerstenth, label="Power (W/10)", color="C3")

        plt.xlabel("Distance (m)")
        plt.title("Optimal Combination Profile")
        plt.xlim(0, distances[-1])
        plt.grid()
        
    return(nP, round(currentTime, 2), aP)
