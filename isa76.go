//
//  isa76.swift
//  uconverter
//
//  Created by Charles Mathieu on 11/16/15.
//  Copyright © 2015 Charles.Mathieu. All rights reserved.
//
package main
//import Foundation

import (
    "fmt"
)

const (
    R_SGC = 287.05               // specific gas constant in J / kg.K
    GAMMA_AIR = 1.4
    EARTHRADIUS = 6356766.0

    GRAVITY_ACC = 9.80665       // accel. of GRAVITY_ACC
    MOL_WT  = 28.9644           // molecular weight of air
    R_GAS   = 8.31432           // gas constant

    // standard sealevel
    STD_SEALEVEL_TEMP = 288.15  // in Kelvin
    seaLevelDensity  = 1.2250 // kg/m^3
    seaLevelPressure = 101325.0 // in pa

)

type layer struct {
    name           string 
    base_altitude  float64 
    top_altitude   float64 
    baseTemp       float64 
    pressRatio     float64 
    lapseRate      float64
}

type isa76 struct {
    //var GMR     = 0.0      // hydrostatic constant
    layers []layer
}

func (isa *isa76) set() {
    seaLevelTemp := 0.0
    layersInitialized := false
    layersDisplay := ""
    customAltitude := ""
    lastAltitude := 0.0001
    lastUnit := ""

    GMR = GRAVITY_ACC * (MOL_WT / R_GAS) // hydrostatic constant
        
    seaLevelTemp = STD_SEALEVEL_TEMP
    
    isa.layers = append(isa.layers, layer{name: "Negative Troposphere", base_altitude:-5.004, top_altitude:0.0,    baseTemp: 320.676,  pressRatio: 1.7543,      lapseRate:-6.5 })
    isa.layers = append(isa.layers, layer{name: "Troposphere",        base_altitude:0.0,    top_altitude:11.0,    baseTemp: 288.15,  pressRatio: 1.0,          lapseRate:-6.5 })

    isa.layers = append(isa.layers, layer{name: "Tropopause",         base_altitude:11.0,   top_altitude:20.0,    baseTemp: 216.65,  pressRatio: 2.233611E-1,  lapseRate: 0.0 })
    isa.layers = append(isa.layers, layer{name: "Low Stratosphere", base_altitude:20.0,   top_altitude:32.0,    baseTemp: 216.65,  pressRatio: 5.403295E-2,  lapseRate: 1.0 })
    isa.layers = append(isa.layers, layer{name: "High Stratosphere", base_altitude:32.0,   top_altitude:47.0,    baseTemp: 228.65,  pressRatio: 8.5666784E-3, lapseRate: 2.8 })
    isa.layers = append(isa.layers, layer{name: "Stratopause",        base_altitude:47.0,   top_altitude:51.0,    baseTemp: 270.65,  pressRatio: 1.0945601E-3, lapseRate: 0.0 })
    isa.layers = append(isa.layers, layer{name: "Low Mesophere",    base_altitude:51.0,   top_altitude:71.0,    baseTemp: 270.65,  pressRatio: 6.6063531E-4, lapseRate:-2.8 })
    isa.layers = append(isa.layers, layer{name: "High Mesophere",    base_altitude:71.0,   top_altitude:84.852,  baseTemp: 214.65,  pressRatio: 3.9046834E-5, lapseRate:-2.0 })
    isa.layers = append(isa.layers, layer{name: "Mesopause",          base_altitude:84.852, top_altitude:86.0,    baseTemp: 186.946, pressRatio: 3.68501E-6,   lapseRate: 0.0 })

}

class isa76 {

    //let g0 = 9.80665    // (at sea level)
    let R_SGC = 287.05               // specific gas constant in J / kg.K
    let GAMMA_AIR = 1.4
    let EARTHRADIUS = 6356766.0

    let GRAVITY_ACC = 9.80665  // accel. of GRAVITY_ACC
    let MOL_WT  = 28.9644  // molecular weight of air
    let R_GAS   = 8.31432  // gas constant
    var GMR     = 0.0      // hydrostatic constant

    var layers = [(name: String, base_altitude:Double, top_altitude:Double, baseTemp:Double, pressRatio: Double, lapseRate: Double)]()

    // standard sealevel
    let STD_SEALEVEL_TEMP = 288.15  // in Kelvin
    let seaLevelDensity  = 1.2250 // kg/m^3
    let seaLevelPressure = 101325.0 // in pa
    
    var seaLevelTemp = 0.0
    
    //var path = NSURL()
    //var baseURL : NSURL = NSURL()
    


    var layersInitialized = false
    var layersDisplay = ""
    var customAltitude = ""
    var lastAltitude = 0.0001
    var lastUnit = ""

    
    let output = outputResult()
    
    var params : (  layerName: String,
                    geopotentialAltMeters: Double,
                    geopotentialAltFeet: Double,
                    geometricAltMeters: Double,
                    geometricAltFeet: Double,
                    temperatureK: Double,
                    temperatureC: Double,
                    temperatureF: Double,
                    temperatureR: Double,
                    pressure: Double,
                    density: Double,
                    soundSpeed:Double,
                    pressureRatio: Double,
                    tempRatio: Double,
                    densityRatio: Double) = (   layerName: "",          geopotentialAltMeters: 0.0, geopotentialAltFeet: 0.0,
                                                geometricAltMeters:0.0, geometricAltFeet: 0.0,      temperatureK: 0.0,
                                                temperatureC: 0.0,      temperatureF: 0.0,          temperatureR: 0.0,
                                                pressure:0.0,           density:0.0,                soundSpeed:0.0,
                                                pressureRatio: 0.0,     tempRatio: 0.0,             densityRatio: 0.0)

    

    init() {

        GMR = GRAVITY_ACC * (MOL_WT / R_GAS) // hydrostatic constant
        
        seaLevelTemp = STD_SEALEVEL_TEMP
        
        layers.append((name: "Negative Troposphere",base_altitude:-5.004, top_altitude:0.0,    baseTemp: 320.676,  pressRatio: 1.7543,      lapseRate:-6.5 ))
        layers.append((name: "Troposphere",        base_altitude:0.0,    top_altitude:11.0,    baseTemp: 288.15,  pressRatio: 1.0,          lapseRate:-6.5 ))

        layers.append((name: "Tropopause",         base_altitude:11.0,   top_altitude:20.0,    baseTemp: 216.65,  pressRatio: 2.233611E-1,  lapseRate: 0.0 ))
        layers.append((name: "Low Stratosphere", base_altitude:20.0,   top_altitude:32.0,    baseTemp: 216.65,  pressRatio: 5.403295E-2,  lapseRate: 1.0 ))
        layers.append((name: "High Stratosphere", base_altitude:32.0,   top_altitude:47.0,    baseTemp: 228.65,  pressRatio: 8.5666784E-3, lapseRate: 2.8 ))
        layers.append((name: "Stratopause",        base_altitude:47.0,   top_altitude:51.0,    baseTemp: 270.65,  pressRatio: 1.0945601E-3, lapseRate: 0.0 ))
        layers.append((name: "Low Mesophere",    base_altitude:51.0,   top_altitude:71.0,    baseTemp: 270.65,  pressRatio: 6.6063531E-4, lapseRate:-2.8 ))
        layers.append((name: "High Mesophere",    base_altitude:71.0,   top_altitude:84.852,  baseTemp: 214.65,  pressRatio: 3.9046834E-5, lapseRate:-2.0 ))
        layers.append((name: "Mesopause",          base_altitude:84.852, top_altitude:86.0,    baseTemp: 186.946, pressRatio: 3.68501E-6,   lapseRate: 0.0 ))
        
        // set path to local images...
        //path = NSBundle.mainBundle().URLForResource("Comparison_US_standard_atmosphere_1962", withExtension: "png")!
        //baseURL = NSURL(fileURLWithPath: path)
      
    }

    func setSealevelTemp(newTemp: String) {
        var temp = Double(newTemp)
        
        // temperature must be in kelvin
        if temp < 0.0 {
            temp = 0.0
        }
        self.seaLevelTemp = temp!
    }
    
    func setPressure(inout tuple:(name: String, base_altitude:Double, top_altitude:Double, baseTemp:Double, pressRatio: Double, lapseRate: Double)) {
        tuple.baseTemp = 0 //bs
    
    }
    
    func customAltitude(alt: String, unit: String) {

        var altitude = Double(alt)

        if unit == "f" {
            // convert altitude in metrics
            altitude = altitude! / 3.2808
        }

        // calculate geopotential altitude and perform sanity check
        var geopoAltitude = altitude! * (EARTHRADIUS/(EARTHRADIUS+altitude!))
        if geopoAltitude > 86000 {
            altitude = 86000
            geopoAltitude = altitude! * (EARTHRADIUS/(EARTHRADIUS+altitude!))
            
            //geopoAltitude = 86000 //-altitude!
            //altitude = EARTHRADIUS * (geopoAltitude / (EARTHRADIUS - geopoAltitude))
        }
        else if geopoAltitude < -5004 {
            altitude = -5000
            geopoAltitude = altitude! * (EARTHRADIUS/(EARTHRADIUS+altitude!))
            
            //geopoAltitude = -5004
            //altitude = EARTHRADIUS * (geopoAltitude / (EARTHRADIUS - geopoAltitude))
        }
        
        // convert in km
        altitude = altitude!/1000
        geopoAltitude = geopoAltitude / 1000
        
        if geopoAltitude != self.lastAltitude {
            self.lastAltitude = geopoAltitude
            getParameters(altitude!, geoPotentialAlt:geopoAltitude)
        }

    }
    
//func (isa *isa76) getParameters(geometricAltitude float64, geoPotentialAlt geopoAltitude:Double) {
func (isa *isa76) getParameters(geometricAltitude float64, geoPotentialAlt, geopoAltitude float64) {
    
    for _, layer := range isa.layers {
        
        if geopoAltitude >= layer.base_altitude && geopoAltitude <= layer.top_altitude {

            delta_h := geopoAltitude - layer.base_altitude
            localTemp := layer.baseTemp + (layer.lapseRate * delta_h)
            
            isa.params.tempRatio = localTemp / seaLevelTemp // -> aka theta
            isa.params.geopotentialAltMeters = geopoAltitude * 1000
            isa.params.geopotentialAltFeet = params.geopotentialAltMeters * 3.28084
            
            isa.params.geometricAltMeters = geometricAltitude * 1000
            isa.params.geometricAltFeet = isa.params.geometricAltMeters * 3.28084
            
            if layer.lapseRate == 0.0 {
                // isothermal
                isa.params.pressureRatio = layer.pressRatio * pow(M_E, -GMR*(delta_h/layer.baseTemp))   // -> aka delta
            }
            else {
                isa.params.pressureRatio = layer.pressRatio * pow(layer.baseTemp/localTemp, GMR/layer.lapseRate) // -> aka delta
            }
            
            isa.params.densityRatio = params.pressureRatio / params.tempRatio // -> aka sigma

            isa.params.layerName = layer.name
            isa.params.densityRatio = isa.params.pressureRatio / isa.params.tempRatio
            isa.params.temperatureK = isa.params.tempRatio * seaLevelTemp
            isa.params.temperatureC = isa.params.temperatureK - 273.15
            isa.params.temperatureF = ((isa.params.temperatureK - 273.15) * 1.8) + 32
            isa.params.temperatureR = (isa.params.temperatureK) * 1.8
            
            isa.params.pressure     = isa.params.pressureRatio * seaLevelPressure
            isa.params.density      = isa.params.densityRatio * seaLevelDensity
            isa.params.soundSpeed   = pow((GAMMA_AIR * isa.params.temperatureK * R_SGC), 0.5)
            
        }
    }
}
    
func showEarthAcceleration() string {
        
        //         Re^2         where Re is the radius of the earth
        // g = -------------    hg is the geometric altitude above sea level
        //     go (Re + hg)2
        
        glocal := GRAVITY_ACC * pow(EARTHRADIUS, 2) / pow(EARTHRADIUS + params.geometricAltMeters, 2)
        
        //postProcessed := output.normalizeNumber(String(format:"%2.6f",round(glocal*1000000)/1000000))
        //return postProcessed
        return fmt.Sprintf("%2.6f",round(glocal*1000000)/1000000))
}
 
/*
    private func formatTitle(title:String, abrvFrom: String) -> String {
    
        let additional = "<font size=-2>(Standard sea level)</font>"

        if abrvFrom == "f" {
            return output.startTable(additional, title: title+" at<br>"+output.normalizeNumber("\(round(params.geometricAltFeet))")+" ft ("+output.normalizeNumber("\(round(params.geometricAltMeters))")+" m)<br>")
        } else {
            return output.startTable(additional, title: title+" at<br>"+output.normalizeNumber("\(round(params.geometricAltMeters))")+" m ("+output.normalizeNumber("\(round(params.geometricAltFeet))")+" ft)<br>")
        }

    }
*/    
    func showTemperature(abrvFrom: String) -> String {
        var postProcessed = ""
        var localValue = ""
        localValue = formatTitle("Temperature", abrvFrom: abrvFrom)
        
        postProcessed = output.normalizeNumber(String(format:"%4.2f",round(self.params.temperatureK*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "K", comment: "Kelvin")
        
        postProcessed = output.normalizeNumber(String(format:"%4.2f",round(self.params.temperatureC*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "C", comment: "Celcius")
        
        postProcessed = output.normalizeNumber(String(format:"%4.2f",round(self.params.temperatureF*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "F", comment: "Fahrenheit")
        
        postProcessed = output.normalizeNumber(String(format:"%4.2f",round(self.params.temperatureR*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "R", comment: "Rankine")
        
        postProcessed = output.normalizeNumber(String(format:"%4.2f",round((self.params.tempRatio)*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "", comment: "[θ] (Ratio t/t0)")
        localValue += output.endTable()
        
        return localValue
        
    }
  
    
    func showPressure(abrvFrom: String) -> String {
        var postProcessed = ""
        var localValue = ""
        
        localValue = formatTitle("Pressure", abrvFrom: abrvFrom)

        
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round(self.params.pressure*100000)/100000))
        localValue += output.addLineToTable(postProcessed, unit: "pa", comment: "pascals")
   
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressure * 1E-03)*1E+9)/1E+9))
        localValue += output.addLineToTable(postProcessed, unit: "kpa", comment: "kilopascals")

        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressure * 1E-05)*1E+9)/1E+9))
        localValue += output.addLineToTable(postProcessed, unit: "bar", comment: "bars")
        
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressure * 1E-02)*1E+6)/1E+6))
        localValue += output.addLineToTable(postProcessed, unit: "mbar", comment: "millibars")
        
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressure/6.8947572798677E+03)*1E+5)/1E+5))
        localValue += output.addLineToTable(postProcessed, unit: "lb/in²", comment: "pound/inch²")
        
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressure/0.47880258888E+02)*1E+5)/1E+5))
        localValue += output.addLineToTable(postProcessed, unit: "lb/ft²", comment: "pound/foot²")
        
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressure/1.013250E+05)*1E+9)/1E+9))
        localValue += output.addLineToTable(postProcessed, unit: "atm", comment: "atmosphere")
        
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressure/1.3332239E+02)*1E+5)/1E+5))
        localValue += output.addLineToTable(postProcessed, unit: "mmhg", comment: "mm of mercury")
        
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round((self.params.pressureRatio)*1E+8)/1E+8))
        localValue += output.addLineToTable(postProcessed, unit: "", comment: "[δ] (Ratio p/p0)")
        
        localValue += output.endTable()
        
        return localValue
    }
    

    func showSoundSpeed(abrvFrom: String) -> String {
        var postProcessed = ""
        var localValue = ""
        
        localValue = formatTitle("Speed of Sound", abrvFrom: abrvFrom)

        
        let metricValue = pow((GAMMA_AIR * self.params.temperatureK * R_SGC), 0.5)
        
        postProcessed = output.normalizeNumber(String(format:"%4.3f",round(metricValue*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "m/s", comment: "meters/second")

        postProcessed = output.normalizeNumber(String(format:"%4.3f",round((metricValue/340.37713655487805)*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "mach", comment: "mach number")
        
        postProcessed = output.normalizeNumber(String(format:"%4.3f",round((metricValue*3.281)*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "ft/s", comment: "feet/second")
        
        postProcessed = output.normalizeNumber(String(format:"%4.3f",round((metricValue*3.6)*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "km/h", comment: "kilometers/hour")

        postProcessed = output.normalizeNumber(String(format:"%4.3f",round((metricValue*(0.000621371*3600))*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "mi/h", comment: "miles/hour")
        
        postProcessed = output.normalizeNumber(String(format:"%4.3f",round((metricValue*1.94384)*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "kts", comment: "knots")
        localValue += output.endTable()

        return localValue
    }
    
    
    func showDensity(abrvFrom: String) -> String {
        var postProcessed = ""
        var localValue = ""
        
        localValue = formatTitle("Density", abrvFrom: abrvFrom)

        
        postProcessed = output.normalizeNumber(String(format:"%3.9f",round(self.params.density*1000000000)/1000000000))
        localValue += output.addLineToTable(postProcessed, unit: "kg/m³", comment: "kilograms/meter³")
        
        postProcessed = output.normalizeNumber(String(format:"%3.9f",round((self.params.density * 1.94024058983E-03)*1000000000)/1000000000))
        localValue += output.addLineToTable(postProcessed, unit: "slug/ft³", comment: "slug/foot³")
        
        postProcessed = output.normalizeNumber(String(format:"%3.9f",round((self.params.density/1000)*1000000000)/1000000000))
        localValue += output.addLineToTable(postProcessed, unit: "g/cm³", comment: "gram/centimeter³")
        
        postProcessed = output.normalizeNumber(String(format:"%3.9f",round((self.params.density * 1.94024058983E-03/53.7056928034)*1000000000)/1000000000))
        localValue += output.addLineToTable(postProcessed, unit: "lb/in³", comment: "pound/inch³")

        postProcessed = output.normalizeNumber(String(format:"%3.9f",round((self.params.density * 1.94024058983E-03/3.10848616723E-02)*1000000000)/1000000000))
        localValue += output.addLineToTable(postProcessed, unit: "lb/ft³", comment: "pound/foot³")
        
        postProcessed = output.normalizeNumber(String(format:"%3.9f",round((self.params.densityRatio)*1000000000)/1000000000))
        localValue += output.addLineToTable(postProcessed, unit: "", comment: "[σ] (Ratio ρ/ρ0)")
        localValue += output.endTable()
 
        return localValue
    }
    
    func showAllMetrics(abrvFrom: String) -> String {
        var postProcessed = ""
        var localValue = ""
        
        localValue = formatTitle("All Parameters", abrvFrom: abrvFrom)

        // layer name
        localValue += output.addHeadlineToTable(params.layerName, unit: "", comment: "Layer Name")
        
        // geometric
        postProcessed = output.normalizeNumber(String(format:"%5.2f",round(self.params.geometricAltMeters*100)/100))
        localValue += output.addLineToTable(postProcessed, unit: "m", comment: "Geometric Alt.")
        
        // geopotential
        postProcessed = output.normalizeNumber(String(format:"%5.2f",round(self.params.geopotentialAltMeters*100)/100))
        localValue += output.addLineToTable(postProcessed, unit: "m", comment: "Geopotential Alt.")
        
        // density
        postProcessed = output.normalizeNumber(String(format:"%3.9f",round(self.params.density*1000000000)/1000000000))
        localValue += output.addLineToTable(postProcessed, unit: "kg/m³", comment: "Density")
        
        // speed of sound
        let metricValue = pow((GAMMA_AIR * self.params.temperatureK * R_SGC), 0.5)
        postProcessed = output.normalizeNumber(String(format:"%4.3f",round(metricValue*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "m/s", comment: "Speed of sound")
        
        // pressure
        postProcessed = output.normalizeNumber(String(format:"%7.6f",round(self.params.pressure*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "pa", comment: "Pressure")

        // temperature
        postProcessed = output.normalizeNumber(String(format:"%4.2f",round(self.params.temperatureK*1000)/1000))
        localValue += output.addLineToTable(postProcessed, unit: "K", comment: "Temperature")
        
        // gravity acceleration
        // temperature
        localValue += output.addLineToTable(showEarthAcceleration(), unit: "m/s²", comment: "Accel. of gravity")

        localValue += output.endTable()
        return localValue
   
    }
    
    func showCustomAltitude(fromValue: String) -> String {
        var postProcessed = ""
        if lastAltitude > -610 {
            
            customAltitude = "Geometric Alt:".uppercaseString+"\n"
            postProcessed = output.normalizeNumber("\(round(self.params.geometricAltMeters*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" m ("
            postProcessed = output.normalizeNumber("\(round(self.params.geometricAltFeet*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" ft)\n"
            
            customAltitude = customAltitude+"\nGeopotential Alt:".uppercaseString+"\n"
            postProcessed  = output.normalizeNumber("\(round(self.params.geopotentialAltMeters*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" m ("
            postProcessed  = output.normalizeNumber("\(round(self.params.geopotentialAltFeet*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" ft)\n"
            
            customAltitude = customAltitude+"\nTemperature:".uppercaseString+"\n"
            postProcessed = output.normalizeNumber("\(round(self.params.temperatureK*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" Kelvin\n"
            
            postProcessed = output.normalizeNumber("\(round(self.params.temperatureC*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" Celcius\n"
            postProcessed = output.normalizeNumber("\(round(self.params.temperatureF*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" Fahrenheit\n"
            
            
            
            customAltitude = customAltitude+"\nPressure:".uppercaseString+"\n"
            postProcessed = output.normalizeNumber("\(round(self.params.pressure*10000000)/10000000)")
            customAltitude = customAltitude+postProcessed+" pa\n"
            
            customAltitude = customAltitude+"\nDensity:".uppercaseString+"\n"
            postProcessed = output.normalizeNumber("\(round(self.params.density*100000000)/100000000)")
            customAltitude = customAltitude+postProcessed+" Kg/m³\n"
            
            customAltitude = customAltitude+"\nSpeed of sound:".uppercaseString+"\n"
            postProcessed = output.normalizeNumber("\(round(pow((GAMMA_AIR * self.params.temperatureK * R_SGC), 0.5)*1000)/1000)")
            customAltitude = customAltitude+postProcessed+" m/s\n"
            
            customAltitude = customAltitude+"\nLayer Name:".uppercaseString+"\n\(self.params.layerName)\n"
            
        }
        return customAltitude
    }
    
    func showAbout() -> String {
        //let path = NSBundle.mainBundle().pathForResource("Comparison_US_standard_atmosphere_1962", ofType: "svg")
        //let baseURL = NSURL(fileURLWithPath: path!)
        
        //let myUrl = NSURL(string: "Comparison_US_standard_atmosphere_1962.svg");
        let help = "<div style='text-align:justify;width:100%;'>The U.S. Standard Atmosphere is an atmospheric model of how the pressure, temperature, density, and viscosity of the Earth's atmosphere change over a wide range of altitudes or elevations. The model, based on an existing international standard, was first published in 1958 by the U.S. Committee on Extension to the Standard Atmosphere, and was updated to its most recent version in 1976. It is largely consistent in methodology with the International Standard Atmosphere, differing mainly in the assumed temperature distribution at higher altitudes.<br><br>This USSA calculator will determine the temperature, pressure, density, speed of sound, geopotential altitude, acceleration of gravity for any altitude between -5,000m (-16,404ft) and +86,000m (+282,152ft) given a standard sea-level temperature of 288.15 Kelvin and pressure of 101,325 pascals.<br><br><img src='USstandardAtmosphere1976.jpg' style='width:100%;'></div>"
        
        return help
        
    }
    
    func showAllLayers() -> String {
        return "whoops!"
        /*
        if !layersInitialized {
            var postProcessed = ""
            for layer in atmLayers {
                layersDisplay = layersDisplay+"\(layer.name)".uppercaseString+" layer:\n"
                postProcessed = normalizeNumber("\(Int(layer.base_altitude))")
                
                layersDisplay = layersDisplay+"["+postProcessed
                postProcessed = normalizeNumber("\(Int(layer.top_altitude))")
                
                layersDisplay = layersDisplay+" m ➡︎ "+postProcessed+" m]\n"
                layersDisplay = layersDisplay+"Lapse Rate    = \(layer.lapse_rate) K/m\n"
                postProcessed = normalizeNumber("\(Int(layer.base_altitude))")
                layersDisplay = layersDisplay+"Base Altitude = "+postProcessed+" m\n"
                postProcessed = normalizeNumber("\(Int(layer.top_altitude))")
                layersDisplay = layersDisplay+"Top  Altitude = "+postProcessed+" m\n"
                layersDisplay = layersDisplay+"Base Temp (K) = \(round(layer.base_temperature_K*1000)/1000) K\n"
                layersDisplay = layersDisplay+"Base Temp (C) = \(round(layer.base_temperature_C*1000)/1000) C\n"
                layersDisplay = layersDisplay+"Base Temp (F) = \(round(layer.base_temperature_F*1000)/1000) F\n"
                print("before-> \(round(layer.base_pressure*1000)/1000)")
                postProcessed = normalizeNumber("\(round(layer.base_pressure*1000)/1000)")
                print("after-> "+postProcessed)
                layersDisplay = layersDisplay+"Base Pressure = "+postProcessed+" pa\n"
                layersDisplay = layersDisplay+"Base Density  = \(round(layer.base_density*1000000)/1000000) kg/m³\n"
                layersDisplay = layersDisplay+"Top  Density  = \(round(layer.top_density*1000000)/1000000) kg/m³\n\n"
                
            }
            layersInitialized = true
        }
        return layersDisplay*/
    }
    
}