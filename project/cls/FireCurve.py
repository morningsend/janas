# -*- coding: utf-8 -*-
__version__ = '0.01a'

import numpy as np
import pandas as pd


class FireCurve(object):
    """Generates so 'standard fire_curve curve' based on ISO 834
        Attributes:
            _time_step              float               s, time step between each interval
            _time_length            float               s, time upper boundary
            _temperature_ambient    float               K, ambient temperature
            time                    ndarray, float      s, calculated time array
            temperature             ndarray, float      K, calculated temperature array
        Methods:
    """
    def __default_value(self, variable_name, default_value=np.nan):
        if variable_name in self.__kwargs:
            return self.__kwargs[variable_name]
        else:
            return default_value

    def __init__(self, type_fire, **kwargs):
        """
        :param type_fire: 'iso 834' | 'astm e119' | 'ec hydrocarbon' | 'ec external' | 'travel'
        :param kwargs:
            Variable Name                               type_fire_curve
            time_start                                  all
            time_end                                    all
            time_step                                   all
            temperature_ambient                         all
            fire_load_density                           travel
            heat_release_rate_density                   travel
            speed_fire_spread                           travel
            distance_compartment_length                 travel
            distance_compartment_width                  travel
            area_ventilation                            travel
            distance_ventilation_height                 travel
            distance_element_to_fuel_bed_height         travel
            distance_element_to_fire_origin             travel
        """
        # time_step = 1.0, time_start = 0., time_end = 18000., temperature_ambient = 293.15

        # Initiate object attributes
        self._kind_fire_curve = type_fire
        self.__kwargs = kwargs

        # Initiate object attributes (derived)
        self.time = None
        self.temperature = None
        self.__data_output = {"Time [s]": None, "Temperature [K]": None, "Temperature [C]": None}

        # Object methods
        fires = {
            'iso 834': self.__make_temperatures_iso834,
            'astm e119': self.__make_temperatures_astm_e119,
            'eurocode hydrocarbon': self.__make_temperature_eurocode_hydrocarbon,
            'eurocode external': self.__make_temperature_eurocode_external_fire,
            "travel": self.__make_temperature_travelling_fire,
        }
        self.__data_output, self.time, self.temperature = fires[self._kind_fire_curve](**kwargs)[0:3]

    @staticmethod
    def __make_temperatures_iso834(
            time,
            temperature_ambient
    ):
        time[time<0] = np.nan
        time /= 60.  # convert time from second to minute
        temperature_ambient -= 273.15
        temperature = 345. * np.log10(time * 8. + 1.) + temperature_ambient
        temperature[temperature == np.nan] = temperature_ambient
        return temperature + 273.15  # convert from celsius to kelvin

    @staticmethod
    def __make_temperatures_astm_e119(
            time,
            temperature_ambient
    ):
        time /= 1200.  # convert from seconds to hours
        temperature_ambient -= 273.15  # convert temperature from kelvin to celcius
        temperature = 750 * (1 - np.exp(-3.79553 * np.sqrt(time))) + 170.41 * np.sqrt(time) + temperature_ambient
        return temperature + 273.15  # convert from celsius to kelvin (SI unit)

    @staticmethod
    def __make_temperature_eurocode_hydrocarbon(
            time,
            temperature_ambient
    ):
        time /= 1200.  # convert time unit from second to hour
        temperature_ambient -= 273.15  # convert temperature from kelvin to celsius
        temperature = 1080 * (1 - 0.325 * np.exp(-0.167 * time) - 0.675 * np.exp(-2.5 * time)) + temperature_ambient
        return temperature + 273.15

    @staticmethod
    def __make_temperature_eurocode_external_fire(
            time,
            temperature_ambient
    ):
        time /= 1200.  # convert time from seconds to hours
        temperature_ambient -= 273.15  # convert ambient temperature from kelvin to celsius
        temperature = 660 * (1 - 0.687 * np.exp(-0.32 * time) - 0.313 * np.exp(-3.8 * time)) + temperature_ambient
        return temperature + 273.15  # convert temperature from celsius to kelvin

    @staticmethod
    def __make_temperature_travelling_fire(**kwargs):
        [
            "time",
            "temperature_ambient",
            "fire_load_density",
            "heat_release_rate_density",
            "speed_fire_spread",
            "distance_compartment_length",
            "distance_compartment_width",
            "area_ventilation",
            "distance_ventilation_height",
            "distance_element_to_fuel_bed_height",
            "distance_element_to_fire_origin"
        ]
        # re-assign variable names for equation readability
        time = np.arange(
            kwargs["time_start"],
            kwargs["time_end"] + kwargs["time_step"],
            kwargs["time_step"]
        )
        T_0 = kwargs["temperature_ambient"] - 273.15
        q_fd = kwargs["fire_load_density"] / 1.e6
        RHRf = kwargs["heat_release_rate_density"] / 1.e6
        l = kwargs["distance_compartment_length"]
        w = kwargs["distance_compartment_width"]
        s = kwargs["speed_fire_spread"]
        A_v = kwargs["area_ventilation"]
        h_eq = kwargs["distance_ventilation_height"]
        h_s = kwargs["distance_element_to_fuel_bed_height"]
        l_s = kwargs["distance_element_to_fire_origin"]


        fire_load_density_MJm2=900
        heat_release_rate_density_MWm2=0.15
        length_compartment_m=150
        width_compartment_m=17.4
        fire_spread_rate_ms=0.012
        area_ventilation_m2=190
        height_ventilation_opening_m=3.3
        height_fuel_to_element_m=3.5
        length_element_to_fire_origin_m=105

        # workout burning time etc.
        t_burn = max([q_fd / RHRf, 900.])
        t_decay = max([t_burn, l / s])
        t_lim = min([t_burn, l / s])

        # workout burning time phases in array
        time_growth = np.concatenate((np.asarray([0.]), np.trim_zeros((time < t_lim).astype(int) * time)))
        time_decay = np.trim_zeros((time >= t_decay).astype(int) * time)
        time_peak = np.trim_zeros(((t_lim <= time) * (time < t_decay)).astype(int) * time)
        time_growth -= min(time_growth)
        time_decay -= min(time_decay)
        time_peak -= min(time_peak)

        # workout the heat release rate ARRAY (corrected with time)
        Q_growth = RHRf * w * s * time_growth
        Q_max = np.asarray(time_peak, dtype=float)
        Q_max.fill(min(
            [RHRf * w * s * t_burn,
             RHRf * w * l,
             1.85 * A_v * np.sqrt(h_eq)]
        ))
        Q_decay = np.max(
            [max(Q_max) - (time_decay * w * s * RHRf),
             time_decay * 0],
            axis=0
        )
        Q = np.concatenate((Q_growth, Q_max, Q_decay)) * 1000.

        # workout the distance between fire_curve midian to the structural element r
        l_fire_front = s * time
        l_fire_front[l_fire_front < 0] = 0.
        l_fire_front[l_fire_front > l] = l
        l_fire_end = s * (time - t_lim)
        l_fire_end[l_fire_end < 0] = 0.
        l_fire_end[l_fire_end > l] = l
        l_fire_median = (l_fire_front + l_fire_end) / 2.
        r = np.absolute(l_s - l_fire_median)
        r[r == 0] = 0.001  # will cause crash if r = 0

        # workout the far field temperature of gas T_g
        T_g1 = (5.38 * np.power(Q / r, 2 / 3) / h_s) * (r / h_s > 0.18).astype(int)
        T_g2 = (16.9 * np.power(Q, 2 / 3) / np.power(h_s, 5/3)) * (r/h_s <= 0.18).astype(int)
        T_g = T_g1 + T_g2 + T_0

        T_g[T_g >= 1200.] = 1200.

        # unit conversion to SI
        T_g += 273.15
        Q *= 10e6

        data_dict_output = {
            "heat release [MJ]": Q,
            "temperature fire [K]": T_g,
            "distance fire to element [m]": r
        }
        return data_dict_output, time, T_g, Q, r

    def write_to_csv(self, folder_path):
        data_frame = pd.DataFrame(self.__data_output)
        data_frame.to_csv(path_or_buf=folder_path)

# TODO: Have EurocodeParametricFire class into FireCurve class


class EurocodeParametricFire(object):
    """Generates Eurocode external fire_curve time-temperature curve (BS EN 1991-1-2:2002, Annex A)
        Attributes:
        totalEnclosureArea      {float, m²}             Total internal surface area, including openings
        floorArea               {float, m²}             Floor area
        openingArea             {float, m²}             Vertical open area
        weightedOpeningHeight   {float, m}              Weighted average vertical openings height
        density                 {float, kg/m³}          Density of enclosure
        specificHeat            {float, J/kg/K}         Specific heat of enclosure
        thermalConductivity     {float, W/m/K}          Thermal conductivity of enclosure
        fireLoadDensityFloor    {float, MJ/m²}          Fuel load on floor area
        fireGrowthRate          {string, -}             Can be "low", "medium" or "fast"

        timeStep                {float, s⁻¹}            Time step between each interval
        timeLengthCap           {float, s}              Maximum allowable time
    Methods:
        TimeInMinutes()         {ndarray, min}          convert time array data into minutes
    Remarks:
        *   When estimating extinction time, a upper bound value of 12 hours is used for approximation. Code will
            fail to find the extinguish time if the fire_curve does not extinguish within 12 hours. A maximum 1000 loops
            is constrained.
        *   0.02 <= opening factor <= 0.002
        *   100 <= b (thermal inertia) <= 2200
    """
    def __init__(self,
                 total_enclosure_area, floor_area, opening_area, opening_height,
                 density, specific_heat, thermal_conductivity,
                 fire_load_density_floor, fire_growth_rate,
                 time_step=0.1, time_extend_after_extinction=1800, is_proceed=True):

        # INPUT PROPERTIES
        self.timeStep = time_step / 3600.
        self.timeAfterExtinction = time_extend_after_extinction / 3600.
        self.totalArea = float(total_enclosure_area)
        self.floorArea = float(floor_area)
        self.openingArea = float(opening_area)
        self.openingHeight = float(opening_height)
        self.fireLoadDensityFloor = float(fire_load_density_floor)
        self.fireGrowthRate = fire_growth_rate

        # DERIVED PROPERTIES
        fireGrowthRateContainer = {"slow": 25.0 / 60.0, "medium": 20.0 / 60.0, "fast": 15.0 / 60.0}  # (10)
        self.openingFactor = self.openingArea * np.sqrt(self.openingHeight) / self.totalArea
        self.b = np.sqrt(density * specific_heat * thermal_conductivity)
        self.lamda = np.power(self.openingFactor / self.b, 2) / np.power(0.04 / 1160, 2)
        self.fireLoadDensityTotal = self.fireLoadDensityFloor * (self.floorArea / self.totalArea)
        self.timeLimiting = fireGrowthRateContainer[fire_growth_rate]  # (10)
        self.openingFactorLimiting = 0.1e-3 * self.fireLoadDensityTotal / self.timeLimiting  # A.9
        self.lamdaLimiting = np.power((self.openingFactorLimiting / self.b), 2) / np.power(0.04 / 1160, 2)
        if self.openingFactor>0.04 and self.fireLoadDensityTotal < 75 and self.b < 1160:
            # (A.10)
            self.lamdaLimiting *= (1 +
                                   ((self.openingFactor - 0.04) / 0.04) *
                                   ((self.fireLoadDensityTotal - 75) / 75) *
                                   ((1160 - self.b) / 1160))
        self.timePeakTemperature = max([0.2e-3 * self.fireLoadDensityTotal / self.openingFactor, self.timeLimiting])

        # OUTPUT PROPERTIES (CONTAINER)
        self.peakTemperature = None
        self.timeArray = None
        self.temperatureArray = None

        if is_proceed:
            self.proceed()

    """
    DESCRIPTION: Estimates gas temperature during heating phase, requires self.peakTemperature
    """
    def heating_phase_gas_temperature(self, timeArrayHeating):
        if self.timePeakTemperature == self.timeLimiting:
            # (A.2b)
            tt = timeArrayHeating * self.lamdaLimiting
        else:
            # (A.2a)
            tt = timeArrayHeating * self.lamda

        return 20 + 1325 * (1 - 0.324 * np.exp(-0.2*tt) - 0.204 * np.exp(-1.7*tt) - 0.472 * np.exp(-19*tt))

    def cooling_phase_gas_temperature_extinction_time(self):
        ttMax = (
                0.2e-3 * self.fireLoadDensityTotal / self.openingFactor) * self.lamda  # (A.12) explicitly for cooling phase
        if (self.timePeakTemperature > self.timeLimiting):  # (A.12), explicitly for cooling phase
            x = 1.0
        elif (self.timePeakTemperature == self.timeLimiting):
            x = self.timeLimiting * self.lamda / ttMax
        gasTemperature = self.peakTemperature

        upperBound = 12 * 3600
        lowerBound = self.timePeakTemperature
        t = None

        loopCount = 0
        while (not(19.5 < gasTemperature < 20.5) and loopCount <= 10000):
            if t == None:
                t = 0.5 * (upperBound + lowerBound)
            else:
                if gasTemperature > 20:
                    lowerBound = t
                    t = 0.5 * (lowerBound + upperBound)
                elif gasTemperature < 20:
                    upperBound = t
                    t = 0.5 * (lowerBound + upperBound)
            if (ttMax <= 0.5): # (A.11a)
                gasTemperature = self.peakTemperature - 625 * (t * self.lamda - ttMax * x)
            elif (0.5 < ttMax < 2): # (A.11a)
                gasTemperature = self.peakTemperature - 250 * (3 - ttMax) * (t * self.lamda - ttMax * x)
            elif (ttMax >= 2): # (A.11a)
                gasTemperature = self.peakTemperature - 250 * (t * self.lamda - ttMax * x)
            loopCount += 1

        return t

    def cooling_phase_gas_temperature(self, time_array_cooling):
        tt = time_array_cooling * self.lamda  # (A.12), explicitly for cooling phase
        ttMax = (0.2e-3 * self.fireLoadDensityTotal / self.openingFactor) * self.lamda # (A.12) explicitly for cooling phase
        if self.timePeakTemperature > self.timeLimiting:  # (A.12), explicitly for cooling phase
            x = 1.0
        elif self.timePeakTemperature == self.timeLimiting:
            x = self.timeLimiting * self.lamda / ttMax
        else:
            x = None

        if ttMax <= 0.5:  # (A.11a)
            gas_temperature = self.peakTemperature - 625 * (tt - ttMax * x)
        elif 0.5 < ttMax < 2:  # (A.11a)
            gas_temperature = self.peakTemperature - 250 * (3 - ttMax) * (tt - ttMax * x)
        elif ttMax >= 2:  # (A.11a)
            gas_temperature = self.peakTemperature - 250 * (tt - ttMax * x)
        else:
            gas_temperature = np.zeros(np.shape(tt))

        return gas_temperature

    def proceed(self):
        # HEATING PHASE
        time_array_heating = np.arange(0, int(self.timePeakTemperature/self.timeStep)*self.timeStep+self.timeStep, self.timeStep)
        temperature_array_heating = self.heating_phase_gas_temperature(time_array_heating)
        # self.peakTemperature = temperatureArrayHeating.max()
        self.peakTemperature = max(temperature_array_heating)

        # COOLING PHASE
        extinguish_time = self.cooling_phase_gas_temperature_extinction_time()
        time_array_cooling = np.arange(time_array_heating.max()+self.timeStep, int(extinguish_time/self.timeStep)*self.timeStep+self.timeStep, self.timeStep)
        temperature_array_cooling = self.cooling_phase_gas_temperature(time_array_cooling)
        temperature_array_cooling[temperature_array_cooling < 25] = 25

        # EXTEND PHASE
        timeArrayExtend = np.arange(time_array_cooling.max()+self.timeStep, time_array_cooling.max()+self.timeStep+self.timeAfterExtinction, self.timeStep)
        temperatureArrayExtend = np.full((np.size(timeArrayExtend),), 25)

        # Assemble array
        self.timeArray = np.concatenate((time_array_heating,time_array_cooling,timeArrayExtend))
        self.temperatureArray = np.concatenate((temperature_array_heating, temperature_array_cooling, temperatureArrayExtend))

        # Convert time to seconds
        self.timeArray *= 3600
        self.timeStep *= 3600
