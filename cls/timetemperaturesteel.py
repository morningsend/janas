# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate


def steel_property_(property_type='c_s, k, reduction_factor', steel_type='carbon steel'):
    def steel_specific_heat_carbon_steel(temperature):
        """
            DESCRIPTION:
                [BS EN 1993-1-2:2005, 3.4.1.2]
                Calculate steel specific heat according to temperature
            PARAMETERS:
                temperature     {float, K}          Given temperature

                __return__      {float, J/kg/K}     Specific heat
        """
        temperature -= 273.15
        if 20 <= temperature < 600:
            return 425 + 0.773 * temperature - 1.69e-3 * np.power(temperature, 2) + 2.22e-6 * np.power(temperature, 3)
        elif 600 <= temperature < 735:
            return 666 + 13002 / (738 - temperature)
        elif 735 <= temperature < 900:
            return 545 + 17820 / (temperature - 731)
        elif 900 <= temperature <= 1200:
            return 650
        else:
            return 0

    def thermal_conductivity_carbon_steel(temperature):
        """
            DESCRIPTION:
                [BS EN 1993-1-2:2005, 3.4.1.3]
            PARAMETERS:
            OUTPUTS:
            REMARKS:
        """
        temperature += 273.15
        if 20 <= temperature <= 800:
            return 54 - 0.0333 * temperature
        elif 800 <= temperature <= 1200:
            return 27.3
        else:
            return 0

    def reduction_factor_carbon_steel(steel_temperature):
        """
        DESCRIPTION:
            [BS EN 1993-1-2:2005, Table 3.1]
            Return reductions factors given temperature.
        PARAMETERS:
            temperature     {float, K}  Steel temperature

            self            {tuple, -}  (k_y_theta, k_p_theta, k_E_theta)
            k_y_theta       {float, /}  reduction factor for effective yield strength
            k_p_theta       {float, /}  reduction factor for proportional limit
            k_E_theta       {float, /}  reduction factor for the slope of the linear elastic range
        REMARKS:
            1.  293 [K] < steel_temperature < 1473 [K]
        """
        steel_temperature -= 273.15  # convert from [K] to [C], as the data below is in [C]
        temperature = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
        k_y_theta = [1, 1, 1, 1, 1, 0.78, 0.47, 0.23, 0.11, 0.06, 0.04, 0.02, 0]
        k_p_theta = [1, 1, 0.807, 0.613, 0.42, 0.36, 0.18, 0.075, 0.05, 0.0375, 0.025, 0.0125, 0]
        k_E_theta = [1, 1, 0.9, 0.8, 0.7, 0.6, 0.31, 0.13, 0.09, 0.0675, 0.045, 0.0225, 0]

        k_y_theta = np.interp(steel_temperature, temperature, k_y_theta)
        k_p_theta = np.interp(steel_temperature, temperature, k_p_theta)
        k_E_theta = np.interp(steel_temperature, temperature, k_E_theta)

        return k_y_theta, k_p_theta, k_E_theta

    function_dictionary = {
        "specific heat": steel_specific_heat_carbon_steel,
        "thermal conductivity": thermal_conductivity_carbon_steel,
        "reduction factor": reduction_factor_carbon_steel
    }

    return function_dictionary[property_type]


class SteelTemperature(object):
    """
    DESCRIPTION:
    ATTRIBUTES:
            timeArray                           {ndarray, s}
            temperatureArray                    {ndarray, K}
            sectionPerimeter                    {m}
            sectionArea                         {m2}
            densitySteel                        {kg/m3}K
            convectiveHeatTransferCoefficient   {W/m2/K}
            resultantEmissivity                 {-}
    REMARKS:
            1. Steel temperature is limited within 20C <= T_s <= 1200, due to specific heat
    """
    def __init__(self, fire_curve, type_steel, **kwargs):
        # Attributes, compulsory
        self._fire_curve = fire_curve           # FireCurve object
        self._type_steel = type_steel
        self.__kwargs = kwargs                  # keyword arguments, varies for different 'fire_curve'

        # Attributes, ouput container
        self._data_output = None
        self.temperature = None

        self.time = self._fire_curve.time
        self.temperature_fire = self._fire_curve.temperature

        self.__check_parameters_for_type()

        self.__make_temperature_eurocode_protected_steel(
            time=self.__kwargs["time"],
            temperature_fire=self.__kwargs["temperature_fire"],
            density_steel=self.__kwargs["density_steel"],
            c_steel_T=self.__kwargs["c_steel_T"],
            area_steel_section=self.__kwargs["area_steel_section"],
            k_protection=self.__kwargs["k_protection"],
            density_protection=self.__kwargs["density_protection"],
            c_protection=self.__kwargs["c_protection"],
            thickness_protection=self.__kwargs["thickness_protection"],
            perimeter_protected=self.__kwargs["perimeter_protected"]
        )

    # def __default_value(self, variable_name, default_value=np.nan):
    #     if variable_name in self.__kwargs:
    #         return self.__kwargs[variable_name]
    #     else:
    #         return default_value

    def __check_parameters_for_type(self):
        # Set required parameters for different options/functions
        params_dict = {
            "eurocode unprotected":
                [
                    "time",
                    "temperature_fire",
                    "perimeter_section",
                    "area_section",
                    "perimeter_box",
                    "density_steel",
                    "h_conv",
                    "emissivity_resultant",
                ],
            "eurocode protected":
                [
                    "time",
                    "temperature_fire",
                    "density_steel",
                    "c_steel_T",
                    "area_steel_section",
                    "k_protection",
                    "density_protection",
                    "c_protection",
                    "thickness_protection",
                    "perimeter_protected",
                ],
        }

        # Error handler when parameters are missing
        if self._type_steel not in params_dict:
            raise ValueError("Unexpected steel type: {0}".format(self._type_steel))
        else:
            params_missing, params_required = [], params_dict[self._type_steel]

            for each_param in params_required:
                if each_param not in self.__kwargs:
                    params_missing.append(each_param)

            if len(params_missing) > 0:
                raise ValueError("Missing required parameter: {0}.".format(", ".join(params_missing)))

    @staticmethod
    def __make_temperature_eurocode_unprotected_steel(
            time,
            temperature_fire,
            perimeter_section,
            area_section,
            perimeter_box,
            density_steel,
            h_conv,
            emissivity_resultant
    ):

        # Create steel temperature change array s
        temperature_rate_steel = time * 0.
        temperature_steel = time * 0.
        heat_flux_net = time * 0.
        c_s = time * 0.

        k_sh = 0.9 * (perimeter_box / area_section) / (perimeter_section / area_section)  # BS EN 1993-1-2:2005 (e4.26a)
        F = perimeter_section
        V = area_section
        rho_s = density_steel
        h_c = h_conv
        sigma = 56.7e-9  # todo: what is this?
        epsilon = emissivity_resultant
        c_s_ = steel_property_('specific heat')

        time_, temperature_steel[0], c_s[0] = iter(time), temperature_fire[0], 0.
        next(time_)
        for i, v in enumerate(time_):
            i += 1
            T_f, T_s_ = temperature_fire[i], temperature_steel[i-1]  # todo: steel specific heat
            c_s[i] = c_s_(temperature_steel[i-1])

            # BS EN 1993-1-2:2005 (e4.25)
            a = h_c * (T_f - T_s_)
            b = sigma * epsilon * (np.power(T_f,4) - np.power(T_s_,4))
            c = k_sh * F / V / rho_s / c_s[i]
            d = time[i] - time[i-1]

            heat_flux_net[i] = a + b

            temperature_rate_steel[i] = c * (a + b) * d
            temperature_steel[i] = temperature_steel[i-1] + temperature_rate_steel[i]

        return temperature_steel, temperature_rate_steel, heat_flux_net, c_s

    @staticmethod
    def __make_temperature_eurocode_protected_steel(
            time,
            temperature_fire,
            density_steel,
            c_steel_T,
            area_steel_section,
            k_protection,
            density_protection,
            c_protection,
            thickness_protection,
            perimeter_protected
    ):
        params_required = [
            "time",
            "temperature_fire",
            "density_steel",
            "c_steel_T",
            "area_steel_section",
            "k_protection",
            "density_protection",
            "c_protection",
            "thickness_protection",
            "perimeter_protected",
        ]
        
        V = area_steel_section
        rho_a = density_steel
        lambda_p = k_protection
        rho_p = density_protection
        d_p = thickness_protection
        A_p = perimeter_protected
        c_p = c_protection

        temperature_steel = time * 0.
        temperature_rate_steel = time * 0.
        specific_heat_steel = time * 0.

        temperature_steel[0] = temperature_fire[0]  # initially, steel temperature is equal to ambient
        temperature_steel_ = iter(temperature_steel)  # skip the first item
        next(temperature_steel_)  # skip the first item
        for i, v in enumerate(temperature_steel_):
            i += 1  # actual index since the first item had been skipped.
            specific_heat_steel[i] = c_steel_T(temperature_steel[i - 1] + 273.15)

            # Steel temperature equations are from [BS EN 1993-1-2:2005, Clauses 4.2.5.2]
            phi = (c_p * rho_p / specific_heat_steel[i] / rho_a) * d_p * A_p / V

            a = (lambda_p*A_p/V) / (d_p * specific_heat_steel[i] * rho_a)
            b = (v-temperature_steel[i-1]) / (1+phi/3)
            c = (np.exp(phi/10)-1) * (v-temperature_fire[i-1])
            d = time[i] - time[i-1]

            temperature_rate_steel[i] = a * b * d - c

            if temperature_rate_steel[i] < 0 < v-temperature_steel[i - 1]:
                temperature_rate_steel[i] = 0

            temperature_steel[i] = temperature_steel[i-1] + temperature_rate_steel[i]

        temperature_steel = integrate.cumtrapz(temperature_rate_steel, time, initial=0.)
        return temperature_steel, temperature_rate_steel, specific_heat_steel
