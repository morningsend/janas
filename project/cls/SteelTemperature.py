# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy import integrate, gradient


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
    def __init__(self, type_steel, **kwargs):
        # Attributes, compulsory
        self._type_steel = type_steel
        self.__kwargs = kwargs                  # keyword arguments, varies for different 'fire_curve'

        # Attributes, output container
        self._df_output = None
        self.time = None
        self.temperature = None

        self.__check_parameters()

        self._df_output, self.time, self.temperature = self.__make_temperature_eurocode_protected_steel(
            time=self.__kwargs["time"],
            temperature_fire=self.__kwargs["temperature_fire"],
            rho_steel_T=self.__kwargs["density_steel"],
            c_steel_T=self.__kwargs["c_steel_T"],
            area_steel_section=self.__kwargs["area_steel_section"],
            k_protection=self.__kwargs["k_protection"],
            density_protection=self.__kwargs["density_protection"],
            c_protection=self.__kwargs["c_protection"],
            thickness_protection=self.__kwargs["thickness_protection"],
            perimeter_protected=self.__kwargs["perimeter_protected"]
        )[0:3]

    def __check_parameters(self):
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
                    "c_steel_T",
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

    def write_to_csv(self, folder_path):
        data_frame = pd.DataFrame(self._df_output)
        data_frame.to_csv(path_or_buf=folder_path, index=False)

    @staticmethod
    def __make_temperature_eurocode_unprotected_steel(
            time,
            temperature_fire,
            perimeter_section,
            area_section,
            perimeter_box,
            density_steel,
            c_steel_T,
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

        time_, temperature_steel[0], c_s[0] = iter(time), temperature_fire[0], 0.
        next(time_)
        for i, v in enumerate(time_):
            i += 1
            T_f, T_s_ = temperature_fire[i], temperature_steel[i-1]  # todo: steel specific heat
            c_s[i] = c_steel_T(temperature_steel[i - 1])

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
            rho_steel_T,
            c_steel_T,
            area_steel_section,
            k_protection,
            density_protection,
            c_protection,
            thickness_protection,
            perimeter_protected
    ):
        # todo: 4.2.5.2 (2) - thermal properties for the insulation material
        # todo: revise BS EN 1993-1-2:2005, Clauses 4.2.5.2

        params_required = [
            "time",
            "temperature_fire",
            "rho_steel_T",
            "c_steel_T",
            "area_steel_section",
            "k_protection",
            "density_protection",
            "c_protection",
            "thickness_protection",
            "perimeter_protected",
        ]
        
        V = area_steel_section
        rho_a = rho_steel_T
        lambda_p = k_protection
        rho_p = density_protection
        d_p = thickness_protection
        A_p = perimeter_protected
        c_p = c_protection

        temperature_steel = time * 0.
        temperature_rate_steel = time * 0.
        specific_heat_steel = time * 0.

        # Check time step <= 30 seconds. [BS EN 1993-1-2:2005, Clauses 4.2.5.2 (3)]
        time_change = gradient(time)
        # if np.max(time_change) > 30.:
            # raise ValueError("Time step needs to be less than 30s: {0}".format(np.max(time)))

        temperature_steel[0] = temperature_fire[0]  # initially, steel temperature is equal to ambient
        temperature_fire_ = iter(temperature_fire)  # skip the first item
        next(temperature_fire_)  # skip the first item
        for i, T_g in enumerate(temperature_fire_):
            i += 1  # actual index since the first item had been skipped.
            specific_heat_steel[i] = c_steel_T(temperature_steel[i - 1] + 273.15)

            # Steel temperature equations are from [BS EN 1993-1-2:2005, Clauses 4.2.5.2, Eq. 4.27]
            phi = (c_p * rho_p / specific_heat_steel[i] / rho_a(T_g)) * d_p * A_p / V

            a = (lambda_p*A_p/V) / (d_p * specific_heat_steel[i] * rho_a(T_g))
            b = (T_g-temperature_steel[i-1]) / (1.+phi/3.)
            c = (np.exp(phi/10.)-1.) * (T_g-temperature_fire[i-1])
            d = time[i] - time[i-1]

            temperature_rate_steel[i] = (a * b * d - c) / d  # deviated from e4.27, converted to rate (s-1)

            temperature_steel[i] = temperature_steel[i-1] + temperature_rate_steel[i] * d

            T_range_u = max([temperature_steel[i-1], T_g])
            T_range_l = min([temperature_steel[i-1], T_g])
            if temperature_steel[i] < T_range_l:
                temperature_steel[i] = T_range_l
                temperature_rate_steel[i] = (temperature_steel[i] - temperature_steel[i-1]) / d
            elif temperature_steel[i] > T_range_u:
                temperature_steel[i] = T_range_u
                temperature_rate_steel[i] = (temperature_steel[i] - temperature_steel[i-1]) / d

        data_dict_out = {
            "time [s]": time,
            "temperature fire [K]": temperature_fire,
            "temperature steel [K]": temperature_steel,
            "temperature rate steel [K/s]": temperature_rate_steel,
            "specific heat steel [J/kg/K]": specific_heat_steel
        }

        return data_dict_out, time, temperature_steel, temperature_rate_steel, specific_heat_steel