# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def thermal(property_name):

    property_file = {
        "density": "rho_1_T_steelc_ec.csv",  # BS EN 1993-1-2:2005, 3.2.2
        "thermal conductivity": "k_1_T_steelc_ec.csv",  # BS EN 1993-1-2:2005, 3.4.1.3
        "specific heat capacity": "c_1_T_steelc_ec.csv",  # BS EN 1993-1-2:2005, 3.4.1.2
        "reduction factor for the slope of the linear elastic range": "kE_1_T_steelc_ec.csv",  # BS EN 1993-1-2:2005, Table 3.1
        "reduction factor for proportional limit": "kp_1_T_steelc_ec.csv",  # BS EN 1993-1-2:2005, Table 3.1
        "reduction factor for effective yield strength": "ky_1_T_steelc_ec.csv"  # BS EN 1993-1-2:2005, Table 3.1
    }

    # read file
    data = pd.read_csv("./dat/"+property_file[property_name], delimiter=",", dtype=float)
    x, y = data.values[:, 0], data.values[:, 1]

    return interp1d(x, y)

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



temperature = np.arange(20.,1200.+0.5,0.5) + 273.15
prop = temperature * 0
for i,v in enumerate(temperature):
    prop[i] = steel_specific_heat_carbon_steel(v)
df = pd.DataFrame(
    {
        "Temperature [K]": temperature,
        "Specific Heat Capacity [J/kg/K]": prop
    }
)
df.to_csv("c_1_T_steelc_ec.csv", index=False, columns=["Temperature [K]", "Specific Heat Capacity [J/kg/K]"])

