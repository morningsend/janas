# -*- coding: utf-8 -*-
import os
import sys
import matplotlib.pyplot as plt
from dat.property_steel import thermal as SteelProperty
from cls.timetemperaturefires import FireCurve
from cls.timetemperaturesteel import SteelTemperature

folder_path = os.path.realpath(__file__)
folder_path = os.path.dirname(folder_path)
print(folder_path)

# The standard fire_curve curve ISO 834
fire_curve = FireCurve(
    "travel",
    time_start=0.,
    time_end=22080,
    time_step=60.,
    temperature_ambient=20+273.15,
    fire_load_density=900e6,
    heat_release_rate_density=0.15e6,
    speed_fire_spread=0.012,
    distance_compartment_length=150,
    distance_compartment_width=17.4,
    area_ventilation=190,
    distance_ventilation_height=3.3,
    distance_element_to_fuel_bed_height=3.5,
    distance_element_to_fire_origin=105
)

fire_curve.write_to_csv("/".join([folder_path, "fire.csv"]))
plt.plot(fire_curve.time, fire_curve.temperature-273.15)
plt.savefig('fire.png')

protected_steel = SteelTemperature(
    "eurocode protected",
    time=fire_curve.time,
    temperature_fire=fire_curve.temperature,
    density_steel=SteelProperty("density"),
    c_steel_T=SteelProperty("specific heat capacity"),
    area_steel_section=0.002,
    k_protection=0.2,
    density_protection=800,
    c_protection=1700,
    thickness_protection=0.01712,
    perimeter_protected=0.3
)
protected_steel.write_to_csv("/".join([folder_path, "steel.csv"]))
plt.plot(protected_steel.time, protected_steel.temperature-273.15)
plt.savefig("steel.png")