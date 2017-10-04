# -*- coding: utf-8 -*-
import os
import sys
import matplotlib.pyplot as plt
from cls.timetemperaturefires import FireCurve
from cls.timetemperaturesteel import SteelTemperature

folder_path = os.path.realpath(__file__)
folder_path = os.path.dirname(folder_path)
print(folder_path)

# The standard fire_curve curve ISO 834
fire_curve = FireCurve(type_fire_curve='iso 834', time_end=0.5 * 60 * 60, time_step=12)

fire_curve.write_to_csv("/".join([folder_path, "test data.csv"]))
plt.plot(fire_curve.time, fire_curve.temperature)
plt.savefig('test figure.png')

heated_steel = SteelTemperature(fire_curve, "eurocode protected", zhidaole='a', haha='b')
# raise ValueError("error raising testing", "what?")