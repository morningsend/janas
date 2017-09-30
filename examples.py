# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt

from cls.timetemperature import FireCurve

folder_path = os.path.realpath(__file__)
folder_path = os.path.dirname(folder_path)
print(folder_path)

# The standard fire curve ISO 834
fire = FireCurve(kind_fire_curve='eurocode hydrocarbon', time_end=3*60*60)

fire.write_to_csv("/".join([folder_path, "test data.csv"]))

plt.plot(fire.time, fire.temperature)
plt.savefig('test figure.png')