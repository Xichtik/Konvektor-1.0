# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 20:47:29 2021

@author: Admin
"""

import pandas

data = pandas.read_csv("raw.csv", header=0)
col_a = list(data.Pressure)
col_b = list(data.Altitude)
col_c = list(data.Temperature)
col_d = list(data.Dew_point)
col_e = list(data.Wind_direction)
col_f = list(data.Wind_speed)
print(str(col_a))
print(str(col_b))
print(str(col_c))
print(str(col_d))
print(str(col_e))
print(str(col_f))