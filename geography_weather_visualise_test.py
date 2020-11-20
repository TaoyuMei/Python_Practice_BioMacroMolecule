#!/usr/bin/python
#the file exam_test.py contains test code for the individual questions

import numpy as np
import matplotlib.pyplot as plt

#Part 2: Loading the data into Python

import exam
country_to_coordinates_map = exam.create_country_to_coordinates_map("country-capitals-final.csv")
print(country_to_coordinates_map["Denmark"])

continent_to_countries_map = exam.create_continent_to_countries_map("country-capitals-final.csv")

exam.plot_continent_countries(continent_to_countries_map)
plt.savefig("continent_countries.png")

temperature_data = exam.read_temperature_data("gistemp1200_ERSSTv5.csv")
print(temperature_data[:3,:])

domain = exam.get_domain(temperature_data[:,:3])

#Part 3: Looking up temperature values.

temperature_data_3d = exam.create_3d_array(temperature_data, domain)

print(exam.lookup_by_coordinates(temperature_data_3d,
                           domain,
                           latitude=56.0,
                           longitude=13.0,
                           year=2017))

print(exam.lookup_by_country(temperature_data_3d,
                             domain,
                             country_to_coordinates_map,
                             country="Denmark"
                             ))

from exam import TemperatureDataFinder                             
temperature_data_finder = TemperatureDataFinder('gistemp1200_ERSSTv5.csv',
                                                'country-capitals-final.csv')
print(temperature_data_finder.lookup_by_coordinates(56, 13, 2017))
print(temperature_data_finder.lookup_by_country("Denmark", 2017))
print(temperature_data_finder.lookup_by_coordinates(56, 13))
print(temperature_data_finder.lookup_by_country("Denmark"))

temperature_data_3d_fast = exam.create_3d_array_fast(temperature_data, domain)

#part 4: Visualizing temperature anomalies over time.

country_to_index_map = exam.create_country_to_index_map(continent_to_countries_map)

(time_indices,
 country_indices,
 temperature_anomaly_values
 ) = exam.create_scatter_plot_input(temperature_data_finder,country_to_index_map)
 
exam.generate_scatter_plot(time_indices,country_indices,temperature_anomaly_values)
plt.savefig("scatter_plot1.png")

exam.generate_scatter_plot(time_indices,country_indices,temperature_anomaly_values)
xtick_inputs = exam.generate_xtick_inputs(min_year=domain[0][0],
                           max_year=domain[0][1])
plt.xticks(xtick_inputs[0], xtick_inputs[1])
plt.savefig("scatter_plot2.png")

exam.generate_scatter_plot(time_indices,country_indices,temperature_anomaly_values)
xtick_inputs = exam.generate_xtick_inputs(min_year=domain[0][0],
                           max_year=domain[0][1])
plt.xticks(xtick_inputs[0], xtick_inputs[1])
ytick_inputs = exam.generate_ytick_inputs(continent_to_countries_map)
plt.yticks(ytick_inputs[0], ytick_inputs[1])
plt.savefig("scatter_plot3.png")
