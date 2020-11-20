#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

#Part 2: Loading the data into Python

def create_country_to_coordinates_map(filename):
    """read the file country-capitals-final.csv and return a dictionary
    with country names as keys and a tuple of (latitude, longitude)
    as values"""
    filehandle = open(filename, 'r')
    map_dict = {}
    filehandle.readline()   #skip the header of the csv table
    for line in filehandle:
        line = line.rstrip('\n')
        line_list = line.split(',')
        map_dict[line_list[0]] = (float(line_list[2]), float(line_list[3]))
    return map_dict

def create_continent_to_countries_map(filename):
    """read the file country-capitals-final.csv and return a dictionary
    with continent names as keys and a list of all countries on the continents
    as values"""
    filehandle = open(filename, 'r')
    map_dict = {}
    filehandle.readline()   #skip the header of the csv table
    for line in filehandle:
        line = line.rstrip('\n')
        line_list = line.split(',')
        continent = line_list[5]
        country = line_list[0]
        if not map_dict.get(continent):  #if the continent name hasn't been collected before
            map_dict[continent] = [country, ]
        else:
            if country in map_dict[continent]:  #if the country name has been collected before
                continue
            else:
                map_dict[continent].append(country)
    return map_dict

def plot_continent_countries(list_dict):
    """create a barplot showing the number of countries
    on each continent"""
    plt.figure()
    continent_list = list(list_dict.keys())
    country_number_list = []
    for continent in continent_list:
        country_number_list.append(len(list_dict[continent]))
    position = list(range(len(continent_list)))
    #the y axis should show more number
    plt.bar(x = position, height = country_number_list)
    plt.xticks(position, continent_list)
    
def read_temperature_data(temperature_filename):
    """read the GISSTEMP data, store the temperature in a numpy array"""
    temperature_array = np.genfromtxt(temperature_filename, delimiter=",", skip_header=2251000)
    #change to skip_header=1 before submitting
    return temperature_array

def get_domain(np_array):
    """calculate the the minimum, the maximum and the number of unique values of 
    each colunm of the given numpy array"""
    tuple_list = []
    for i in range(np_array[0]):    #range(np_array[0]) is the column number of the numpy array
        column_i = np_array[:,i]
        tuple_list.append((np.min(column_i), np.max(column_i), len(np.unique(column_i))))
    return tuple_list

#Part 3: Looking up temperature values.

def get_indices(data, domain):
    """For each column in data, translate values into indices,
       using the provided domain information"""

    # If only a 1d sequence is provided, convert it to 2d
    array = np.atleast_2d(data)

    # Create a container for indices
    indices = np.zeros(array.shape, dtype='int')

    # Iterate over column
    for j in range(array.shape[1]):
        # Use the (min, max, nvalues) information provided for each column
        bins = np.linspace(domain[j][0], domain[j][1], domain[j][2])

        # The bin values refer to the center of the bin - move to refer to
        # left boundary instead
        bins += 0.5*(bins[1]-bins[0])

        # Calculate index values for all values in column
        indices[:,j] = np.digitize(array[:,j], bins=bins)

    # If we have a 2D array consisting of only 1 column,
    # transform to 1D
    return indices.squeeze()
    
def create_3d_array(temperature_array, domain):
    """put the data in a 3D numpy array to support fast lookup;
    look up the temperature anomaly value for a particular coordinate;"""
    shape1 = np.array(domain,dtype='int')[:,2]  #decide the size of each dimention which is the number of unique values
    three_d_array = np.zeros(shape=shape1)
    indices = get_indices(temperature_array[:,:3], domain)
    for i in range(len(indices)):
        x = indices[i][0]
        y = indices[i][1]
        z = indices[i][2]
        three_d_array[x][y][z] = temperature_array[i][3]
    return three_d_array
    
def lookup_by_coordinates(three_d_array,domain,latitude,longitude,year=0):
    """get the corresponding temperature anomaly
    by specifying 3d temperature data, domain, latitude, longitude and the year"""
    if not year==0:
        index = get_indices(np.array([year,latitude,longitude]), domain)
        tempanomaly =three_d_array[index[0]][index[1]][index[2]]
        return tempanomaly
    else:   #if the year is not specified
        all_year_values = []
        for each_year in range(domain[0][2]):   #the total number of year
            index = get_indices(np.array([0,latitude,longitude]), domain)
            tempanomaly = three_d_array[each_year][index[1]][index[2]]
            all_year_values.append(tempanomaly)
        return all_year_values

def lookup_by_country(three_d_array,domain,country_to_coordinates_map,country,year=0):
    """lookup the temperature anomaly by country name"""
    lat_lon_tuple = country_to_coordinates_map[country]
    tempanomaly = lookup_by_coordinates(three_d_array,
                                        domain,
                                        latitude=lat_lon_tuple[0],
                                        longitude=lat_lon_tuple[1],
                                        year=year
                                        )
    return tempanomaly

class TemperatureDataFinder:
    def __init__(self, temperature_filename, country_capitals):
        temperature_array = read_temperature_data(temperature_filename)
        self.domain = get_domain(temperature_array[:,:3])
        self.three_d_array = create_3d_array(temperature_array, self.domain)
        self.country_to_coordinates_map = create_country_to_coordinates_map(country_capitals)

    def lookup_by_coordinates(self,latitude,longitude,year=0):
        tempanomaly = lookup_by_coordinates(self.three_d_array,self.domain,latitude,longitude,year)
        return tempanomaly

    def lookup_by_country(self,country,year=0):
        tempanomaly = lookup_by_country(self.three_d_array,
                                        self.domain,
                                        self.country_to_coordinates_map,
                                        country,
                                        year)
        return tempanomaly

def create_3d_array_fast(temperature_array, domain):
    """a faster way to put the data in a 3D numpy array
    without for loop but only use numpy operations"""
    shape1 = np.array(domain,dtype='int')[:,2]  #decide the size of each dimention which is the number of unique values
    three_d_array = np.zeros(shape=shape1)
    indices = get_indices(temperature_array[:,:3], domain)
    indices_transpose = [indices[:,0],indices[:,1],indices[:,2]]
    indices_raveled = np.ravel_multi_index(indices_transpose,dims=shape1)
    temperature = temperature_array[:,3]
    index_not_used = np.zeros(shape=(1,np.min(indices_raveled)))
    fill_in_values = np.concatenate((index_not_used[0],temperature))    #values in order
    three_d_array.flat = fill_in_values #traverse the empty array and assign value to each element at once
    return three_d_array
    
#part 4: Visualizing temperature anomalies over time.

def create_country_to_index_map(continent_to_countries_map):
    """creat a dictionary that maps each country to an index, 
    the indices are natural number beginning with 0"""
    country_to_index_map = {}
    country_index = 0
    continent_list = list(continent_to_countries_map.keys())
    continent_list.sort()
    for continent in continent_list:
        country_list = continent_to_countries_map[continent]
        country_list.sort()
        for country in country_list:
            country_to_index_map[country] = country_index
            country_index += 1
    return country_to_index_map
    
def create_scatter_plot_input(TemperatureDataFinder,country_to_index_map):
    """create three lists, i.e. year indices, country indices and temperature anomaly
    to facilitate drawing the scatter plot"""
    year_index_list = []
    country_index_list = []
    tempanomaly_list = []
    year_number = TemperatureDataFinder.domain[0][2]
    year_index_each_country = list(np.array(range(year_number)))
    for country in list(country_to_index_map.keys()):
        year_index_list += year_index_each_country
        country_index_repeat = [country_to_index_map[country]]*year_number
        country_index_list += country_index_repeat
        tempanomaly_every_year = TemperatureDataFinder.lookup_by_country(country)
        tempanomaly_list += tempanomaly_every_year
    return (year_index_list,country_index_list,tempanomaly_list)

def generate_scatter_plot(year_index_list,country_index_list,tempanomaly_list):
    """generate a scatter plot with year as x-axis and country as y-axis;
    temperature anomaly above 0 is shown in red dot, below zero blue, zero white"""
    plt.figure()
    plt.scatter(x=year_index_list,
                y=country_index_list,
                c=tempanomaly_list,
                cmap="RdBu_r")
                
def generate_xtick_inputs(min_year, max_year):
    """generate the list of positions and the list of labels
    for the plt.xticks() function; get a tick every 20 years"""
    number_of_20 = (max_year-min_year+1)//20
    position_list = [0]
    year_list = [int(min_year)]
    year = int(min_year)
    i=1
    while(i <= number_of_20):
        position_list += [20*i]
        year += 20
        year_list += [year]
        i += 1
    return (position_list, year_list)

def generate_ytick_inputs(continent_to_countries):
    """generate two list,
    i.e. the center y-value for each continent
    and the corresponding continent names"""
    center_y_list = []
    y_value = 0 #begin with the zero point of y-axis
    continent_list = list(continent_to_countries.keys())
    continent_list.sort()
    for continent in continent_list:
        country_number = len(continent_to_countries[continent])
        y_value += country_number//2
        center_y_list += [y_value]
        y_value += country_number//2 + country_number%2
        #make sure that y_value corresponds to the first country on the next continent
    return (center_y_list, continent_list)



    
