'''
Subtractive Clustering Algorithm according to "Learning of fuzzy rules by mountain clustering" by Ronald R. Yager and Dimitar P. Filev.
'''
__author__ = 'Alex Povod'
__version__ = '0.1.1'

import numpy
import numpy.typing as npt
import scipy as sp


class SubtractiveClustering:

    def __init__(self, r_1:int, r_2:float, eps = 0.5):
        '''
        Implemetation Subtractive Clustering Method for computing the cluster centroids, 
        which belongs to unsupervised learning and can quickly determine the number of clusters and 
        cluster centroids based on the raw data.
        :param: r_1:    is the naighbor radius, which influences the scope of a cluster centroid.
        :param: r_2:    is the influencig weight of the last cluster centroid.
                        To avoid getting close cluster centroids, in general r_2 = r_1 * 1.5
        :param: eps:    is the accuracy ratio, positive constant less than 1.
        '''
        
        self.r_1 = r_1; self.r_2 = r_2; self.eps = eps    
        
    def fit(self, data: npt.NDArray):
        
        # Step 1: Compute the Mountain Function
        mountain_function = self.compute_mountain_function(data)
        self.potential = mountain_function

        # Init array of number centroid clusters
        self.centers = numpy.empty((0, data.shape[1]))
        # Step 2: Select maximum mountain function
        max = numpy.max(mountain_function)
        max_star = 0
        while(max >= self.eps * max_star):
            
            center = data[numpy.argmax(mountain_function)]
            self.centers = numpy.vstack((self.centers, center))
            # Step 3: Update the Mountain Function
            mountain_function = self.update_mountain_function(mountain_function, max, center, data)

            # Step 4: Update maximum mountain fuction and check in loop
            max_star = max
            max = numpy.max(mountain_function)
            
        return self
        

    def compute_mountain_function(self,data: npt.NDArray) -> npt.NDArray:
        '''
        Step 1 of algoritm compute the potential (mountain function) for each point or sample.
        returns: the potential for that given row given no previous potential known
        '''
        potential = numpy.array([])
        i = 0
        while(i < len(data)):
            ls = numpy.array([])
            for ii in range(0, len(data)):
                numerator = numpy.power((data[i] - data[ii]), 2)
                denomerator = numpy.power((self.r_1 / 2),2)
                ratio = numpy.exp(numpy.negative(numerator.sum() / denomerator))
                ls = numpy.append(ls, ratio).sum()
            i += 1
            potential = numpy.append(potential, ls)
            
        # alternative...
        # potential = numpy.sum(numpy.exp(numpy.negative( sp.spatial.distance_matrix(data, data) / numpy.power((self.r_1 / 2), 2) )), axis=1)   
        return potential

    
    def update_mountain_function(self, mountain_function, max_mountain_function, center, data):
        '''
        Step 3 of algorithm update the mountain function of each data.
        '''
        numerator = numpy.linalg.norm(data - center, axis = 1)
        denomerator = numpy.power((self.r_2 / 2), 2)
        ratio = numpy.exp(numpy.negative(numerator / denomerator))
        return mountain_function - max_mountain_function * ratio

        