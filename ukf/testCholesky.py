import scipy
import time
import numpy as np
import timeit


# bigAverage = 0.0
#     trials = 10

#     for x in range(trials):
#         bigAverage += run_ukf_textfile(start, cov, r, q, filename)

    
#         bigAverage = (bigAverage / trials)
#     print("difference between data and state over {} trials: {}".format(trials, bigAverage))
#     #what we have in dev currently: 7.56
#     #starting cov * .1: 8.51
#     #unit test ukf function and switch over to numpy arrays. python generators?


# def run_ukf_textfile(start, cov, r, q, filename):
#     '''
#     run_ukf_textfile
#         runs and visualizes UKF algorithm on input data file

#     @params
#         start: initialized state (1 x n)
#         cov: initialized covariance matrix (n x n)
#         r: noise vector for predictions (1 x n)
#         q: noise vector for sensors (1 x m)
#         filename: text file of cleaned sensor data to read from (any length)
#     '''
#     f = open(filename, "r")
#     data = f.readline()
#     splitData = data.split(",")
#     splitData = [float(x) for x in splitData]
#     # start[0] = splitData[0]
#     # start[1] = splitData[1]
#     # start[2] = splitData[2]
#     # start[3] = splitData[3]
#     i = 1
#     u_k = []
#     estimate_difference = 0.0
#     prevAve = 0.0
#     while(data):
#         # u_k = dataFromGPS()
#         start, cov = UKF_algorithm.UKF(start, cov, r, q, u_k, splitData)
#         # game_visualize(np.array([start[:4]]), i)

#         data = f.readline()
#         if(data == ''):
#             break
#         splitData = data.split(",")
#         splitData = [float(x) for x in splitData]
#         i+=1

#         aveDiff = 0.0

#         for x in range(6):
#             stepDiff = 0.0
#             if(splitData[x] == 0):
#                 stepDiff = prevAve
#             else:
#                 stepDiff = (start[x + 1] - splitData[x]) / abs((splitData[x]))
#             aveDiff += stepDiff

#         aveDiff = (aveDiff / 6)
#         prevAve = aveDiff
#         # print("difference at step {}: {}%".format(i, aveDiff))

#         estimate_difference += aveDiff

#     estimate_difference = (estimate_difference / i)
#     print("difference between data and state: {}".format(estimate_difference))
#     f.close()
#     return estimate_difference


# tot = 0
# iterations = 1000
# n = 10
# scaling = 3-n
# #100 iterations for original sigma: .00064
# for i in range(iterations):
#     start=time.time()
#     means = np.random.rand(1, n)
#     cov = np.random.rand(10,10)
#     sigma(means, cov, n , scaling)
#     tot += (time.time()-start)

# print(tot/iterations)
# print(timeit.timeit("sigma(means, cov, n, scaling)", number=100))


test_code = '''
import numpy as np
import scipy
n = 10
scaling = 3-n
means = np.random.rand(1, n)
cov = np.random.rand(10,10)

sigmaMatrix = np.zeros((2*n+1,n))
temp = np.zeros((n, n))

# for i in range(len(cov)):
#     for j in range(len(cov[i])):
#         temp[i][j] = cov[i][j]  * (n + scaling)
temp = np.multiply(cov, (n + scaling))

temp = scipy.linalg.sqrtm(temp)
sigmaMatrix[0] = means

# for i in range(n):
#     sigmaMatrix[i + 1] = np.add(means, temp[i])  # mean + covariance
#     sigmaMatrix[i + 1 + n] = np.subtract(means, temp[i])  # mean - covariance
sigmaMatrix[1:(n+1)] = np.add(means, temp)
sigmaMatrix[(n+1):(2*n+1)] = np.subtract(means, temp)

'''

print(timeit.repeat(test_code, number=2000, repeat=10))
# print(timeit.timeit(test_code, number=2000))
#2000 runs: 1.1 ish for original


''' testing for equivalence between two methods '''
# n = 10
# scaling = 3-n
# means = np.random.rand(1, n)
# cov = np.random.rand(10,10)


# sigmaMatrix = np.zeros((2*n+1,n))
# temp = np.zeros((n, n))

# for i in range(len(cov)):
#   for j in range(len(cov[i])):
#     temp[i][j] = cov[i][j]  * (n + scaling)

# temp = scipy.linalg.sqrtm(temp)
# # first column of sigma matrix is means
# sigmaMatrix[0] = means

# # traverse N (10) dimensions, calculating all sigma points
# for i in range(n):
#     sigmaMatrix[i + 1] = np.add(means, temp[i])  # mean + covariance
#     sigmaMatrix[i + 1 + n] = np.subtract(means, temp[i])  # mean - covariance
# print("First sigma matrix: ", sigmaMatrix)


# sigmaMatrix2 = np.zeros((2*n+1,n))
# temp = np.zeros((n, n))

# # sigma point formula from website
# temp = np.multiply(cov, (n + scaling))

# temp = scipy.linalg.sqrtm(temp)
# # first column of sigma matrix is means
# sigmaMatrix2[0] = means

# # traverse N (10) dimensions, calculating all sigma points
# sigmaMatrix2[1:(n+1)] = np.add(means, temp)
# sigmaMatrix2[(n+1):(2*n+1)] = np.subtract(means, temp)
# print("Second sigma matrix: ", sigmaMatrix2)
# print(np.array_equal(sigmaMatrix2, sigmaMatrix))

''' np.array testing '''
# test = np.array([[]])
# test = np.zeros((5, 2))
# test = np.array([[1, 2], [3, 4], [5,6], [7, 8], [9, 10]])
# temp = np.array([[1, 2], [3, 4]])
# add = np.array([1, 2])
# print(np.add(add, temp))
# # 1:n+1
# print(test)
# test[1:3] = np.add(add, temp)
# test[3:5] = np.subtract(add, temp)
# print(test)
# print(np.add(add, temp[0]))



