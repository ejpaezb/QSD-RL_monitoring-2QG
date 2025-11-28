import numpy as np

####

pulse_1_3_13 = 0.59736 * np.array([[0.147556151],[0.305590349],[0.458738592],[0.577224436],[0.707687278],[0.905632086],[1.0],[0.904234417],[0.70681344],[0.578407297],[0.459733025],[0.304708503],[0.146785477]])
delta_1_3_13 = 19.031604864060135

####

pulse_2_4_13 = 0.5747 * np.array([[0.154927247],[0.295281936],[0.462245408],[0.635467347],[0.761570425],[0.894149959],[1.0],[0.893655757],[0.761434094],[0.635123678],[0.462055112],[0.294863711],[0.15412843]])
delta_2_4_13 = 18.702549032429342

####

pulse_3_4_13_1 = 0.55624 * np.array([[0.132105096],[0.279775948],[0.44561626],[0.600517022],[0.744229638],[0.911799097],[1.0],[0.911624874],[0.743515531],[0.599913095],[0.445008202],[0.279236752],[0.132014885]])
delta_3_4_13_1 = 18.7051357863

pulse_3_4_13_2 = 0.55624 * np.array([[0.132105096],[0.279775948],[0.44561626],[0.600517022],[0.744229638],[0.911799097],[1.0],[0.911624874],[0.743515531],[0.599913095],[0.445008202],[0.279236752],[0.132014885]])
delta_3_4_13_2 = 18.609943722499885

pulse_3_4_11 = 0.43964 * np.array([[0.984490],[0.844742],[0.893510],[-1.192426],[-0.394975],[-1.535114],[-0.392193],[-1.219164],[0.865962],[0.831183],[1.000000]])
delta_3_4_11 = 18.609943722499885

pulse_3_4_9 = 0.55818 * np.array([[0.078327],[0.283847],[0.594395],[0.873738],[1.000000],[0.874213],[0.594913],[0.284201],[0.079001]])
delta_3_4_9 = 18.609943722499885

pulse_3_4_7 = 0.49873 * np.array([[0.519721],[0.785753],[0.791546],[1.000000],[0.791272],[0.785170],[0.519869]])
delta_3_4_7 = 18.609943722499885

####

pulse_3_5_13 = 0.61278 * np.array([[0.134832035],[0.295798819],[0.458590608],[0.583709332],[0.719155517],[0.915128836],[1.0],[0.913087886],[0.715032714],[0.580618733],[0.455101202],[0.291577366],[0.133032294]])
delta_3_5_13 = 18.7077281952

pulse_3_5_11 = 0.5539 * np.array([[1.000000],[0.784355],[0.658355],[-1.332440],[-0.303626],[-1.175303],[-0.357502],[-1.364903],[0.604666],[0.749940],[0.948044]])
delta_3_5_11 = 18.612618617192144

pulse_3_5_9 = 0.70235 * np.array([[0.124090],[0.334539],[0.611812],[0.858269],[1.000000],[0.865024],[0.623616],[0.344282],[0.128906]])
delta_3_5_9 = 18.612618617192144

####

def detuning_calc(mu1, mu2, r):
    
    det = (r * mu2 - mu1) / (r - 1.0)
    
    return det


#### Archive

# from scipy.optimize import rosen, differential_evolution

# bounds = [(1.0, 20.0)]*7
# result = differential_evolution(objective, bounds)
# result.x, result.fun

####

# def objective(Mass):
    
#     global M
#     M = Mass
        
#     print("M =", M)
        
#     ## Part 1
#     while True:
#         try:
#             ans = newton_krylov(equations, (np.random.rand(N) - 0.5) * 2.0)

#             # print("ans", ans)
#             # print("sum(equations(ans))", sum(equations(ans)))

#             break
#         except:
#             pass

#     ## Part 2
#     A_nm = np.zeros([N, N])

#     for n in range(N):
#         for m in range(N):
#             if (n == m):
#                 A_nm[n, m] = M[m] + 2.0 * sum([1.0 / np.abs(ans[m] - ans[p])**3 for p in [cnt for cnt in range(N) if cnt != m]])
#             else:
#                 A_nm[n, m] = -2.0 / np.abs(ans[m] - ans[n])**3

#     w, v = la.eig(A_nm)

#     # sho_mode = 0
#     # print(w[sho_mode])
#     # print(v[:,sho_mode])
    
#     ## Part 3
#     error = 0.0
#     for idx in range(N):
#         error += np.sum(np.abs(np.abs(v[:,idx][1:-1]) - gki[idx][1:-1]))
    
#     print("error", error)
#     return error

# objective([1,1,1,1,1,1,1])

####

# from scipy.optimize import fsolve, anderson, newton_krylov
# from math import exp

# N = 7
# M = [0.56, 1.0, 1.8466539956891497, 4.053690674600955, 8.217777777777778, 9.729227630271952, 11]


# def equations(variables):
    
#     equations = []
#     for m in range(N):
#         equations.append(M[m] * variables[m] - sum([1.0 / (variables[m] - variables[n])**2 for n in range(0, m)]) + sum([1.0 / (variables[m] - variables[n])**2 for n in range(m + 1, N)]))
            
#     return equations

# while True:
#     try:
#         ans = newton_krylov(equations, (np.random.rand(N) - 0.5) * 2.0)
    
#         print("ans", ans)
#         print("sum(equations(ans))", sum(equations(ans)))
        
#         break
#     except:
#         pass

####

# from numpy import linalg as la

# # Axial modes
# A_nm = np.zeros([N, N])

# for n in range(N):
#     for m in range(N):
#         if (n == m):
#             A_nm[n, m] = M[m] + 2.0 * sum([1.0 / np.abs(ans[m] - ans[p])**3 for p in [cnt for cnt in range(N) if cnt != m]])
#         else:
#             A_nm[n, m] = -2.0 / np.abs(ans[m] - ans[n])**3

# w_A, v_A = la.eig(A_nm)

# # Transverse modes
# K_nm = np.zeros([N, N])

# for n in range(N):
#     for m in range(N):
#         if (n == m):
#             K_nm[n, m] = (60 + 0.5) -0.5 * A_nm[n, m]
#         else:
#             K_nm[n, m] = -0.5 * A_nm[n, m]

# w_K, v_K = la.eig(K_nm)

# # print(np.argmin(w_A), w_A)
# # print(np.argmax(w_K), w_K)

# # sho_mode = 1
# # print(w[sho_mode])
# # print(v[:,sho_mode])

####

# m_0_scale =  1.0 * np.array([ 1, 1, 1, 1, 1]) # 1
# m_1_scale =  1.0 * np.array([ 1, 1, 1, 1, 1]) # 2
# m_2_scale =  1.0 * np.array([ 1, 1, 1, 1, 1]) # 3
# m_3_scale =  1.0 * np.array([ 1, 1, 1, 1, 1]) # 4
# m_4_scale =  1.0 * np.array([ 1, 1, 1, 1, 1]) # 5
# m_5_scale =  1.0 * np.array([ 1, 1, 1, 1, 1]) # 6
# m_6_scale =  1.0 * np.array([ 1, 1, 1, 1, 1]) # 0

# mode  = 1
# scale = m_1_scale
# bias  = 0.0

# sho_mode = np.argmax(w_K) + 1

# # plt.plot([2,3,4,5,6], v_A[:,sho_mode][1:-1], "o--", label="SHO")
# # plt.plot([2,3,4,5,6], v_K[:,sho_mode][1:-1], "o--", label="SHO")
# plt.plot([2,3,4,5,6], g_ki_exp[mode][1:-1] * scale + bias, "o--", label="Exp.")
# plt.legend(loc=0)