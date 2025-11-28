import gymnasium as gym
from gymnasium import spaces
from gymnasium.utils import seeding
import numpy as np
from os import path

class IonTrapEnv(gym.Env):
    metadata = {
        'render.modes': ['human', 'rgb_array'],
        'video.frames_per_second': 30
    }

    def __init__(self, max_xp, max_step, max_rabi):
        self.max_xp = max_xp
        self.max_step = max_step
        self.max_rabi = max_rabi
        self.viewer = None
        self.state = []

        self.lastmod = 0

        self.action_space = spaces.Box(
            low = -self.max_rabi,
            high = self.max_rabi,
            shape = (1,),
            dtype = np.float64
        )

        self.observation_space = spaces.Box(
            low = np.array([-self.max_xp, -self.max_xp, -self.max_xp, -self.max_xp,-self.max_xp, -self.max_xp, 0]),  # x_1, x_2, p_1, p_2, ..., step
            high = np.array([self.max_xp, self.max_xp, self.max_xp, self.max_xp, -self.max_xp, -self.max_xp, max_step]),  # x_1, x_2, p_1, p_2, ..., step
            dtype = np.float64
        )

        self.seed()

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]

    def step(self, state, fidelity=None):

        #### State ####
        x_0_op = state[0]
        x_1_op = state[1]
        x_2_op = state[2]

        p_0_op = state[3]
        p_1_op = state[4]
        p_2_op = state[5]

        idx_AM = int(state[6])
        if idx_AM == 1:
            self.lastmod = 0

        ###############

    # Two ions 100us (0,1), idx: 15
        # x_1_ideal = np.array([0.0,  0.0298697,  0.0596972, -0.0601572, -0.0295971,  9.60455e-05])
        # p_1_ideal = np.array([0.0, -0.0215547,  0.0703746,  0.0699092, -0.0217742,  1.80296e-05])
        # x_2_ideal = np.array([0.0,  0.1245390, -0.1926620,  0.1827030, -0.0935813,  0.000166019])
        # p_2_ideal = np.array([0.0, -0.0101025,  0.0677253, -0.0916269,  0.0824470, -0.000117070])
        # x_1_ideal = np.array([0, 0.00707896, 0.0192188, -0.0205744, -0.0370672, -0.0369609, 0.00721014, 4.09614e-06])
        # p_1_ideal = np.array([0, 0.0291264, 0.0413966, 0.0166446, 0.0164077, 0.0222981, -0.023131, 9.83529e-05])
        # x_2_ideal = np.array([0, 0.302222, 0.312891, 0.601597, 0.531422, 0.488772, 0.196258, -0.00041997])
        # p_2_ideal = np.array([0, -0.165335, -0.151002, 0.0649408, 0.101322, 0.211805, -0.0583048, -0.00281042])
        # x_1_ideal = np.array([0.0,  0.189874,  0.161518,   0.160948,   0.183087,  0.186927,  0.192969,  0.214025,  0.220543,  0.225959,  -1.86379e-05])
        # p_1_ideal = np.array([0.0,  0.097032, -0.0839933, -0.0834242, -0.0798346, -0.0722106, -0.0839703, -0.0872065, -0.0805947, -0.113448,  -0.000224995])
        # x_2_ideal = np.array([0.0,  0.304603,  0.603452,   0.602537,   0.561576,  0.55419,   0.528147,  0.516541,  0.535665,  0.528616,  -0.00161046])
        # p_2_ideal = np.array([0.0, -0.328006, -0.0941261, -0.0927328, -0.115465, -0.0992601, -0.108493, -0.065276, -0.0616529, 0.00754669, -0.00280923])

        # Three ions 100us (0,1), idx: 15. Ideal pulse: [-2.4111659, 1.53511, 1.21170287, -0.71186317, 1.20462246, 1.52825255, -2.42575168]
        # x_0_ideal = np.array([0, -0.0986125, -0.128232, -0.113253, -0.1039, -0.143181, -0.091569, 0.000103509])
        # p_0_ideal = np.array([0, -0.00936201, -0.0600901, -0.0155677, 0.0187963, 0.0430775, 0.00887957, 9.52618e-05])
        # x_1_ideal = np.array([0, -0.00200403, -0.00119132, -3.80114e-06, -0.00406654, -0.00356735, 0.00163083, 2.35733e-05])
        # p_1_ideal = np.array([0, 0.00593558, 0.00280489, -0.00046461, 0.00696915, 0.00726099, -0.00181514, 8.68313e-06])
        # x_2_ideal = np.array([0, 0.0895835, 0.00734111, 0.0255816, -0.0151648, -0.00542519, -0.0882542, -1.82573e-06])
        # p_2_ideal = np.array([0, -0.117316, -0.162957, -0.234533, -0.232153, -0.15785, -0.120987, -1.91697e-06])

        # Three ions 100us (0,1), 15-seg external pulse
        # x_0_ideal = np.array([0, 0.0089613, 0.0260773, 0.0282566, 0.0393904, 0.0236099, 0.00206759, 0.0282761, 0.04068, -0.00291574, -0.000923718, 0.0368168, 0.0336259, -0.00606589, 0.00253346, -0.000808727])
        # x_1_ideal = np.array([0, -0.0215542, 0.00793864, 0.0107448, -0.0126203, 0.0191821, -0.0127677, 0.00638934, -0.013339, -0.00325739, 0.000256483, -0.00331271, 0.00998749, -0.0172886, -0.00297816, 0.000125831])
        # x_2_ideal = np.array([0, -0.00640074, -0.0258108, -0.0237631, -0.00533599, 0.00240717, -0.00974109, -0.0375863, -0.073418, -0.116731, -0.12151, -0.0966121, -0.057018, -0.012264, 0.000278913, -0.000510825])

        # p_0_ideal = np.array([0, 0.0110738, -0.00887049, -0.00707625, -0.0224813, -0.0387389, 0.00200464, 0.00879178, -0.0135557, -0.0230845, 0.0142336, 0.0186242, -0.0305694, -0.0181521, -0.00162269, 0.000109795])
        # p_1_ideal = np.array([0, -0.000207575, -0.00669714, -0.00761737, 0.00586205, -0.0118461, 0.0212053, -0.013849, 0.020034, -0.0254097, 0.0287982, -0.0271989, 0.0276033, -0.0174508, -0.00252283, 1.52326e-05])
        # p_2_ideal = np.array([0, 0.0182516, 0.031878, 0.0322252, 0.0493105, 0.0832185, 0.0972573, 0.115748, 0.119181, 0.107896, 0.0559835, 0.0105775, -0.014737, -0.00664191, 0.00417213, 0.00112394])

        # Three ions 100us (0,1), 7-seg SOTA pulse
        x_0_ideal = np.array([0, 0.0344467, -0.00280103, 0.000215045, -0.00226058, 0.0124853, -0.0465305, 0.000123612])
        x_1_ideal = np.array([0, -0.00231988, -0.000502387, -0.000286582, 0.00056158, 0.0030073, 0.00185794, 0.00103087])
        x_2_ideal = np.array([0, -0.0533966, -0.142134, -0.138021, -0.139198, -0.139715, -0.0449758, -7.92084e-05])

        p_0_ideal = np.array([0, 0.0224145, -0.0302847, -0.0204367, -0.0217914, -0.0196347, 0.00426835, -0.000188727])
        p_1_ideal = np.array([0, 0.0111055, 0.00440646, 0.00384291, 0.00223324, -0.00318431, -0.000440032, 0.000501741])
        p_2_ideal = np.array([0, 0.0606039, 0.0128337, -0.00304387, -0.00378951, -0.0231765, -0.0615063, 4.02413e-05])

        # three ions idx: 50
        # x_0_ideal = np.array([0, -0.104418, -0.106074, -0.0575035, 0.068174, 0.11111, 0.109232, -1.69518e-05])
        # x_1_ideal = np.array([0, 0.0411166, 0.0362819, 0.135886, -0.029886, 0.0437944, 0.0415275, -2.13958e-06])
        # x_2_ideal = np.array([0, -0.0257709, -0.0251767, -0.00353069, -0.00195549, -0.0226942, -0.0236309, 2.97348e-06])
        #
        # p_0_ideal = np.array([0, -0.0509401, -0.0559934, -0.114343, -0.107297, -0.0449781, -0.0390344, 7.74698e-07])
        # p_1_ideal = np.array([0, -0.153187, -0.148437, -0.174492, -0.219875, -0.145812, -0.153449, 5.08192e-06])
        # p_2_ideal = np.array([0, 0.0201565, 0.0217199, 0.0176963, -0.018103, -0.024376, -0.0230504, 5.35318e-07])

        # print(1.0*x_0_op - 0.0*x_0_ideal[idx_AM], 1.0*x_1_op - 0.0*x_1_ideal[idx_AM], 1.0*x_2_op - 0.0*x_2_ideal[idx_AM])


# scale_list = [1.0, 1.0, 0.8, 0.5, 0.3, 0.5, 0.8, 1.0]
        # scale_factor = scale_list[idx_AM]

        # Three ions 200 us (0, 1), 7-seg SOTA pulse
        # x_0_ideal = np.array([0, -0.000242308, 0.0115521, 0.0115744, 0.0114022, 0.00820617, -0.000897929, 9.38797e-05])
        # x_1_ideal = np.array([0, -0.000704725, -0.00104002, 0.000687623, -0.00028563, -0.000422492, 0.0011323, 6.0275e-05])
        # x_2_ideal = np.array([0, 0.0553217, 0.120757, 0.111861, 0.108335, 0.11608, 0.0553354, 8.35302e-06])

        # p_0_ideal = np.array([0, 0.00693735, 0.0072649, 0.00476843, -0.00646198, -0.00584595, -0.00611867, -0.000334783])
        # p_1_ideal = np.array([0, 0.00128848, 0.00186075, -0.000821723, -6.3767e-05, -5.26411e-05, -0.000968606, -0.000145995])
        # p_2_ideal = np.array([0, -0.0140135, -0.0616403, -0.0432874, 0.0491735, 0.0705771, 0.0179935, -2.62138e-05])

        # modulus_1 = np.sqrt((x_0_op - x_0_ideal[idx_AM]) ** 2 + (p_0_op - p_0_ideal[idx_AM]) ** 2)
        # modulus_2 = np.sqrt((x_1_op - x_1_ideal[idx_AM]) ** 2 + (p_1_op - p_1_ideal[idx_AM]) ** 2)
        # modulus_3 = np.sqrt((x_2_op - x_2_ideal[idx_AM]) ** 2 + (p_2_op - p_2_ideal[idx_AM]) ** 2)
        # mod = 10*(1.0*np.sqrt(x_0_op**2 + p_0_op**2) + 1.0*np.sqrt(x_1_op**2 + p_1_op**2) + 1.0*np.sqrt(x_2_op**2 + p_2_op**2))
        mod = np.pi * (1.0*(x_0_op**2 + p_0_op**2)**1 + 9.0*(x_1_op**2 + p_1_op**2)**1 + 9.0*(x_2_op**2 + p_2_op**2)**1)
        # self.lastmod += mod

        # mod = 10.0 * (1.0 * modulus_1 + 1.0 * modulus_2 + 1.0 * modulus_3)
        costs = 1*np.abs(mod - np.pi/4/idx_AM)
        # print(self.lastmod)

        # costs = 100*abs_alpha*1/(7 - idx_AM + 1e-2)
        # costs = 100*(1-fid_rew)

        # if fidelity is not None:
        #     costs -= 20.0 * (fidelity)

        if fidelity is not None:

            if fidelity < 0.5:
                costs += 2.0 * (1 - fidelity)  # 5.0 * fidelity
            elif fidelity < 0.6:
                costs -= 1.0 * fidelity  # 1.0 * fidelity
            elif fidelity < 0.7:
                costs -= 2.0 * fidelity  # 2.0 * fidelity
            elif fidelity < 0.8:
                costs -= 3.0 * fidelity  # 3.0 * fidelity
            elif fidelity < 0.9:
                costs -= 4.0 * fidelity
            elif fidelity < 0.95:
                costs -= 5.0 * fidelity
            elif fidelity < 0.97:
                costs -= 6.0 * fidelity
            elif fidelity < 0.98:
                costs -= 7.0 * fidelity
            elif fidelity < 0.985:
                costs -= 9.0 * fidelity
            elif fidelity < 0.99:
                costs -= 10.0 * fidelity
            else:
                costs -= 12.0 * fidelity
        # if fidelity is not None:
        #
        #     if np.abs(fidelity-0.5) < 0.1:
        #         costs += 2.0 * (1 - fidelity)  # 5.0 * fidelity
        #     elif np.abs(fidelity-0.5) < 0.2:
        #         costs -= 1.0 * fidelity  # 1.0 * fidelity
        #     elif np.abs(fidelity-0.5) < 0.3:
        #         costs -= 2.0 * fidelity  # 2.0 * fidelity
        #     elif np.abs(fidelity-0.5) < 0.4:
        #         costs -= 3.0 * fidelity  # 3.0 * fidelity
        #     elif np.abs(fidelity-0.5) < 0.45:
        #         costs -= 4.0 * fidelity
        #     elif np.abs(fidelity-0.5) < 0.48:
        #         costs -= 7.0 * fidelity
        #     elif np.abs(fidelity-0.5) < 0.49:
        #         costs -= 9.0 * fidelity
        #     else:
        #         costs -= 11.0 * fidelity

        return self.state, -costs, True, {}

    def reset(self):
        self.state = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        return self.state

    def render(self, mode='human'):
        return None

    def close(self):
        if self.viewer:
            self.viewer.close()
            self.viewer = None
