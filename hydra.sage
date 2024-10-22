from sage.all import *
from math import *

# might need to install pycryptodomex into sage:
# sage -sh
# pip3 install pycryptodomex
# exit
from Cryptodome.Hash import SHAKE128

class Hydra:
    def __init__(self, p, t, kappa):
        self.p = p
        self.F = GF(p)
        self.kappa = kappa
        self.perms = Hydra.num_perms(t)
        self.d = Hydra.get_d(p)
        self.Re_1, self.Re_2, self.Ri, self.Rh = Hydra.get_rounds(p, kappa, self.d)

        self.Me = matrix.circulant([self.F(3), self.F(2), self.F(1), self.F(1)])
        # generate matrices and vectors with conditions
        shake = Hydra.init_shake(p, "Matrices")
        self.Mi = self.gen_mi(shake)
        self.Mh = self.gen_mh(shake)

        # from here, just random values (without zero)
        shake = Hydra.init_shake(p, "Constants")
        self.rc_b = self.gen_rc(self.Re_1 + self.Re_2 + self.Ri, 4, shake)
        self.rc_h = self.gen_rc(self.Rh, 8, shake)

        # from here, random constants for the heads
        shake = Hydra.init_shake(p, "Rolling")
        self.rc_r = self.gen_rc(self.perms - 1, 8, shake)


    @staticmethod
    def get_R_star(kappa):
        assert(kappa >= 80)
        assert(kappa <= 256)
        R_star = [
            19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 23,
            23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26,
            26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 28, 28, 29, 29, 29, 29, 29,
            30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32,
            33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 34, 35, 35, 36, 36, 36, 36,
            36, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39,
            39, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42,
            43, 43, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 46, 46, 46, 46,
            46, 46, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49,
            49, 50, 50, 50, 50, 50, 50, 51, 51, 52, 52, 52, 52, 52, 53, 53, 53,
            53, 53, 53, 54, 54, 54, 54
        ]
        return R_star[kappa - 80]

    @staticmethod
    def get_round_num_head(p, kappa):
        R_star = Hydra.get_R_star(kappa)

        x0 = kappa / 24 + log(12, 2)
        x1 = (kappa - log(p, 2)) / 22 + log(11, 2)
        R_hat = 3 + ceil(max(x0, x1))
        R = ceil(1.25 * ceil(max(24, R_hat, R_star + 2)))
        return R

    @staticmethod
    def get_round_num_internal(p, kappa, d):
        x0 = kappa / float(16)
        x1 = (kappa - log(p, 2)) / float(12)
        R_hat = 4 - floor(log(d, 2)) + ceil(max(x0, x1))
        R = ceil(1.125 * ceil(max(kappa/4 - log(d, 2) + 6, R_hat)))
        return R

    @staticmethod
    def get_rounds(p, kappa, d):
        Re_1 = 2
        Re_2 = 4
        Ri = Hydra.get_round_num_internal(p, kappa, d)
        Rh = Hydra.get_round_num_head(p, kappa)
        return Re_1, Re_2, Ri, Rh

    @staticmethod
    def init_shake(p, string):
        bitlen = len(p.bits())
        num = ceil(bitlen / 64)

        shake = SHAKE128.new()
        shake.update('Hydra'.encode('ascii'))
        shake.update(string.encode('ascii'))
        for i in range(num):
            prime_block = (p >> (i * 64)) & ((0x1 << 64) - 1)
            shake.update(int(prime_block).to_bytes(8, byteorder = 'little'))
        return shake

    def field_element_from_shake(self, shake):
        bitlen = len(p.bits())
        byte = ceil(bitlen / 8)
        word = ceil(byte / 8)

        while True:
            word_buf = [0] * word
            buf = shake.read(byte)
            for i in range(word):
                byte_array = [0] * 8
                for j in range(i * 8, min((i + 1) * 8, byte)):
                    byte_array[j - i * 8] = buf[j];
                word_buf[i] = int.from_bytes(byte_array, byteorder = 'little')
            res = 0
            for el in reversed(word_buf):
                res = (res << 64) + el
            if res < self.p:
                return self.F(res)

    def field_element_from_shake_without_0(self, shake):
        while True:
            el = self.field_element_from_shake(shake)
            if el != self.F(0):
                return el

    @staticmethod
    def check_minpoly_condition(M, size):
        max_period = 2 * size
        M_temp = M
        for i in range(1, max_period):
            if not ((M_temp.minimal_polynomial().degree() == size) and (M_temp.minimal_polynomial().is_irreducible() == True)):
                return False
            M_temp = M * M_temp
        return True

    def check_conditions(self, lambdas, matrix, size):
        if not matrix.is_invertible() or not self.check_minpoly_condition(matrix, size):
            return False

        sum_ = self.F(0)
        for j in range(size):
            inner = self.F(0)
            for l in range(size):
                inner = inner + matrix[j, l]
            sum_ = sum_ + lambdas[j] * inner

        if sum_ == self.F(0):
            return False

        for j in range(size):
            sum_ = self.F(0)
            for l in range(size):
                sum_ = sum_ + lambdas[l] * matrix[l, j]
            if sum_ == self.F(0):
                return False

        return True

    def gen_matrix(self, shake, size):
        M = matrix(self.F, size, size, lambda _i, _j: self.F(1))
        M[0,0] = self.field_element_from_shake_without_0(shake)
        for i in range(1, size):
            M[i, 0] = self.field_element_from_shake_without_0(shake)
            M[i, i] = self.field_element_from_shake_without_0(shake)
        return M

    def gen_mi(self, shake):
        lambda1 = [self.F(1), self.F(-1), self.F(1), self.F(-1)]
        lambda2 = [self.F(1), self.F(1), self.F(-1), self.F(-1)]
        while True:
            M = self.gen_matrix(shake, 4)
            if self.check_conditions(lambda1, M, 4) and self.check_conditions(lambda2, M, 4):
                return M

    def gen_mh(self, shake):
        lambdas = [self.F(1), self.F(1), self.F(1), self.F(1), self.F(-1), self.F(-1), self.F(-1), self.F(-1)]
        while True:
            M = self.gen_matrix(shake, 8)
            if self.check_conditions(lambdas, M, 8):
                return M

    def gen_no_zero_sum(self, shake, size):
        v = vector(self.F, size)
        sum_ = self.F(0)
        while sum_ == self.F(0):
            for i in range(size):
                v[i] = self.field_element_from_shake_without_0(shake)
                sum_ = sum_ + v[i]
        return v

    def gen_lm(self, shake, size):
        lm = vector(self.F, size)
        sum_ = self.F(0)
        while sum_ == self.F(0):
            for i in range(size - 1):
                lm[i] = self.field_element_from_shake_without_0(shake)
                sum_ = sum_ + lm[i]
        lm[size - 1] = -sum_
        return lm

    @staticmethod
    def linear_independent(v1, v2):
        if (len(v1) != len(v2)) or len(v1) < 2:
            return False
        factor = v1[0] * v2[0].inverse_of_unit()
        for i in range(1, len(v1)):
            if v2[i] * factor != v1[i]:
                return True
        return False

    def gen_rc(self, R, t, shake):
        round_constants = []
        for _ in range(R):
            rc = vector(self.F, t)
            for i in range(t):
                rc[i] = self.field_element_from_shake(shake)
            round_constants.append(rc)
        return round_constants

    @staticmethod
    def get_d(p):
        for d in range(3, p):
            if gcd(d, p - 1) == 1:
                break
        return d

    @staticmethod
    def num_perms(t):
        t_ = t // 8
        t__ = t % 8

        perms = t_
        if t__ > 0:
            perms = perms + 1
        return perms

    def non_linear_e(self, state):
        for i in range(len(state)):
            state[i] = state[i] ** self.d
        return state

    @staticmethod
    def get_lm_dot(state):
        assert(len(state) == 4)
        tmp = state[0] - state[3]
        dot1 = tmp - state[1] + state[2]
        dot2 = tmp + state[1] - state[2]

        return dot1, dot2

    @staticmethod
    def non_linear_i(state):
        dot1, dot2 = Hydra.get_lm_dot(state)

        dot1 = dot1 * dot1
        sum_ = dot1 + dot2
        prod = sum_ * sum_

        for i in range(len(state)):
            state[i] = state[i] + prod
        return state

    @staticmethod
    def non_linear_h(state):
        assert(len(state) == 8)
        dot = state[0] + state[1] + state[2] + state[3] \
            - state[4] - state[5] - state[6] - state[7]
        dot = dot * dot

        out = vector([s + dot for s in state])
        return out

    @staticmethod
    def non_linear_r(y, z):
        vy, wy = Hydra.get_lm_dot(y)
        wz, vz = Hydra.get_lm_dot(z)

        v = vy * vz
        w = wy * wz

        y = [yi + v for yi in y]
        z = [zi + w for zi in z]

        return vector(y), vector(z)


    def R(self, state, i):
        assert(len(state) == 8)
        assert(len(self.rc_r) >= i)

        if i == 0:
            return state

        y = state[:4]
        z = state[4:]

        y, z = self.non_linear_r(y, z)
        y_perm = self.Mi * y
        z_perm = self.Mi * z

        state = self.concat_vec(y_perm, z_perm) + self.rc_r[i - 1]
        return state

    @staticmethod
    def concat(a, b):
        return vector(flatten([a, b]))

    @staticmethod
    def concat_vec(a, b):
        return vector(flatten([list(a), list(b)]))

    def permutation_b(self, state):
        sum_ = vector(self.F, 4)

        state = self.Me * state
        for i in range(self.Re_1):
            state = self.non_linear_e(state)
            state = self.Me * state + self.rc_b[i]
            sum_ = sum_ + state
        for i in range(self.Ri):
            state = self.non_linear_i(state)
            state = self.Mi * state + self.rc_b[i + self.Re_1]
            sum_ = sum_ + state
        for i in range(self.Re_1, self.Re_1 + self.Re_2):
            state = self.non_linear_e(state)
            state = self.Me * state + self.rc_b[i + self.Ri]
            if i < self.Re_1 + self.Re_2 - 1:
                sum_ = sum_ + state
        return state, sum_

    def permutation_h(self, state, K):
        for i in range(self.Rh):
            state = self.non_linear_h(state)
            state = self.Mh * state + self.rc_h[i]
            state = state + K
        return state

    def gen_ks(self, t,  K, IV, N):
        assert(len(IV) == 3)
        assert(len(K) == 4)
        state = Hydra.concat(N, IV)
        assert(len(state) == 4)
        K = vector(K)

        # first step
        state = state + K
        state, z = self.permutation_b(state)
        state = state + K

        # second step
        K_mat = self.Me * K
        K_ = Hydra.concat_vec(K, K_mat)

        keystream = []
        perm_counter = -1
        roll = Hydra.concat_vec(state, z)
        perm = None
        for i in range(t):
            off = i % 8
            if off == 0:
                perm_counter = perm_counter + 1
                roll = self.R(roll, perm_counter)
                perm = self.permutation_h(roll, K_)
                perm = perm + roll # feed forward
            keystream.append(perm[off])

        return keystream

    def encrypt(self, plains, K, IV, N):
        t = len(plains)
        keystream = self.gen_ks(t, K, IV, N)
        ciphers = []
        for i in range(t):
            ciphers.append(plains[i] + keystream[i])

        return ciphers

    def decrypt(self, ciphers, K, IV, N):
        t = len(ciphers)
        keystream = self.gen_ks(t, K, IV, N)
        plains = []
        for i in range(t):
            plains.append(ciphers[i] - keystream[i])

        return plains

    def get_f(self):
        return self.F

    def print_mat(self, mat, prefix):
        print("{} = [".format(prefix))
        for (count, row) in enumerate(mat):
            print("    [", end="")
            for i in range(len(row)):
                el = row[i]
                if el == self.F(-1):
                    el = -1
                elif el == self.F(-2):
                    el = -2
                if i == len(row) - 1:
                    print(el, end="")
                else:
                    print(el, end=", ")
            if count == mat.nrows() - 1:
                print("]")
            else:
                print("],")
        print("]")

    def print_lm_mat(self, mat, prefix):
        new_mat = mat
        new_mat[0, 0] = new_mat[0, 0] - 1
        for i in range(1, mat.nrows()):
            new_mat[i, 0] = new_mat[i, 0] - 1
            new_mat[i, i] = new_mat[i, i] - 1
        self.print_mat(new_mat, prefix)

    def print_conditionals(self):
        self.print_lm_mat(self.Mi, "self.Mi")
        self.print_lm_mat(self.Mh, "self.Mh")

p = 170141183460469231731687303715884105773
t = 2
sec = 128


hydra = Hydra(p, t, sec)

F = hydra.get_f()

MK = [F(0), F(0), F(0), F(0)]
IV = [F(1), F(1), F(1)]
N = F(2)

state_in  = [F(i) for i in range(t)]
state_out = hydra.encrypt(state_in, MK, IV, N)
state_out2 = hydra.decrypt(state_out, MK, IV, N)

assert(state_out2 == state_in)

for el in state_out:
    print((el))

# hydra.print_conditionals()
