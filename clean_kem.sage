from sage.all import Matrix
from fpylll import LLL, IntegerMatrix
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
from sage.matrix.constructor import random_unimodular_matrix
import hashlib, random
import os
import sys
import copy

# PARAMS
#dimension of lattices
if(len(sys.argv[1:]) == 0):
    ndim = 75
else:
    ndim = int(sys.argv[1])
print("lattice dimensions (ndim): ", ndim)

rho = 0.475 # rho < lambda_1(S)/2 = 1/2
print("PARAMS:\nrho: ", rho)

# Define the standard deviation parameter 's', must be larger than smoothing parameter
# realistically, s is a value from 0.5 in the lowest dimension, and almost 1 in dim 200
# s is either larger than lambda_n(S)=1 or larger than ||B*||*ln(2n+4). since the second case is lower than 1 for dim <200, lambda_n(S) is the max value. 
s = 100*rho/sqrt(ndim)#2**(ndim) # Choose a positive real value # from corollary 3.3 HuckBennett's paper on rotations of Zn
print("s: ", s)

q = ceil((s * ndim / rho) * sqrt(log(2*ndim + 4)/pi)) #from Ducas' KEM TODO: delete q * 100000 from sampler line, just for DEBUG
print("q: ", q) 

# Define the MatrixSpace over RR (Real Field with 53 bits of precision)
MS_RR_c = MatrixSpace(RR, 1, ndim)
MS_RR = MatrixSpace(RR, ndim, ndim)
MS_ZZ = MatrixSpace(ZZ, ndim, ndim)
MS_ZZ_c = MatrixSpace(ZZ, 1, ndim)
MS_ZZ_c_transposed = MatrixSpace(ZZ, ndim, 1)
MS_QQ = MatrixSpace(QQ, ndim, ndim)
MS_QQ_c = MatrixSpace(QQ, 1, ndim)

# END PARAMS


# AUX FUNCTIONS

def matrix_subtract(m1, m2):
    print("begin matrix subtraction")
    print(m1, " M1!")
    print(m2, " M2!")
    res = m1 - m2
    print(res, " RES!")
    print("end matrix subtraction")
    return res

#Let x be a vector in Zn. The norm of x over the quadratic form Q is ||x||_Q = sqrt(x^t * Q * x)
def quadratic_norm(x, Q):
    return sqrt((x.transpose() * Q * x)[0][0])

def restart_program():
    """Restarts the current program."""
    print("Retrying...")
    python = sys.executable  # Path to the Python interpreter
    os.execl(python, python, *sys.argv)  # Replace the current process with a new one

#    UNIVERSAL HASH FUNCTION VERSION W/OUT RANDOM ORACLE
#    # define random seed z
#    Z = random.randint(0, s)
#
#    # vec_e has dim = ndim, rand_z needs to be a computationally secure random binary sequence
#    def extractor(vec_e, rand_z):
#        print("inputs: vec_e, rand_z ::: ", vec_e, rand_z)
#        hash_obj = hashlib.sha512()
#        for i in range(len(vec_e)):
#            rand_z = rand_z + vec_e[i] # doubt: i decided summing the values, not sure there is a better alternative for a deterministic way of computing this
#        hash_obj.update(str(rand_z).encode('utf-8'))
#        return hash_obj.hexdigest()

#RANDOM ORACLE VERSION
# vec_e has dim = ndim
# hashes vector e
def extractor(vec_e):
    print("extractor input vector: ::: ", vec_e)
    hash_obj = hashlib.sha256()
    hash_obj.update(str(vec_e).encode('utf-8'))
    return hash_obj.hexdigest()

# Discrete Gaussian Sampler for vector e from KEM (step 2 Encaps)
## s = 500 works well for ndim=2; 
def discrete_n_dim_vector_sampler(q, rho, P):
    sampler = DiscreteGaussianDistributionLatticeSampler(P, (q  * rho) / sqrt(ndim))
    sample = sampler()
    e_vec = [0] * ndim

    if(type(sample) == type(None)):
        return restart_program()
    if(sample.is_zero()):
        return restart_program()
    print("sample ", sample)

    #copy to itself so the vector becomes mutable
    sample = copy.copy(sample)

    for k in range(ndim):
        e_vec[k] = sample[k] / q
    
    return e_vec

#returns a vector from S such that ||y - Uc||_S <= rho
#decoding in Z_n is equivalent to rounding to the nearest integer
def decode(S, U, c):
    Uc = U * c
    print("Uc:::::: ", Uc.transpose())

    y = Uc.apply_map(lambda x: round(x)) #if(x>=0):ceil(x) else:floor(x))
    print("y at decode: ", y.transpose())

    quad_norm = quadratic_norm((y.transpose() - Uc.transpose()).transpose(), S)
    print("quad_norm of y-(U*c) :  ", quad_norm)
    if(quad_norm <= rho):
        print("GOOD DECODING, DECODED: ", y.transpose())
        return y

    print("BAD QUADRATIC NORM FOR DECODING, RESTART")
    restart_program()

# END AUX FUNCTIONS

#GEN function from KEM
#Outputs:
#   S and P, quadratic forms/pub keys
#   U, unimodular matrix that acts as secret key
def gen():
    print("\n\n\nBEGIN KEYGEN::::\n")
    S = identity_matrix(ndim)

    # Obtain a random unimodular matrix in GL, reduce with BKZ or LLL
    matrix_space = sage.matrix.matrix_space.MatrixSpace(ZZ, ndim)
    U = random_unimodular_matrix(matrix_space).BKZ()

    print("QUADRATIC FORM S=I_n")
    print("PRIV/SECRET KEY, U")
    print(U)
    print("SECRET KEY DETERMINANT: ", U.determinant())

    # This is the P from the KEM
    # If S=I_n, then P =U^t * I_n * U = U^t * U = I_n, regardless of U 
    #P = U.transpose() * U
    P = S

    print("PUBLIC KEY/QUADRATIC FORM P = U^t * S * U = I_n")

    #secret_key = U
    #public_key = P
    return S, U, P


#ENCAPS function from KEM
#inputs: quadratic form P, which is the public key
#outputs: 
#   encaps_symmetric_k, which is the common secret 
#   c, clue for decapsulation
def encaps(P):
    print("\n\n\nBEGIN ENCAPS::::\n")
    ## sampling e vector
    e = discrete_n_dim_vector_sampler(q, rho, P) 
    print("e vector: ", e)

    # c vector definition: vector e over the discretized torus Tq
    c = copy.copy(e)
    for i in range(ndim):
        c[i] = c[i] - floor(c[i])
    print("c vector ::: ", c)

    encaps_symmetric_k = extractor(e)
    print("encaps_symmetric_k: ", encaps_symmetric_k)

    #old: from universal hash function model
    #encaps_ciphertext = [c, Z]
    
    #random oracle model
    encaps_ciphertext = [c]

    output_encaps = [encaps_ciphertext, encaps_symmetric_k]

    #print("((c, Z), k) :::: ", output_encaps)
    print("((c), k) :::: ", output_encaps)

    return encaps_symmetric_k, c


#DECAPS function from KEM
# inputs: secret_key = U; encaps_ciphertext = (c, Z)
# outputs: decaps_symmetric_k, a decapsulated symmetric key.
def decaps(S, U, c):
    print("\n\n\nBEGIN DECAPS::::\n")
    
    print("inputs: secret_key = U =\n", U, "\n; encaps_ciphertext = (c := (c)) : ", c, "\n")
    
    #Step 1: obtain y from S, Uc
    y = decode(S, MS_QQ(U), MS_QQ_c(c).transpose())

    #Step 2: compute k
    #U^-1 * y is the integer part of vector "e" from encaps
    U_invert_times_y = U.inverse() * y
    print("U_invert_times_y::::: ", U_invert_times_y.transpose())

    c_minus_Uy = MS_QQ_c(c) - U_invert_times_y.transpose()
    #convert c_minus_Uy to list for extractor
    c_minus_Uy_list = [0] * ndim
    for i in range(ndim):
        c_minus_Uy_list[i] = c_minus_Uy[0, i]
    print("c_minus_Uy ::::: ", c_minus_Uy_list)

    #OLD: #decaps_symmetric_k = extractor(c_minus_Uy_list, Z)
    decaps_symmetric_k = extractor(c_minus_Uy_list)

    return decaps_symmetric_k


def main():
    ### BEGIN KEYGEN
    S, U, P = gen()
    ### END KEYGEN
    
    ### BEGIN ENCAPS
    encaps_symmetric_k, c = encaps(P)
    ### END ENCAPS

    ### BEGIN DECAPS
    decaps_symmetric_k = decaps(S, U, c)
    ### END DECAPS

    print("decaps_symmetric_k: ", decaps_symmetric_k)
    print("is encaps_symmetric_k equal to decaps_symmetric_k?: ", decaps_symmetric_k == encaps_symmetric_k)

main()