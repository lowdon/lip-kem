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
ndim = 3


rho = 0.475 # rho < lambda_1(S)/2 = 1/2
print("PARAMS:\nrho: ", rho)

# Define the standard deviation parameter 's', must be larger than smoothing parameter
# realistically, s is a value from 0.5 in the lowest dimension, and almost 1 in dim 200
# s is either larger than lambda_n(S)=1 or ||B*||*ln(2n+4). since the second case is lower than 1 for dim <200, lambda_n(S) is the max value. 
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
    #x = MS_ZZ_c_transposed(x)
    return sqrt((x.transpose() * Q * x)[0][0])

#Outputs the lambda_n_of_quadratic Q
def lambda_n_of_quadratic(Q):
    #x = [0] * ndim
    max_norm = sys.maxsize
    #smallest_norm_index = -1
    for i in range(ndim):
        qnorm = quadratic_norm(MS_ZZ_c(Q[i]).transpose(), Q)
        if(qnorm < max_norm):
            max_norm = qnorm
            #smallest_norm_index = i
            #x = Q[i]
    return max_norm

def restart_program():
    """Restarts the current program."""
    print("Retrying...")
    python = sys.executable  # Path to the Python interpreter
    os.execl(python, python, *sys.argv)  # Replace the current process with a new one

# END AUX FUNCTIONS

def main():
    print("\n\n\nBEGIN KEYGEN::::\n")
    def discrete_n_dim_sampler():
        sampler = DiscreteGaussianDistributionIntegerSampler(s)
        initial_basis = [0] * ndim


        for i in range(ndim):
            ivec = [0] * ndim
            for k in range(ndim):
                ivec[k] = sampler()
            initial_basis[i] = ivec

        initial_basis = MS_ZZ(initial_basis)
        
        return initial_basis


    #sampledmatrix = discrete_n_dim_sampler(int(ndim))

    #print("sampledmatrix\n", sampledmatrix)

    # Define a basis matrix for the lattice
    #basis_matrix = IntegerMatrix.from_matrix(sampledmatrix)

    # Do a LLL reduction - This is the lattice basis
    #reduced = LLL.reduction(basis_matrix)
    reduced = identity_matrix(ndim)

    
    
    #U = discrete_n_dim_sampler()
    ##output sample has to be invertible!
    #while(U.determinant() != 1 or U.determinant() != -1):
    #    U = discrete_n_dim_sampler()

    # Obtain a random unimodular matrix in GL
    matrix_space = sage.matrix.matrix_space.MatrixSpace(ZZ, ndim)
    U = random_unimodular_matrix(matrix_space).LLL()


    print("S=I_n", reduced)
    print()
    print("PRIV/SECRET KEY, U\n", U)
    print("SECRET KEY DETERMINANT: ", U.determinant())

    # This is the quadratic form S from the KEM
    #S = Matrix(ZZ, reduced.transpose() * reduced)
    #S = S.LLL()
    S = reduced

    print("QUADRATIC FORM S: \n", S)

    # This is the P from the KEM
    # If S=I_n, then P =U^t * I_n * U = U^t * U
    P = U.transpose() * U
    #P = P.LLL()

    print("PUBLIC KEY/QUADRATIC FORM P: \n", P)

    #secret_key = U
    #public_key = P



    ### END KEYGEN

    ### BEGIN ENCAPS
    print("\n\n\nBEGIN ENCAPS::::\n")



    ## sampling e vector
    ## i think this gets stuck if the pkey is bad(need to define what bad means... see HAWK spec)
    ## s = 500 works well for ndim=2; 
    small_s = lambda_n_of_quadratic(P)
    def discrete_n_dim_vector_sampler(q, rho, pkey):
        sampler = DiscreteGaussianDistributionLatticeSampler(S, (q  * rho) / sqrt(ndim))#, c=(q/2,q/2,q/2))#5 * math.pow(10, ndim))#
        sample = sampler()
        e_vec = [0] * ndim

        if(type(sample) == type(None)):
            return restart_program()
        if(sample.is_zero()):
            return restart_program()
            #sample = sampler()
        print("sample ", sample)

        #copy to itself so the vector becomes mutable
        sample = copy.copy(sample)

        for k in range(ndim):
            #print("RR(sample[k])::: ", RR(sample[k]), " q::: ", q)
            
            #we want k to have positive values (this helps me with completeness)
            #while(sample[k] < 0):
            #    sample[k] = sample[k] + q
            
            e_vec[k] = sample[k] / q
        
        print("e_vec ", e_vec)
        #check if sample has norm less than rho, if not, repeat
        e_norm = quadratic_norm(MS_QQ_c(e_vec).transpose(), S)
        #if(e_norm > rho):
        #    print("sample for e_vec has norm > rho")
        #    restart_program()

        return e_vec

    e = discrete_n_dim_vector_sampler(q, rho, P) 

    print("e vector: ", e)

    # c vector definition: vector e over the discretized torus Tq
    c = copy.copy(e)
    for i in range(ndim):#if(c[i] >= 0):
        c[i] = c[i] - floor(c[i])
        #else: c[i] = c[i] - ceil(c[i])

    print("c vector ::: ", c)

    #CHECK if c = e. should the program abort in this case?
    #if(c == e):
    #    print("Vectors c and e are equal, restart.")
    #    restart_program()


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

    encaps_symmetric_k = extractor(e)

    print("encaps_symmetric_k: ", encaps_symmetric_k)

    #old: from universal hash function model
    #encaps_ciphertext = [c, Z]

    encaps_ciphertext = [c]

    output_encaps = [encaps_ciphertext, encaps_symmetric_k]

    #print("((c, Z), k) :::: ", output_encaps)
    print("((c), k) :::: ", output_encaps)


    ## END ENCAPS

    ## BEGIN DECAPS



#OLD CODE FOR DECODER
    #returns a vector from S such that ||y - Uc||_S <= rho
#    def decode(S, U, c):
#        Uc = U * c
#        sampler = DiscreteGaussianDistributionLatticeSampler(S, (q  * rho) / sqrt(ndim))
#        y = sampler()
#        if(y.is_zero()):
#            restart_program()
#        print("y sample at decode: ", y)
#        quad_norm = quadratic_norm(matrix_subtract(MS_ZZ_c(y).transpose(), Uc), S)
#        print("quad_norm of y-(U*c) :  ", quad_norm)
#
#        #first_min = first_minimum_S(S)
#        #print("first minimum of S: ", first_min)
#
#        if(quad_norm < 1):
#            print("GOOD SAMPLE, SAMPLE: ", y)
#            return y
#
#        restart_program()


    #returns a vector from S such that ||y - Uc||_S <= rho
    #decoding in Zn is equivalent to rounding to the nearest integer
    def decode(S, U, c):
        
        Uc = U * c
        print("Uc:::::: ", Uc)

        y = Uc.apply_map(lambda x: round(x)) #if(x>=0):ceil(x) else:floor(x))

        #if(y.is_zero()):
        #    restart_program()
        print("y at decode: ", y)
        quad_norm = quadratic_norm(matrix_subtract(y, Uc), S)
        print("quad_norm of y-(U*c) :  ", quad_norm)

        #first_min = first_minimum_S(S)
        #print("first minimum of S: ", first_min)

        if(quad_norm <= rho):
            print("GOOD SAMPLE, SAMPLE: ", y)
            return y

        restart_program()


    
    def begin_decaps():
        print("\n\n\nBEGIN DECAPS::::\n")
        # inputs: secret_key = U; encaps_ciphertext = (c, Z)
        print("inputs: secret_key = U =\n", U, "\n; encaps_ciphertext = (c := (c, Z) : ", encaps_ciphertext, "\n")
        
        #Step 1: obtain y from S, Uc
        y = decode(S, MS_QQ(U), MS_QQ_c(c).transpose())

        #Step 2: compute k

        #U^-1 * y is the integer part of vector "e" from encaps
        U_invert_times_y = U.inverse() * y # MUST BE U times y and not U^-1 as in the scheme!?
        print("U_invert_times_y::::: ", U_invert_times_y)

        c_minus_Uy = matrix_subtract(MS_QQ_c(c), U_invert_times_y.transpose())


        #convert c_minus_Uy to list for extractor
        c_minus_Uy_list = [0] * ndim
        for i in range(ndim):
            #since we're working in the Torus, this next while is needed
            #while(c_minus_Uy[0, i] < 0):
            #    c_minus_Uy[0, i] = c_minus_Uy[0, i] + 1
            c_minus_Uy_list[i] = c_minus_Uy[0, i]
        print("c_minus_Uy_list ::::: ", c_minus_Uy_list)


        

        #decaps_symmetric_k = extractor(c_minus_Uy_list, Z)
        decaps_symmetric_k = extractor(c_minus_Uy_list)

        print("decaps_symmetric_k: ", decaps_symmetric_k)

        print("is encaps_symmetric_k equal to decaps_symmetric_k?: ", decaps_symmetric_k == encaps_symmetric_k)

    begin_decaps()

main()