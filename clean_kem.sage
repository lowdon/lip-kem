from sage.all import Matrix
from fpylll import LLL, IntegerMatrix
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
from sage.matrix.constructor import random_unimodular_matrix
import hashlib, random
import os
import sys

# PARAMS
#dimension of lattices
ndim = 2

# Define the standard deviation parameter 's', must be larger than smoothing parameter
s = 2**(ndim) # Choose a positive real value # from corollary 3.3 HuckBennett's paper on rotations of Zn

# Define the MatrixSpace over RR (Real Field with 53 bits of precision)
MS_RR_c = MatrixSpace(RR, 1, ndim)
MS_RR = MatrixSpace(RR, ndim, ndim)
MS_ZZ = MatrixSpace(ZZ, ndim, ndim)
MS_ZZ_c = MatrixSpace(ZZ, 1, ndim)
MS_ZZ_c_transposed = MatrixSpace(ZZ, ndim, 1)

# END PARAMS



def restart_program():
    """Restarts the current program."""
    print("Retrying...")
    python = sys.executable  # Path to the Python interpreter
    os.execl(python, python, *sys.argv)  # Replace the current process with a new one



def main():
    print("\n\n\nBEGIN KEYGEN::::\n")
    def discrete_n_dim_sampler(n):
        sampler = DiscreteGaussianDistributionIntegerSampler(s)
        initial_basis = [0] * n


        for i in range(n):
            ivec = [0] * n
            for k in range(n):
                ivec[k] = sampler()
            initial_basis[i] = ivec
        
        return initial_basis


    sampledmatrix = discrete_n_dim_sampler(int(ndim))

    print("sampledmatrix\n", sampledmatrix)

    # Define a basis matrix for the lattice
    basis_matrix = IntegerMatrix.from_matrix(sampledmatrix)

    # Do a LLL reduction - This is the lattice basis
    reduced = LLL.reduction(basis_matrix)

    # Obtain a random unimodular matrix in GL
    matrix_space = sage.matrix.matrix_space.MatrixSpace(ZZ, ndim)
    U = random_unimodular_matrix(matrix_space)


    print("LLLreduced_sampledmatrix\n", reduced)
    print()
    print("PRIV/SECRET KEY, U\n", U)
    print("SECRET KEY DETERMINANT: ", U.determinant())

    # This is the quadratic form S from the KEM
    S = Matrix(ZZ, reduced.transpose() * reduced)
    S = S.LLL()

    print("QUADRATIC FORM S: \n", S)

    # This is the P from the KEM
    P = U.transpose() * S * U
    P = P.LLL()

    print("PUBLIC KEY/QUADRATIC FORM P: \n", P)

    secret_key = U
    public_key = P



    ### END KEYGEN

    ### BEGIN ENCAPS
    print("\n\n\nBEGIN ENCAPS::::\n")

    ro = 0.999999999999999 #0.99 #has to be less than 1 because we're working in Zn??? not sure (shortest vector has norm 1), no clue if this is a good value
    q = ceil((s * ndim / ro) * sqrt(log(2*ndim + 4)/pi)) #from Ducas' KEM TODO: delete q * 100000 from sampler line, just for DEBUG

    print("ro ", ro)
    print("q ", q)  

    ## sampling e vector
    ## i think this gets stuck if the pkey is bad(need to define what bad means... see HAWK spec)
    ## s = 500 works well for ndim=2; 
    def discrete_n_dim_vector_sampler(q, ro, pkey):
        sampler = DiscreteGaussianDistributionLatticeSampler(pkey, (q  * ro) / sqrt(ndim))#5 * math.pow(10, ndim))#
        e_vec = [0] * ndim

        sample = sampler()
        if(type(sample) == type(None)):
            return restart_program()
        if(sample.is_zero()):
            return restart_program()
            #sample = sampler()
        print("sample ", sample)

        for k in range(ndim):
            e_vec[k] = sample[k]

        return e_vec

    e = discrete_n_dim_vector_sampler(q, ro, public_key) 

    print("e vector: ", e)

    # c vector definition: vector e over the discretized torus Tq
    c = e
    for i in range(ndim):
        c[i] = Mod(c[i], q)

    # TODO: ADD LOOP FOR CASE WHERE c = e, THIS CANNOT BE THE CASE EVER!!! while(c == e): resample e, reconstruct c


    print("c vector ::: ", c)



    # define random seed z
    Z = random.randint(0, s)

    # vec_e has dim = ndim, rand_z needs to be a computationally secure random binary sequence
    def extractor(vec_e, rand_z):
        print("inputs: vec_e, rand_z ::: ", vec_e, rand_z)
        hash_obj = hashlib.sha512()
        for i in range(len(vec_e)):
            rand_z = rand_z + vec_e[i] # doubt: i decided summing the values, not sure there is a better alternative for a deterministic way of computing this
        hash_obj.update(str(rand_z).encode('utf-8'))
        return hash_obj.hexdigest()

    encaps_symmetric_k = extractor(e, Z)

    print("encaps_symmetric_k: ", encaps_symmetric_k)

    encaps_ciphertext = [c, Z]

    output_encaps = [encaps_ciphertext, encaps_symmetric_k]

    print("((c, Z), k) :::: ", output_encaps)


    ## END ENCAPS

    ## BEGIN DECAPS

    def matrix_subtract(m1, m2):
        print(m1, " M1!")
        print(m2, " M2!")
        res = m1 - m2
        print(res, " RES!")
        return res

    #Let x be a vector in Zn. The norm of x over the quadratic form Q is ||x||_Q = sqrt(x^t * Q * x)
    def quadratic_norm(x, Q):
        x = MS_ZZ_c_transposed(x)
        return sqrt((x.transpose() * Q * x)[0][0])

    #Calculates first minimum of a lattice
    def first_minimum_S(Q):
        #x = [0] * ndim
        max_norm = sys.maxsize
        smallest_norm_index = -1
        for i in range(ndim):
            qnorm = quadratic_norm(Q[i], Q)
            if(qnorm < max_norm):
                max_norm = qnorm
                smallest_norm_index = i
                #x = Q[i]
        return max_norm

    #returns a vector from S such that ||y - Uc||_S <= rho
    def decode(S, U, c):
        Uc = U * c
        sampler = DiscreteGaussianDistributionLatticeSampler(S, (q  * ro) / sqrt(ndim))
        y = sampler()
        if(y.is_zero()):
            restart_program()
        print("y sample at decode: ", y)
        quad_norm = quadratic_norm(matrix_subtract(MS_ZZ_c(y).transpose(), Uc), S)
        print("quad_norm of y-(U*c) :  ", quad_norm)

        #first_min = first_minimum_S(S)
        #print("first minimum of S: ", first_min)

        if(quad_norm < 1):
            print("GOOD SAMPLE, SAMPLE: ", y)
            return y

        restart_program()


    
    def begin_decaps():
        print("\n\n\nBEGIN DECAPS::::\n")
        # inputs: secret_key = U; encaps_ciphertext = (c, Z)
        print("inputs: secret_key = U =\n", U, "\n; encaps_ciphertext = (c := (c, Z) : ", encaps_ciphertext, "\n")
        
        #Step 1: obtain y from S, Uc
        y = decode(S, MS_ZZ(U), MS_ZZ_c(c).transpose())

        #Step 2: compute k
        U_invert_times_y = U.inverse() * y

        c_minus_Uy = matrix_subtract(MS_ZZ_c(c), MS_ZZ_c(U_invert_times_y))


        #convert c_minus_Uy to list for extractor
        #round values close to integer
        c_minus_Uy_list = [0] * ndim
        threshold = 0.000000001
        for i in range(ndim):
            if(c_minus_Uy[0, i] < threshold and c_minus_Uy[0, i] > - threshold):
                c_minus_Uy[0, i] = 0
            c_minus_Uy_list[i] = c_minus_Uy[0, i]
        print("c_minus_Uy_list ::::: ", c_minus_Uy_list)


        

        decaps_symmetric_k = extractor(c_minus_Uy_list, Z)

        print("decaps_symmetric_k: ", decaps_symmetric_k)

        print("is encaps_symmetric_k equal to decaps_symmetric_k?: ", decaps_symmetric_k == encaps_symmetric_k)

    begin_decaps()

main()