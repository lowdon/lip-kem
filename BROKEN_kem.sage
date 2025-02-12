## THIS BROKEN VERSION HAS SOME ATTEMPTS AT CONSTRUCTING THE E AND C VECTORS IN ENCAPS - I'll use it for later (maybe?)

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

# END PARAMS



def restart_program():
    """Restarts the current program."""
    print("Restarting program...")
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

    print("QUADRATIC FORM S: \n", S)

    P = U.transpose() * S * U

    print("PUBLIC KEY/QUADRATIC FORM P: \n", P)

    secret_key = U
    public_key = P



    ### END KEYGEN

    ### BEGIN ENCAPS
    print("\n\n\nBEGIN ENCAPS::::\n")

    ro = 0.999999999999999 #0.99 #has to be less than 1 because we're working in Zn??? not sure (shortest vector has norm 1), no clue if this is a good value
    q = ceil((s * ndim / ro) * sqrt(log(2*ndim + 4)/pi)) #from Ducas' KEM 

    print("ro ", ro)
    print("q ", q)  

    ## sampling e vector
    ## i think this gets stuck if the pkey is bad(need to define what bad means... see HAWK spec)
    ## s = 500 works well for ndim=2; 
    def discrete_n_dim_vector_sampler(q, ro, pkey):
        sampler = DiscreteGaussianDistributionLatticeSampler(pkey, 5 * math.pow(10, ndim))#TODO: should be (q  * ro) / sqrt(ndim)) but is failing to sample with suggested value
        e_vec = [0] * ndim
        c_vec = [0] * ndim

        sample = sampler()
        if(type(sample) == type(None)):
            return restart_program()
        if(sample.is_zero()):
            return restart_program()
            #sample = sampler()
        print("sample ", sample)

        
        for k in range(ndim):
            print("RR(sample[k])::: ", RR(sample[k]), " q::: ", q)
            e_vec[k] = RR(sample[k]) / q
            print("e_vec[k] :::: ", e_vec[k])
            print("e_vec ", e_vec)
            print("MODULAR ARITH ::::" , (RR(sample[k]) % q))
            c_vec[k] = (RR(sample[k]) % q) / q
            

        print("sample e_vec / q = ", e_vec)
        print("sample c_vec / q = ", c_vec)

        return e_vec, c_vec

    e, c = discrete_n_dim_vector_sampler(q, ro, public_key) 

    print("e vector: ", e)

    # c vector definition: fractional parts of vector e
    #c = e
    #for i in range(len(e)):
    #    c[i] = frac(c[i])

    # TODO: ADD LOOP FOR CASE WHERE c = e, THIS CANNOT BE THE CASE EVER!!! while(c == e): resample e, reconstruct c


    print("c vector:\n", c)



    # define random seed z



    Z = random.randint(0, 2**ndim)

    # vec_e has dim: ndim, rand_z is a secure random binary sequence
    def extractor(vec_e, rand_z):
        hash_obj = hashlib.sha512()
        for i in range(len(vec_e)):
            rand_z = rand_z + float(vec_e[i]) # doubt: i decided summing the values, not sure there is a better alternative for a deterministic way of computing this
        hash_obj.update(str(rand_z).encode('utf-8'))
        return hash_obj.hexdigest()

    encaps_symmetric_k = extractor(e, Z)

    print("encaps_symmetric_k: ", encaps_symmetric_k)

    encaps_ciphertext = [c, Z]

    output_encaps = [encaps_ciphertext, encaps_symmetric_k]

    print(output_encaps)


    ## END ENCAPS

    ## BEGIN DECAPS

    def matrix_subtract(m1, m2):
        print(m1, " M1!")
        print(m2, " M2!")
        #m1 = list(m1)
        #m2 = list(m2)
        #res = [0] * ndim
        res = m1 - m2
        #for i in range(ndim):
        #    print(i)
        #    res[i] = float(m1[i] - m2[i])
        # Round each entry to the specified decimal places
        #res = res.apply_map(lambda x: round(x), 10)
        print(res, " RES!")
        return res

    
    def begin_decaps():
        print("\n\n\nBEGIN DECAPS::::\n")
        # inputs: secret_key = U; encaps_ciphertext = (c, Z)

        # we need to define vector y through decoding
        # sample from S a lattice vector: must have norm less than ro, must not be zero vector
        # note (7 feb): i noticed that when y is a zero vector, the key always successfully decodes with U*c. does this mean we don't need to sample a y from S to decode?
        MS_RR_c = MatrixSpace(RR, 1, ndim)

        # Define the MatrixSpace over RR (Real Field with 53 bits of precision)
        MS_RR = MatrixSpace(RR, ndim, ndim)
        real_U = MS_RR(U)
        # TODO:temp/debug
        new_new_e = []

        def decode(S, Uc):
            sampler = DiscreteGaussianDistributionLatticeSampler(S, (q  * ro) / sqrt(ndim))
            new_e = sampler()
            print("dec_e_sample::: ", new_e)
            Ue = real_U * MS_RR_c(new_e).transpose()

            Ue_plus_c = Ue + MS_RR_c(c).transpose()
            print("Ue_plus_c :::::" , Ue_plus_c)

            y = Ue_plus_c
            # CORRECT y = matrix_subtract(Uc, Ue)

            #print("Uc minus Ue = ", y)

            new_new_e = new_e # TODO REMOVEEEE

            return y

            
            
            """
            y = MS_RR_c(list(sampler()))
            #y = MS_RR_c([0] * ndim)
            print(y, " :::: Y VECTOR IS TYPE: ", type(y))
            #print(type(Uc), " Uc type!")
            y_Uc_norm = vector(matrix_subtract(y, Uc.transpose())).norm()

            if(type(sample) == type(None)):
                print("SAMPLE IS NONE; FAIL")
                return begin_decaps() #"SAMPLE IS NONE; FAIL"
            #elif(y.is_zero()): # or y_Uc_norm > ro
                #print("SAMPLE Y IS ZERO; FAIL")
                #return begin_decaps()
                #print("y_Uc_norm ", y_Uc_norm)
                #y = MS_RR_c(list(sampler()))
                #y_Uc_norm = vector(matrix_subtract(y, Uc.transpose())).norm()
            elif(y_Uc_norm > ro):
                print("Y NORM TOO LARGE (larger than ro)")
                return begin_decaps() #"Y NORM TOO LARGE (larger than ro)"
            else:    
                print("SAMPLE Y IS CORRECT!")
            
            return y
            """


        Uc = real_U * MS_RR_c(c).transpose()

        y_vec = decode(S, Uc)

        Uy = real_U.inverse() * y_vec

        c_minus_Uy = matrix_subtract(MS_RR_c(c), Uy.transpose())
        
        print("calculated e vector: ", c_minus_Uy)

        c_minus_Uy_list = [0] * ndim
        for i in range(ndim):
            c_minus_Uy_list[i] = c_minus_Uy[0, i]
        #c_minus_Uy_list = [c_minus_Uy[0,0], c_minus_Uy[0,1]] # for ndim=2

        # CORRECT decaps_symmetric_k = extractor(c_minus_Uy_list, Z)
        print("vector c::: ", c)
        decaps_symmetric_k = extractor(c, Z)

        print("decaps_symmetric_k: ", decaps_symmetric_k)

        retry_decaps_symmetric_k = extractor(new_new_e, Z)

        print("retry_decaps_symmetric_k: ", retry_decaps_symmetric_k)

        print("is encaps_symmetric_k equal to decaps_symmetric_k?: ", decaps_symmetric_k == encaps_symmetric_k)

    begin_decaps()

main()