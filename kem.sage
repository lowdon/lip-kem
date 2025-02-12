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
    q = ceil((s * ndim / ro) * sqrt(log(2*ndim + 4)/pi)) #from Ducas' KEM TODO: delete q * 100000 from sampler line, just for DEBUG

    print("ro ", ro)
    print("q ", q)  

    ## sampling e vector
    ## i think this gets stuck if the pkey is bad(need to define what bad means... see HAWK spec)
    ## s = 500 works well for ndim=2; 
    def discrete_n_dim_vector_sampler(q, ro, pkey):
        sampler = DiscreteGaussianDistributionLatticeSampler(pkey, 5 * math.pow(10, ndim))#(q  * ro) / sqrt(ndim))
        e_vec = [0] * ndim

        sample = sampler()
        if(type(sample) == type(None)):
            return restart_program()
        if(sample.is_zero()):
            return restart_program()
            #sample = sampler()
        print("sample ", sample)

        for k in range(ndim):
            e_vec[k] = RR(sample[k]) / q

        return e_vec

    e = discrete_n_dim_vector_sampler(q, ro, public_key) 

    print("e vector: ", e)

    # c vector definition: fractional parts of vector e
    c = e
    for i in range(len(e)):
        c[i] = frac(c[i])

    # TODO: ADD LOOP FOR CASE WHERE c = e, THIS CANNOT BE THE CASE EVER!!! while(c == e): resample e, reconstruct c


    print("c vector:\n", c)



    # define random seed z



    Z = random.randint(0, s)

    # vec_e has dim: ndim, rand_z is a secure random binary sequence
    def extractor(vec_e, rand_z):
        hash_obj = hashlib.sha512()
        for i in range(len(vec_e)):
            rand_z = rand_z + float(vec_e[i]) # doubt: i decided summing the values, not sure there is a better alternative for a deterministic way of computing this
        hash_obj.update(str(rand_z).encode('utf-8'))
        return hash_obj.hexdigest()

    encaps_symmetric_k = extractor(e, Z)

    print("encaps_symmetric_k: ", encaps_symmetric_k)

    c_times_U = MS_RR_c(c) * MS_RR(U)   # 1x2 times 2x2

    encaps_ciphertext = [c_times_U, Z]

    output_encaps = [encaps_ciphertext, encaps_symmetric_k]

    print("((c_times_U, Z), k) :::: ", output_encaps)


    ## END ENCAPS

    ## BEGIN DECAPS

    def matrix_subtract(m1, m2):
        print(m1, " M1!")
        print(m2, " M2!")
        res = m1 - m2
        #for i in range(ndim - 1):
        #    print(i)
        #    res[i] = m1[i] - m2[i]
        print(res, " RES!")
        return res

    
    def begin_decaps():
        print("\n\n\nBEGIN DECAPS::::\n")
        # inputs: secret_key = U; encaps_ciphertext = (c, Z)
        print("inputs: secret_key = U =\n", U, "\n; encaps_ciphertext = (c := (c_times_U, Z) : ", encaps_ciphertext, "\n")

        # we need to define vector y through decoding
        # sample from S a lattice vector: must have norm less than ro, must not be zero vector
        # note (7 feb): i noticed that when y is a zero vector, the key always successfully decodes with U*c. does this mean we don't need to sample a y from S to decode?
        
        def decode(S, Uc):
            sampler = DiscreteGaussianDistributionLatticeSampler(S, (q  * ro) / sqrt(ndim))
            #y = MS_RR_c(list(sampler()))
            y = MS_RR_c([0] * ndim)
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
            #elif(y_Uc_norm > ro):
            #    print("Y NORM TOO LARGE (larger than ro)")
            #    return begin_decaps() #"Y NORM TOO LARGE (larger than ro)"
            else:    
                print("SAMPLE Y IS CORRECT!")
            
            return y



        
        real_U = MS_RR(U)

        #y_vec = MS_RR_c([0] * ndim) # de facto: decode(S, real_U * MS_RR_c(c).transpose())

        #Uy = real_U.inverse() * y_vec.transpose()

        #c_minus_Uy = MS_RR_c(c_times_U) # de facto: matrix_subtract(MS_RR_c(c), Uy.transpose())

        #c_minus_Uy_list = [0] * ndim
        #for i in range(ndim):
        #    c_minus_Uy_list[i] = c_minus_Uy[0, i]
        ##c_minus_Uy_list = [c_minus_Uy[0,0], c_minus_Uy[0,1]] # for ndim=2
        #print("c_minus_Uy_list ::::: ", c_minus_Uy_list)
        #decaps_symmetric_k = extractor(c_minus_Uy_list, Z)

        c_times_U_times_Uminus = c_times_U * MS_RR(U.inverse()) # 2x2 times 2x1
        print("c_times_U_times_Uminus ::::: ", c_times_U_times_Uminus)

        c_times_U_times_Uminus_list = [0] * ndim
        for i in range(ndim):
            c_times_U_times_Uminus_list[i] = c_times_U_times_Uminus[0, i]
        #c_minus_Uy_list = [c_minus_Uy[0,0], c_minus_Uy[0,1]] # for ndim=2
        print("c_times_U_times_Uminus_list ::::: ", c_times_U_times_Uminus_list)

        decaps_symmetric_k = extractor(c_times_U_times_Uminus_list, Z)

        print("decaps_symmetric_k: ", decaps_symmetric_k)

        print("is encaps_symmetric_k equal to decaps_symmetric_k?: ", decaps_symmetric_k == encaps_symmetric_k)

    begin_decaps()

main()