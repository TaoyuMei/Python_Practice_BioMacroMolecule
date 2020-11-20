import numpy as np
def match_or_not(b1,b2):
    if ((b1 == 'A' and b2 == 'U')
        or (b1 == 'U' and b2 == 'A')
        or (b1 == 'G' and b2 == 'C')
        or (b1 == 'C' and b2 == 'G')
        or (b1 == 'G' and b2 == 'U')
        or (b1 == 'U' and b2 == 'G')
    ):
       return 1
    else:
        return 0

def calculate_bifurcation(i, j, value_matrix):
    sum = []
    sum_k = {}
    k = i + 1
    if(k >= j-1):
        sum = [0]
        sum_k[0] = k
    while (k < j - 1):
        a_sum = value_matrix[i, k] + value_matrix[k + 1, j]
        sum.append(a_sum)
        sum_k[a_sum] = k
        k += 1
    max_sum = np.max(np.array(sum))
    the_k = sum_k[max_sum]
    return (max_sum, the_k)

rna_sequence = "AAACUUUCCCAGGG"
matrix_scale = len(rna_sequence)
value_matrix = np.zeros(shape=(matrix_scale,matrix_scale))

# fill in the matrix
j = 2   # matrix column
i = 0   # matrix raw
n = 2
while(n < matrix_scale):
    while(j < matrix_scale):
        s_i_j = match_or_not(rna_sequence[i], rna_sequence[j])
        diagonal = value_matrix[i+1,j-1] + s_i_j
        horizontal = value_matrix[i,j-1]
        vertical = value_matrix[i+1,j]
        sum_key_tuple = calculate_bifurcation(i, j, value_matrix)
        bifurcation = sum_key_tuple[0]
        value_matrix[i,j] = np.max(np.array([diagonal, horizontal, vertical, bifurcation]))
        i += 1
        j += 1
    n += 1
    j = n
    i = 0

print(value_matrix)
max_num_base_pair = value_matrix[0, matrix_scale-1]

# backtrack
matched_index_pair = {}
i = 0
j = matrix_scale - 1
def backtrack(i, j, value_matrix):

    while(i < j):
        sum_key_tuple = calculate_bifurcation(i, j, value_matrix)
        s_i_j = match_or_not(rna_sequence[i], rna_sequence[j])
        if(value_matrix[i+1,j] == value_matrix[i,j]):
            (i, j) = (i+1,j)
        elif(value_matrix[i,j-1] == value_matrix[i,j]):
            (i, j) = (i,j-1)
        elif(value_matrix[i,j] == value_matrix[i+1,j-1] + s_i_j):
            matched_index_pair[i] = j
            (i, j) = (i+1,j-1)
        elif(value_matrix[i,j] == sum_key_tuple[0]):
            k = sum_key_tuple[1]
            return (backtrack(i, k, value_matrix), backtrack(k+1, j, value_matrix))

backtrack(i, j, value_matrix)
parenthesis = ["."]
parenthesis *= matrix_scale
for i in range(matrix_scale):
    if(matched_index_pair.get(i)):
        parenthesis[i] = "("
        parenthesis[matched_index_pair[i]] = ")"
parenthesis = ''.join(parenthesis)

print(rna_sequence)
print(parenthesis)
