import copy
from pprint import pprint

import numpy as np


def global_alignment(sequence1, sequence2):
    # Initialize matrices in the size of the input arrays plus 1 on each axis.
    # We will keep a value matrix and a origin matrix, to store the cell the value came from.
    val_row = [0 for i in range(len(sequence2) + 1)]
    dir_row = [(0, 0) for i in range(len(sequence2) + 1)]

    val_matrix = [copy.deepcopy(val_row) for i in range(len(sequence1) + 1)]
    dir_matrix = [copy.deepcopy(dir_row) for i in range(len(sequence1) + 1)]

    # Initialize first row
    for i in range(1, len(val_matrix[0])):
        val_matrix[0][i] = -i
        dir_matrix[0][i] = [0, i - 1]

    # Initialize the first column
    for i in range(1, len(val_matrix)):
        val_matrix[i][0] = -i
        dir_matrix[i][0] = [i - 1, 0]

    for i in range(1, len(val_matrix)):
        for j in range(1, len(val_matrix[i])):
            val_matrix, dir_matrix = get_value(sequence1, sequence2, val_matrix, dir_matrix, i, j)

    pprint(val_matrix)
    pprint(dir_matrix)

    # Run backtrace algorithm and get an optimal path
    path = backtrace(val_matrix, dir_matrix)
    print("Optimal alignment path: " + str(path))

    # Get the alignment from the path
    get_alignment(sequence1, sequence2, path)


'''
Sets the optimal value of a cell and the cell it originated from.
'''
def get_value(
        sequence1: str, sequence2: str,
        val_matrix: list, dir_matrix: list,
        row: int, column: int):

    # Calculate match/mismatch:
    match_val = -1
    if sequence1[row - 1] == sequence2[column - 1]:
        match_val = 1

    # Get the values of the three cells that this cell can get values from
    diag_dir = (row - 1, column - 1)
    left_dir = (row, column - 1)
    up_dir = (row - 1, column)

    # Calculate values from each direction
    up_val = val_matrix[up_dir[0]][up_dir[1]] - 1
    diag_val = val_matrix[diag_dir[0]][diag_dir[1]] + match_val
    left_val = val_matrix[left_dir[0]][left_dir[1]] - 1

    # Find the highest value, very naively:
    if diag_val >= up_val and diag_val >= left_val:
        val = diag_val
        origin = diag_dir
    elif up_val <= left_val:
        val = left_val
        origin = left_dir
    else:
        val = up_val
        origin = up_dir

    val_matrix[row][column] = val
    dir_matrix[row][column] = origin

    #print("Val: " + str(val) + " Origin: " + str(origin))
    #print("Sequence 1: " + sequence1[row - 1] + " Sequence two: " + sequence2[column - 1])
    return val_matrix, dir_matrix


'''
Backtraces an optimal alignment sequence.
'''
def backtrace(val_matrix: list, dir_matrix: list):

    # Get the starting cell
    cell = (len(val_matrix) - 1, len(val_matrix[0]) - 1)

    # Get and print the alignment score
    print("Alignment score: " + str(val_matrix[cell[0]][cell[1]]))

    path = [cell]
    while cell[0] != 0 and cell[1] != 0:
        cell = dir_matrix[cell[0]][cell[1]]
        path.append(cell)

    # Reverse the path
    path.reverse()
    return path


'''
Get our alignment from the provided path and sequences.
'''
def get_alignment(sequence1: str, sequence2: str, path: list):

    # Initialize empty alignments
    alignment1 = ""
    alignment2 = ""

    for i in range(1, len(path)):
        seq1_index = path[i][0] - 1
        seq2_index = path[i][1] - 1

        if path[i][0] == path[i - 1][0] + 1 and path[i][1] == path[i - 1][1] + 1:
            alignment1 += sequence1[seq1_index]
            alignment2 += sequence2[seq2_index]

        elif path[i][0] <= path[i - 1][0]:
            alignment1 += "_"
            alignment2 += sequence2[seq2_index]

        else:
            alignment1 += sequence1[seq1_index]
            alignment2 += "_"

    print("Alignment 1: " + alignment1)
    print("Alignment 2: " + alignment2)



'''
Creates a random sequence for testing. Not my code.
Credit: Professor Haynes Heaton, Auburn University
'''
def random_sequence(n):
    return("".join(np.random.choice(["A","C","G","T"], n)))

'''
Mutates a sequence for testing. Not my code. 
Credit: Professor Haynes Heaton, Auburn University
'''
def mutate(s, snp_rate, indel_rate):
    x = [c for c in s]
    i = 0
    while i < len(x):
        if np.random.random() < snp_rate:
            x[i] = random_sequence(1)
        if np.random.random() < indel_rate:
            length = np.random.geometric(0.5)
            if np.random.random() < 0.5: # insertion
                x[i] = x[i] + random_sequence(length)
            else:
                for j in range(i,i+length):
                    if j < len(x):
                        x[j] = ""
                    i += 1
        i += 1
    return("".join(x))


# Create a sequence and mutate it to imitate evolution
s1 = random_sequence(100)
s2 = mutate(s1, 0.1, 0.1)


# Check alignment
global_alignment(s1, s2)
