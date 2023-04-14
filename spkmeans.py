import mykmeanssp
import sys
import numpy as np
np.random.seed(0)

def spk(datapoints):
    pass

def print_spk_result(result):
    pass

def print_jacobi_result(result):
    pass

def print_matrix(matrix):
    for row in matrix:
        print(','.join('{:0.4f}'.format(i) for i in row))

ALGORITHMS = {
    "spk" : spk,
    "wam": mykmeanssp.wam,
    "ddg": mykmeanssp.ddg,
    "gl": mykmeanssp.gl,
    "jacobi": mykmeanssp.jacobi
}

OUTPUT_RESULT = {
    "spk": print_spk_result,
    "jacobi": print_jacobi_result,
    "wam": print_matrix,
    "ddg": print_matrix,
    "gl": print_matrix
}


def main():
    args = sys.argv
    goal = args[-2]
    file_name = args[-1]
    datapoints = np.loadtxt(file_name, dtype=np.double, delimiter=",")
    result = ALGORITHMS[goal](datapoints)
    OUTPUT_RESULT[goal](result)


if __name__ == "__main__":
    main()