import mykmeanssp
import sys
import numpy as np
np.random.seed(0)

def spk(datapoints):
    args = sys.argv
    if len(args) == 3:
        cluster_count = args[0]
    else:
        gl = np.array(mykmeanssp.gl(datapoints))
        eigenvalues = mykmeanssp.jacobi(gl)[0]
        eigenvalues.sort()
        delta = np.abs(eigenvalues[:-1] - eigenvalues[1:])
        cluster_count = np.argmax(delta)
    initial_centroids, initial_centroids_idx = kmeans_pp(cluster_count, datapoints)
    res = mykmeanssp.spk(initial_centroids, datapoints)
    return res, initial_centroids_idx

def kmeans_pp(k, datapoints):
    np.random.seed(0)
    centers = []
    centers_idx = []
    idx = np.random.randint(low=0, high=datapoints.shape[0], size=(1,))[0]
    centers.append(datapoints[idx])
    centers_idx.append(idx)

    while len(centers) < k:
        distances = []
        for i, datapoint in enumerate(datapoints):
            distances.append(distance(centers, datapoint))

        dist_sum = np.sum(distances)
        probabilities = np.array(distances) / dist_sum
        new_center_idx = np.random.choice(datapoints.shape[0], size=(1,), p=probabilities)[0]
        centers.append(datapoints[new_center_idx])
        centers_idx.append(new_center_idx)

    return centers, centers_idx


def distance(centers, datapoint):
    return np.min(np.sqrt([np.sum(np.power(center - datapoint, 2)) for center in centers]))

def print_spk_result(result):
    final_centroids, initial_centroids_idx = result
    print_matrix([initial_centroids_idx])
    print_matrix(final_centroids)

def print_jacobi_result(result):
    eigenvalues, eigenvectors = result
    print_matrix([eigenvalues])
    print_matrix(np.array(eigenvectors).T)

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