import mykmeanssp
import sys
import numpy as np

np.random.seed(0)


def spk(datapoints):
    np.random.seed(0)
    args = sys.argv
    gl = np.array(mykmeanssp.gl(datapoints))
    eigenvalues, eigenvectors = mykmeanssp.jacobi(gl)
    eigenvalues = np.array(eigenvalues)
    eigenvectors = np.array(eigenvectors)
    sort_indices = np.argsort(eigenvalues)
    eigenvectors = np.array(eigenvectors)[:, sort_indices]  # sorts based on eigenvalues
    eigenvalues.sort()
    if len(args) == 4:
        cluster_count = int(args[1])
    else:
        delta = np.abs(eigenvalues[:-1] - eigenvalues[1:])
        cluster_count = np.argmax(delta[:len(eigenvalues) // 2]) + 1

    U = eigenvectors[:, :cluster_count]
    initial_centroids, initial_centroids_idx = kmeans_pp(cluster_count, U)

    res = mykmeanssp.spk(initial_centroids, U)
    return res, initial_centroids_idx


def kmeans_pp(k, datapoints):
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
    print(str(initial_centroids_idx)[1:-1].replace(" ", ""))
    print_matrix(final_centroids)


def print_jacobi_result(result):
    eigenvalues, eigenvectors = result
    print_matrix([eigenvalues])
    print_matrix(np.array(eigenvectors))


def print_matrix(matrix):
    for row in matrix:
        print(','.join('{:0.4f}'.format(i) for i in row))


ALGORITHMS = {
    "spk": spk,
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
