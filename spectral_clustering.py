import numpy as np

import itertools
from scipy.spatial.distance import squareform, pdist
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

from matplotlib import pyplot as plt
from matplotlib import cm

def fit_labels(Y, Y_predicted):
		labels = list(set(Y) | set(Y_predicted))
		label_count = len(labels)

		# Set the cost matrix
		C = np.zeros((label_count, label_count))

		for i, l1 in enumerate(labels):
			for j, l2 in enumerate(labels):
				C[i, j] = len(np.intersect1d(np.where(Y == l1), np.where(Y_predicted == l2)))


		# Run the hungarian algorithm
		_, assignment = linear_sum_assignment(-C)

		asn_dct = {x: y for x, y in zip(assignment, labels)}

		# Reassign the labels
		Y_fit = np.zeros(Y_predicted.shape)

		for label in set(Y_predicted):
			Y_fit[Y_predicted == label] = asn_dct[label]

		return Y_fit


def spectral_clustering(X, max_num, k=0, verbose=False, metric='eucledian'):

	# Dist matrix
	S = squareform(pdist(X,metric=metric))
	sigma = np.percentile(S, 5)
	if sigma == 0: sigma = min([i for i in list(itertools.chain.from_iterable(S)) if i != 0])
	# Affinity matrix
	A = np.exp(- (S ** 2) / (2 * sigma ** 2))
	
	np.fill_diagonal(A, 0)

	m, n = A.shape
	
	# Compute laplacian
	I = np.identity(m)
	
	D = I / (A.sum(axis=1) ** 0.5)
	#print A.sum(axis=1) ** 0.5
	L = I - D.dot(A).dot(D)
	
	# Solve spectral decompositon
	eigen_val, eigen_vec = eigh(L)
	best = None
	if verbose:
		#plt.figure()
		#plt.scatter(range(1, max_num), eigen_val[:max_num - 1])
		#plt.show()
		best = sillhouette(eigen_vec[:, :k],max_num)
	if not best:
		labels = KMeans(n_clusters=k, n_init=100, max_iter=1000).fit_predict(eigen_vec[:, :k])
	elif k == 0:
		raise Exception("Silhouette not used to determine number of clusters and number of clusters not inputted mnanually")
	else:
		labels = KMeans(n_clusters=best, n_init=100, max_iter=1000).fit_predict(eigen_vec[:, :best])
	
	return labels


def sillhouette(X, max_num):

	range_n_clusters = list(range(2, max_num + 1))
	silhouette_avgs = []
	for n_clusters in range_n_clusters:
		# Create a subplot with 1 row and 2 columns
		#fig, (ax1, ax2) = plt.subplots(1, 2)
		#fig.set_size_inches(18, 7)

		# The 1st subplot is the silhouette plot
		# The silhouette coefficient can range from -1, 1 but in this example all
		# lie within [-0.1, 1]
		#ax1.set_xlim([-0.1, 1])
		# The (n_clusters+1)*10 is for inserting blank space between silhouette
		# plots of individual clusters, to demarcate them clearly.
		#ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

		# Initialize the clusterer with n_clusters value and a random generator
		# seed of 10 for reproducibility.
		clusterer = KMeans(n_clusters=n_clusters, random_state=10)
		cluster_labels = clusterer.fit_predict(X)

		# The silhouette_score gives the average value for all the samples.
		# This gives a perspective into the density and separation of the formed
		# clusters
		silhouette_avg = silhouette_score(X, cluster_labels)
		#print("For n_clusters =", n_clusters, "The average silhouette_score is :", silhouette_avg)
		silhouette_avgs.append(silhouette_avg)
		# Compute the silhouette scores for each sample
		sample_silhouette_values = silhouette_samples(X, cluster_labels)

		y_lower = 10
		continue
		for i in range(n_clusters):
			# Aggregate the silhouette scores for samples belonging to
			# cluster i, and sort them
			ith_cluster_silhouette_values = \
				sample_silhouette_values[cluster_labels == i]

			ith_cluster_silhouette_values.sort()

			size_cluster_i = ith_cluster_silhouette_values.shape[0]
			y_upper = y_lower + size_cluster_i

			color = cm.nipy_spectral(float(i) / n_clusters)
			ax1.fill_betweenx(np.arange(y_lower, y_upper),
							  0, ith_cluster_silhouette_values,
							  facecolor=color, edgecolor=color, alpha=0.7)

			# Label the silhouette plots with their cluster numbers at the middle
			ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

			# Compute the new y_lower for next plot
			y_lower = y_upper + 10  # 10 for the 0 samples

		ax1.set_title("The silhouette plot for the various clusters.")
		ax1.set_xlabel("The silhouette coefficient values")
		ax1.set_ylabel("Cluster label")

		# The vertical line for average silhouette score of all the values
		ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

		ax1.set_yticks([])  # Clear the yaxis labels / ticks
		ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

		# 2nd Plot showing the actual clusters formed
		colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
		ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
					c=colors, edgecolor='k')

		# Labeling the clusters
		centers = clusterer.cluster_centers_
		# Draw white circles at cluster centers
		ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
					c="white", alpha=1, s=200, edgecolor='k')

		for i, c in enumerate(centers):
			ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
						s=50, edgecolor='k')

		ax2.set_title("The visualization of the clustered data.")
		ax2.set_xlabel("Feature space for the 1st feature")
		ax2.set_ylabel("Feature space for the 2nd feature")

		plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
					  "with n_clusters = %d" % n_clusters),
					 fontsize=14, fontweight='bold')
	best = range_n_clusters[silhouette_avgs.index(max(silhouette_avgs))]
	return best
	#plt.show()


# from sklearn.datasets.samples_generator import make_blobs
# X, Y = make_blobs(n_samples=400, centers=4, cluster_std=0.60)

# plt.figure()
# plt.scatter(X[:, 0], X[:, 1])
# plt.show()

# Y_predicted = fit_labels(Y, spectral_clustering(X=X, k=4, verbose=True))

# print (sum([1 for x, y in zip(Y, Y_predicted) if x != y]))

