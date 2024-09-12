import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
from sklearn.cluster import KMeans

style.use('dark_background')

x = np.array([[1,2],[1.5,1.8],[5,8],[8,8],[1,0.6],[9,11],[7,10]])

# plt.scatter(x[:,0], x[:,1], s=50, color='w')
# plt.show()

clf = KMeans(n_clusters=2)
clf.fit(x)

centroids = clf.cluster_centers_
labels = clf.labels_
# labels -> similar to yi in SVM [1, 0, -1]

colors = ['c.', 'w.', 'o.', 'r.', 'g.']

for i in range(len(x)):
    plt.plot(x[i][0], x[i][1], colors[labels[i]], markersize = 15)

plt.scatter(centroids[:,0], centroids[:,1], color = 'b', marker = '*', s=75, linewidths=5)
plt.show()








