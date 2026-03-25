import numpy as np
import contextlib
import umap

@contextlib.contextmanager
def temp_seed(seed):
    state = np.random.get_state()
    np.random.seed(seed)

    try :
        yield
    finally :
        np.random.set_state(state)


def get_embedding(data, random_state):
    reducer = umap.UMAP(random_state=random_state, n_neighbors=15, n_components=2, min_dist=0.1)
    data_transformed = reducer.fit_transform(data)

    return data_transformed
