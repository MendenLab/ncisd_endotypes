import numpy as np
from sklearn.model_selection import train_test_split


class partitions:
    def __init__(self, cfg, dataset, seed):
        self._dataset = dataset

        _, _, _, _, indices_train, indices_test = train_test_split(
            dataset.data, dataset.pattern, np.arange(dataset.data.shape[0]),
            test_size=1 - cfg.ratio, random_state=seed, stratify=dataset.pattern)

        self._train_inds = indices_train
        self._test_inds = indices_test

        self._model_status = True

    def set_train(self):
        self._model_status = 'train'

    def set_test(self):
        self._model_status = 'test'

    def set_application(self):
        self._model_status = 'application'

    @property
    def data(self):
        if self._model_status == 'train':
            return self._dataset.data[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.data[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.data

    def sample_names(self):
        if self._model_status == 'train':
            return self._dataset.sample_names[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.sample_names[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.sample_names

    def responders(self):
        if self._model_status == 'train':
            return self._dataset.responders[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.responders[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.responders

    def observations(self):
        if self._model_status == 'train':
            return self._dataset.obs.iloc[self._train_inds, :]
        elif self._model_status == 'test':
            return self._dataset.obs.iloc[self._test_inds, :]
        else:
            # Use all data during application
            return self._dataset.obs

    @property
    def diag(self):
        if self._model_status == 'train':
            return self._dataset.diag[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.diag[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.diag

    @property
    def pattern(self):
        if self._model_status == 'train':
            return self._dataset.pattern[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.pattern[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.pattern

    @property
    def color_pattern(self):
        if self._model_status == 'train':
            return self._dataset.color_patterns[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.color_patterns[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.color_patterns

    @property
    def color_diag(self):
        if self._model_status == 'train':
            return self._dataset.color_diags[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.color_diags[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.color_diags

    @property
    def color_responders(self):
        if self._model_status == 'train':
            return self._dataset.color_diags[self._train_inds]
        elif self._model_status == 'test':
            return self._dataset.color_diags[self._test_inds]
        else:
            # Use all data during application
            return self._dataset.color_diags
