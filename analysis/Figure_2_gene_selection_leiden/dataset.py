import numpy as np
import scanpy as sc
import copy

from collections import OrderedDict


class dataset:
    def __init__(self, cfg):
        self._cfg = cfg
        data = sc.read(cfg.file_path)

        # Read in filtered, raw counts
        # self._data = np.array(data.X).astype(np.float64)
        self._data = np.array(data.layers['counts']).astype(np.int64)

        self.genes = data.var_names
        self.gene_ensemblid = data.var['gene_name']

        # Checking for valid patterns
        self._valid = data.obs.Pattern.isin(self._cfg.patterns_to_keep)

        # Mapping the diag values to integers
        self._diag, self.diag_mapping = self.map_label_str_int(str_labels=data.obs.diag.values)
        self._diag = np.array(self._diag, dtype=np.int64)

        # Mapping the pattern values to integers
        self._pattern, self.pattern_mapping = self.map_label_str_int(str_labels=data.obs.Pattern.values)
        self._pattern = np.array(self._pattern, dtype=np.int64)
        self.pattern_mapping = {key:  self.pattern_mapping[key] for key in self._cfg.patterns_to_keep}

        # Colors
        self.color_patterns = None
        self.color_diags = None
        self.convert_label_color()

        # Responders info
        if 'Response week 12' in data.obs.columns:
            self._responders = np.array(list(data.obs['Response week 12'].values), dtype=str)
        else:
            self._responders = np.array(['unknown'] * data.shape[0], dtype=str)

        # Sample names
        self._sample_names = np.array(list(data.obs.index))

        # Create joint column of pattern and batchID for train test set split
        assert 'batchID' in data.obs.columns, "Please add the column 'batchID' to your observations"
        data.obs['Pattern_batchID'] = data.obs[['Pattern', 'batchID']].apply(lambda x: '_'.join(x), axis=1)
        self._pattern_batchid = np.array(list(data.obs['Pattern_batchID'].values))

        # read out obs
        self._obs = data.obs

    def _color_pattern(self):
        colors = OrderedDict()
        colors['1'] = 'darkorange'
        colors['2a'] = 'darkred'
        colors['2b'] = 'tomato'
        colors['3'] = 'royalblue'
        colors['4a'] = 'mediumvioletred'
        colors['4b'] = 'magenta'
        colors['5'] = 'darkviolet'
        colors['UD'] = 'grey'
        colors['non-lesional'] = 'lightgrey'

        return colors

    def _color_diag(self):
        colors = OrderedDict()

        # Pattern undefined:
        colors['cutaneous lymphoma'] = 'dimgrey'
        colors['cutaneous side effects of biologics'] = 'darkgray'
        colors['darier disease'] = 'silver'
        colors['keratosis lichenoides chronica'] = 'grey'
        colors['erythrodermia'] = 'slategrey'
        colors['parapsoriasis'] = 'lightgrey'
        colors['undefined'] = 'whitesmoke'

        # Pattern 1:
        colors['lichen planus'] = 'orange'
        colors['lupus erythematosus'] = 'goldenrod'
        colors['lichenoid drug reaction'] = 'gold'

        # Pattern 2a:
        colors['eczema'] = 'maroon'
        colors['prurigo simplex subacuta'] = 'firebrick'

        # Pattern 2b:
        colors['bullous pemphigoid'] = 'salmon'

        # Pattern 3:
        colors['psoriasis'] = 'blue'
        colors['pityriasis rubra pilaris'] = 'dodgerblue'

        # Pattern 4a (darkgreen):
        colors['morphea'] = 'forestgreen'
        colors['venous ulcer'] = 'darkolivegreen'
        colors['systemic sclerosis'] = 'seagreen'

        # Pattern 4b (lightgreen):
        colors['granuloma annulare'] = 'springgreen'
        colors['sarcoidosis'] = 'palegreen'

        # Pattern 5 (violet):
        colors['psoriasis pustulosa'] = 'darkviolet'
        colors['pyoderma gangrenosum'] = 'mediumorchid'

        # No Pattern:
        colors['non lesional'] = 'bisque'

        return colors

    def __len__( self ):
        return len(self._valid)

    @property
    def data( self ):
        return copy.deepcopy( self._data[self._valid] )

    @property
    def sample_names(self):
        return copy.deepcopy(self._sample_names[self._valid])

    @property
    def diag( self ):
        return copy.deepcopy( self._diag[self._valid] )

    @property
    def pattern( self ):
        return copy.deepcopy( self._pattern[self._valid] )

    @property
    def pattern_batchid(self):
        return copy.deepcopy(self._pattern_batchid[self._valid])

    @property
    def responders(self):
        return copy.deepcopy(self._responders[self._valid])

    @property
    def obs(self):
        return copy.deepcopy(self._obs.loc[self._valid, :])

    @property
    def diag_colors( self ):
        if self.color_diags is None:
            self._add_color()
        return copy.deepcopy( self.color_diags )

    @property
    def pattern_colors( self ):
        if self.color_patterns is None:
            self._add_color()
        return copy.deepcopy( self.color_patterns )

    def _add_color(self):
        # get Pattern colors
        dict_color_patterns = self._color_pattern()
        self.color_patterns = list(dict_color_patterns.values())

        # get diag colors
        dict_color_diag = self._color_diag()
        self.color_diags = list(dict_color_diag.values())

    def convert_label_color(self):
        # get Pattern colors
        if np.all([[s for s in self.pattern if isinstance(s, int)]]):
            str_pattern_labels = self.map_label_int_str(int_labels=self.pattern, mapping=self.pattern_mapping)
        else:
            str_pattern_labels = self.pattern
        dict_color_patterns = self._color_pattern()
        self.color_patterns = np.asarray([dict_color_patterns[label] for label in str_pattern_labels])

        # get diag colors
        if np.all([[s for s in self.diag if isinstance(s, int)]]):
            str_diag_labels = self.map_label_int_str(int_labels=self.diag, mapping=self.diag_mapping)
        else:
            str_diag_labels = self.diag
        self.color_diags = np.asarray([self._color_diag()[label] for label in str_diag_labels])

    def map_label_str_int(self, str_labels):
        unique_str_labels = np.unique(str_labels)
        # map string labels to integer
        unique_int_labels = np.arange(0, len(unique_str_labels), 1)
        dict_int_str = dict(zip(unique_str_labels, unique_int_labels))
        int_labels = [dict_int_str[label] for label in str_labels]

        return int_labels, dict_int_str

    def map_label_int_str(self, int_labels, mapping):
        # Swap key and values in mapping dict
        mapping = {value: key for key, value in mapping.items()}
        str_labels = [mapping[label] for label in int_labels]

        return str_labels
