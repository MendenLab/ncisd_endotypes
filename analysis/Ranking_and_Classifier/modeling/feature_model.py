
import numpy as np
#from sklearn.linear_model import LinearRegression

class FeatureModel :
    def __init__( self, classifier ):
        self._classifier = classifier
        self.models = None

    def fit( self, X, y ): 
        X = X.to_numpy().reshape([-1,1])
        columns = y.columns
        self.models = []

        for col in columns :
            y_ = y[col].to_numpy().ravel()
            
            m = self._classifier()
            m.fit( X,y_ )

            self.models.append(m)

    def eval( self, X ):
        
        X = X.to_numpy().reshape([-1,1])

        out = []
        for m in self.models :
            r = m.predict(X)
            out.append( r )

        return np.array( out ).T


