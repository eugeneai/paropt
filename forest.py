from paropt import *


class ForestModel(object):
    """
    """

    def __init__(self, ):
        """
        """




def model():
    m=ForestModel()
    ip=ParOptProcess(m)
    I, X, U, Psi, it, _ = ip.optimize(m.t, eps=0.001, iters=200)
    print ("Result:", I, _)



if __name__=="__main__":
    model()
