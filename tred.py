def tred(t,x,vred):
#   Apply a reduction velocity to the time axis
#Usage: t_red = tred(t,x,vred)
#   t: time
#   x: source-receiver ranges
#   vred: reduction velocity to apply
    import numpy as np

    t_red = t - np.abs(x)/vred
    return t_red