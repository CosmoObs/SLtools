import numpy as np


def translation_and_rotation(x, y, x_0, y_0, theta_0, angle='radian'):
    """
     Changes the coordinate system by a translation followed by a rotation.

    Input:
     - x         <list> : x coordinates to be transformed
     - y         <list> : y  coordinates to be transformed
     - x_0      <float> : x coordinate of the new center of reference
     - y_0      <float> : y coordinate of the new center of reference
     - theta_0  <float> : angle of the (counterclockwise) rotation
     - angle      <str> : 'radian' if theta_0 is given in radians and 'degree' if is given in degrees

    Output:
     - x_transf <list> : transformed x coordinates
     - y_transf <list> : transformed x coordinates

    """

    if angle=='degree':
        theta_0 = np.pi / 180. * theta_0

    x = np.array(x)
    y = np.array(y)

    x_transf =  np.cos(theta_0) * (x - x_0) + np.sin(theta_0) * (y - y_0)
    y_transf = -np.sin(theta_0) * (x - x_0) + np.cos(theta_0) * (y - y_0)

    return x_transf, y_transf


