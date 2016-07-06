import new_fargo3d as f
import numpy as np
import matplotlib.pyplot as plt

a = f.SnapShot(1)
nx = a.rho.data.shape[1]
data = a.rho.data.reshape([nx,nx])

diff = np.ndarray([nx/2,nx/2])

one   = np.flipud(np.fliplr(data[0:nx/2,0:nx/2]))
two   = np.flipud(np.fliplr(data[nx/2:,0:nx/2]))
three = data[nx/2:,nx/2:]
four  = data[0:nx/2,nx/2:]

plt.imshow(one-three)
plt.colorbar()
plt.show()
plt.imshow(two-four)
plt.colorbar()
plt.show()
